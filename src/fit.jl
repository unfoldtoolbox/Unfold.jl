
"""
    fit(type::UnfoldModel,d::Vector{Pair},tbl::DataFrame,data::Array)
    fit(type::UnfoldModel,f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times)
    fit(type::UnfoldModel,f::FormulaTerm,tbl::DataFrame,data::Array{T,2},basisfunction::BasisFunction)

Generates Designmatrix & fits model, either mass-univariate (one model per epoched-timepoint) or time-expanded (modeling linear overlap).

- `eventcolumn` (Symbol/String, default :event) - the column in `tbl::DataFrame` to differentiate the basisfunctions as defined in `d::Vector{Pair}`
- `show_progress` (Bool, default true) - show Progress via ProgressMeter

If a `Vector[Pairs]` is provided, it has to have one of the following structures:
`[:A=>(f,basisfunction), :B=>(f2,bf2)]` - for deconvolutin analyses (use `Any=>(f,bf)` to match all rows of `tbl` in one basis functins)
`[:A=>(f,timesvector), :B=>(f2,timesvector)]` - for mass univariate analyses. If multiple rERPs are calculated at the same time, the timesvectors must be the same


## Notes
- The `type` can be specified directly as well e.g. `fit(type::UnfoldLinearModel)` instead of inferred
- The data is reshaped if it is missing one dimension to have the first dimesion then `1` "Channel".

## Examples
Mass Univariate Linear
```julia-repl
julia> data,evts = loadtestdata("testCase1")
julia> data_r = reshape(data,(1,:))
julia> data_e,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10) # cut the data into epochs. data_e is now ch x times x epoch

julia> f  = @formula 0~1+continuousA+continuousB # 1
julia> model = fit(UnfoldModel,f,evts,data_e,times)
```
Timexpanded Univariate Linear
```julia-repl
julia> basisfunction = firbasis(τ=(-1,1),sfreq=10,name="A")
julia> model = fit(UnfoldModel,[Any=>(f,basisfunction],evts,data_r)
```

"""
function StatsModels.fit(
    UnfoldModelType::Type{T},
    f::FormulaTerm,
    tbl::DataFrame,
    data::AbstractArray,
    basisOrTimes::Union{BasisFunction,AbstractArray};
    kwargs...,
) where {T<:UnfoldModel}
    # old input format, sometimes convenient.Convert to list-based one
    fit(UnfoldModelType, [Any => (f, basisOrTimes)], tbl, data; kwargs...)
end

# deprecated Dict call
function StatsModels.fit(uf::Type{UnfoldModel}, design::Dict, args...; kwargs...)
    #@debug keys(design) .=> Tuple.(values(design))
    @warn "using `Dict(:A=>(@formula,times/basisfunction))` is deprecated, please use `[:A=>(@formula,times/basisfunction)]` from now on"
    fit(uf, collect(pairs(design)), args...; kwargs...)
end


function StatsModels.fit(
    UnfoldModelType::Type{<:UnfoldModel},
    design::Vector,
    tbl::DataFrame,
    data::AbstractArray{T};
    kwargs...,
) where {T}
    if UnfoldModelType == UnfoldModel
        @debug "autodetecting UnfoldModelType"
        UnfoldModelType = design_to_modeltype(design)
    end
    @debug "Check Data + Applying UnfoldModelType: $UnfoldModelType {$T}"
    data_r = check_data(UnfoldModelType, data)

    fit(UnfoldModelType{T}(design), tbl, data_r; kwargs...)
end


function StatsModels.fit(uf::UnfoldModel, tbl::DataFrame, data::AbstractArray; kwargs...)
    @debug "adding desigmatrix!"
    designmatrix!(uf, tbl; kwargs...)
    fit!(uf, data; kwargs...)

    return uf
end


# special case where designmatrix and data are provided directly without a design
function StatsModels.fit(
    UnfoldModelType::Type{UF},
    X::Vector{<:AbstractDesignMatrix{T}},
    data::AbstractArray;
    kwargs...,
) where {UF<:UnfoldModel,T}
    if UnfoldModelType == UnfoldModel
        error(
            "Can't infer model automatically, specify with e.g. fit(UnfoldLinearModel...) instead of fit(UnfoldModel...)",
        )
    end
    uf = UnfoldModelType{T}([:empty => ()], X)


    fit!(uf, data; kwargs...)

    return uf
end



# main fitting function
function StatsModels.fit!(
    uf::UnfoldModel{T},
    data::AbstractArray{T};
    solver = (x, y) -> solver_default(x, y),
    kwargs...,
) where {T}
    @debug "fit!: $T, datasize: $(size(data))"

    @assert ~isempty(designmatrix(uf))
    d_first = design(uf)[1]
    d_tuple = last(d_first)
    @assert typeof(first(d_tuple)) <: FormulaTerm "InputError in design(uf) - :key=>(FORMULA,basis/times), formula not found. Maybe formula wasn't at the first place?"
    @assert (typeof(last(d_tuple)) <: AbstractVector) ⊻
            (typeof(uf) <: UnfoldLinearModelContinuousTime) "InputError: Either a basis function was declared, but a UnfoldLinearModel was built, or a times-vector (and no basis function) was given, but a UnfoldLinearModelContinuousTime was asked for."


    if isa(uf, UnfoldLinearModel)
        @assert length(times(uf)) == size(data, length(size(data)) - 1) "Times Vector does not match second last dimension of input data - forgot to epoch?"
    end

    X = modelmatrix(uf)


    if isa(uf, UnfoldLinearModel)
        d = designmatrix(uf)

        if isa(X, Vector)
            # mass univariate with multiple events fitted at the same time

            coefs = []
            for m = 1:length(X)
                # the main issue is, that the designmatrices are subsets of the event table - we have 
                # to do the same for the data, but data and designmatrix dont know much about each other.
                # Thus we use parentindices() to get the original indices of the @view events[...] from desigmatrix.jl
                push!(coefs, solver(X[m], @view data[:, :, parentindices(events(d)[m])[1]]))
            end
            @debug [size(c.estimate) for c in coefs]
            uf.modelfit = LinearModelFit{T,3}(
                cat([c.estimate for c in coefs]..., dims = 3),
                [c.info for c in coefs],
                cat([c.standarderror for c in coefs]..., dims = 3),
            )
            return # we are done here

        elseif isa(d.events, SubDataFrame)
            # in case the user specified an event to subset (and not any) we have to use the view from now on
            data = @view data[:, :, parentindices(d.events)[1]]
        end
    end


    # mass univariate, data = ch x times x epochs
    X, data = zeropad(X, data)
    @debug typeof(uf.modelfit), typeof(T), typeof(X), typeof(data)
    uf.modelfit = solver(X, data)
    return uf

end




@traitfn function check_data(
    uf::Type{UF},
    data::AbstractArray{T,2},
) where {T,UF<:UnfoldModel;!ContinuousTimeTrait{UF}}
    @debug(" !ContinuousTimeTrait: data array is size (X,Y), reshaping to (1,X,Y)")
    data_r = reshape(data, 1, size(data)...)
    return data_r
end
@traitfn check_data(
    uf::Type{UF},
    data::AbstractArray{T,2},
) where {T,UF<:UnfoldModel;ContinuousTimeTrait{UF}} = data

function check_data(uf::Type{<:UnfoldModel}, data::AbstractVector)
    @debug("data array is size (X,), reshaping to (1,X)")
    data = reshape(data, 1, :)
end

check_data(uf::Type{<:UnfoldModel}, data) = data


isa_lmm_formula(f::FormulaTerm) = isa_lmm_formula(f.rhs)
isa_lmm_formula(f::Tuple) = any(isa_lmm_formula.(f))

isa_lmm_formula(f::InteractionTerm) = false
isa_lmm_formula(f::ConstantTerm) = false
isa_lmm_formula(f::Term) = false
#isa_lmm_formula(f::FunctionTerm) = false
function isa_lmm_formula(f::FunctionTerm)
    try
        isa_lmm_formula(f.f)
    catch
        isa_lmm_formula(f.forig) # StatsMoels  <0.7
    end
end
isa_lmm_formula(f::Function) = false  # catch all

isa_lmm_formula(f::typeof(|)) = true


function design_to_modeltype(design)
    #@debug design
    # autoDetect
    tmp = last(design[1])
    f = tmp[1] # formula
    t = tmp[2] # Vector or BasisFunction

    isMixedModel = isa_lmm_formula(f)
    if isMixedModel
        ext = Base.get_extension(@__MODULE__, :UnfoldMixedModelsExt)
        if isnothing(ext)
            throw(
                "MixedModels not loaded. Please use `]add MixedModels` and `using MixedModels` to install it prior to using",
            )
        end
    end
    if typeof(t) <: BasisFunction
        if isMixedModel

            UnfoldModelType = ext.UnfoldLinearMixedModelContinuousTime
        else
            UnfoldModelType = UnfoldLinearModelContinuousTime
        end
    else
        if isMixedModel
            UnfoldModelType = ext.UnfoldLinearMixedModel
        else
            UnfoldModelType = UnfoldLinearModel
        end
    end
    return UnfoldModelType
end


