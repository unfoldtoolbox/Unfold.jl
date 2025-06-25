
"""
    fit(type::UnfoldModel,d::Vector{Pair},tbl::AbstractDataFrame,data::Array)
    fit(type::UnfoldModel,f::FormulaTerm,tbl::AbstractDataFrame,data::Array{T,3},times)
    fit(type::UnfoldModel,f::FormulaTerm,tbl::AbstractDataFrame,data::Array{T,2},basisfunction::BasisFunction)

Generates Designmatrix & fits model, either mass-univariate (one model per epoched-timepoint) or time-expanded (modeling linear overlap).

## keyword arguments
- `fit::Bool` (default: `true`) - fit the model after constructing the designmatrix. Setting this to `false` is sometimes helpful if you only want to inspect the designmatrix.
- `contrasts::Dict`: (default: `Dict()`) contrast to be applied to formula. Example: `Dict(:my_condition=>EffectsCoding())`. More information here: https://juliastats.org/StatsModels.jl/stable/contrasts/
- `eventcolumn::Union{Symbol,String}` (default `:event`) - the column in `tbl` to differentiate the basisfunctions as defined in `d::Vector{Pair}`
- `solver::function`: (default: `solver_default`). The solver used for `y=Xb`, e.g. `(X,y;kwargs...) -> solver_default(X,y;kwargs...)`. There are faster & alternative solvers available, see `solver_predefined` for a list of options, see `solver benchmark` in the online documentation. To use the GPU, you can provide the data as a `CuArray` after `using CUDA`. Please change the solver to e.g. `solver_predef(X,y;solver=:qr)` as lsmr+cuda => crash typically. It's worth though, speed increases >100x possible
- `show_progress::Bool` (default `true`) - show progress via ProgressMeter - passed to `solver`
- `eventfields::Array: (optional, default `[:latency]`) Array of symbols, representing column names in `tbl`, which are passed to basisfunction event-wise. First field of array always defines eventonset in samples.

If a `Vector[Pairs]` is provided, it has to have one of the following structures:
For **deconvolution** analyses (use `Any=>(f,bf)` to match all rows of `tbl` in one basis functions). Assumes `data` is a continuous EEG stream, either a `Vector` or a `ch x time` `Matrix`
```julia
f1 = @formula(0~1+my_condition)
[
 :A=>(f1,firbasis((-0.1,1),128), # sfreq = 128Hz
 :B=>(f2,firbasis((-3,2),128)
]
```
for **mass-univariate** analyses without deconvolution. Assumes `data` to be cut into epochs already (see `Unfold.epoch`). Follows *eeglab* standard `ch x time x trials`:
```julia
timesvector = range(-0.1,3,step=1/100)
[
 :A=>(f1,timesvector),
 :B=>(f2,timesvector)
]
```

## Notes
- The `type` can be specified directly as well e.g. `fit(type::UnfoldLinearModel)` instead of relying on the automatic inference
- The data is reshaped if it is missing one dimension to have the first dimension then `1` "Channel".

## Examples
Mass Univariate Linear
```julia-repl
julia> data,evts = UnfoldSim.predef_eeg()
julia> data_e,times = Unfold.epoch(data=data,tbl=evts,τ=(-1.,1.9),sfreq=100) # cut the data into epochs. data_e is now ch x times x epoch

julia> f  = @formula 0~1+continuousA+continuousB
julia> model = fit(UnfoldModel,f,evts,data_e,times)
# or:
julia> model = fit(UnfoldModel,[Any=>(f,times)],evts,data_e)
```
Timexpanded Univariate Linear
```julia-repl
julia> basisfunction = firbasis(τ=(-1,1),sfreq=10)
julia> model = fit(UnfoldModel,f,evts,data,basisfunction)
# or
julia> model = fit(UnfoldModel,[Any=>(f,basisfunction],evts,data)
```

"""
function StatsModels.fit(
    UnfoldModelType::Type{T},
    f::FormulaTerm,
    tbl::AbstractDataFrame,
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
    tbl::AbstractDataFrame,
    data::AbstractArray{T};
    kwargs...,
) where {T}

    for k = 1:length(design)
        d_first = design[k]
        d_tuple = last(d_first)
        @assert typeof(first(d_tuple)) <: FormulaTerm "InputError in design(uf) - :key=>(FORMULA,basis/times), formula not found. Maybe formula wasn't at the first place?"
    end



    if UnfoldModelType == UnfoldModel
        @debug "autodetecting UnfoldModelType"
        UnfoldModelType = design_to_modeltype(design)
    end

    for k = 1:length(design)
        d_first = design[k]
        d_tuple = last(d_first)
        @assert (typeof(last(d_tuple)) <: AbstractVector) ⊻
                (SimpleTraits.istrait(Unfold.ContinuousTimeTrait{UnfoldModelType})) "InputError: Either a basis function was declared, but a UnfoldLinearModel was built, or a times-vector (and no basis function) was given, but a UnfoldLinearModelContinuousTime was asked for."
    end
    @debug "Check Data + Applying UnfoldModelType: $UnfoldModelType {$T}"
    data_r = check_data(UnfoldModelType, data)

    fit(UnfoldModelType{T}(design), tbl, data_r; kwargs...)
end


function StatsModels.fit(
    uf::UnfoldModel,
    tbl::AbstractDataFrame,
    data::AbstractArray;
    fit = true,
    kwargs...,
)
    @debug "adding desigmatrix!"
    designmatrix!(uf, tbl; kwargs...)
    if fit
        fit!(uf, data; kwargs...)
    end

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


    if isa(uf, UnfoldLinearModel)
        @assert length(times(uf)[1]) == size(data, length(size(data)) - 1) "Times Vector does not match second last dimension of input data - forgot to cut into epochs?"
    end

    X = modelmatrix(uf)


    if isa(uf, UnfoldLinearModel)

        d = designmatrix(uf)

        # mass univariate with multiple events fitted at the same time

        coefs = []
        for m = 1:length(X)
            # the main issue is, that the designmatrices are subsets of the event table - we have
            # to do the same for the data, but data and designmatrix dont know much about each other.
            # Thus we use parentindices() to get the original indices of the @view events[...] from desigmatrix.jl
            @debug typeof(X) typeof(events(d)[m])
            push!(coefs, solver(X[m], @view data[:, :, parentindices(events(d)[m])[1]]))
        end
        @debug [size(c.estimate) for c in coefs]
        uf.modelfit = LinearModelFit{T,3}(
            cat([c.estimate for c in coefs]..., dims = 3),
            [c.info for c in coefs],
            cat([c.standarderror for c in coefs]..., dims = 3),
        )
        return # we are done here

        #        elseif isa(d.events, SubDataFrame)
        # in case the user specified an event to subset (and not any) we have to use the view from now on
        #            data = @view data[:, :, parentindices(d.events)[1]]
        #        end
    end


    X, data = equalize_size(X, data)
    #    @debug typeof(uf.modelfit), typeof(T), typeof(X), typeof(data)
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

@traitfn check_data(
    uf::Type{UF},
    data::AbstractArray{T,3},
) where {T,UF<:UnfoldModel;ContinuousTimeTrait{UF}} = error(
    "A continuous-time model was request, but the data provided have 3 dimensions. Did you maybe provide epoched data instead of continuous data?",
)

function design_to_modeltype(design)
    #@debug design
    # autoDetect
    tmp = last(design[1])
    f = tmp[1] # formula
    t = tmp[2] # Vector or BasisFunction
    design_to_modeltype(f, t)
end

"""
 !!! Important:
    this is an ugly hack dating back to the time where UnfoldMixedModels was still an extension. We are overloading this function in UnfoldMixedModels.jl with a more specific type, to switch between MixedModels-Unfold types and not...
 """
design_to_modeltype(f, t::BasisFunction) = UnfoldLinearModelContinuousTime
design_to_modeltype(f, t::AbstractVector) = UnfoldLinearModel
