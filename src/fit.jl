
"""
fit(type::UnfoldModel,f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times)
fit(type::UnfoldModel,f::FormulaTerm,tbl::DataFrame,data::Array{T,2},basisfunction::BasisFunction)
fit(type::UnfoldModel,d::Dict,tbl::DataFrame,data::Array)

Generates Designmatrix & fits model, either mass-univariate (one model per epoched-timepoint) or time-expanded (modeling linear overlap).


# Examples
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
julia> model = fit(UnfoldModel,Dict(Any=>(f,basisfunction),evts,data_r)
```

"""
function StatsModels.fit(
    UnfoldModelType::Type{T},
    f::FormulaTerm,
    tbl::DataFrame,
    data::AbstractArray,
    basisOrTimes::Union{BasisFunction,AbstractArray};
    kwargs...,
) where {T<:Union{<:UnfoldModel}}
    # old input format, sometimes convenient.Convert to dict-based one
    fit(UnfoldModelType, Dict(Any => (f, basisOrTimes)), tbl, data; kwargs...)
end


function StatsModels.fit(
    UnfoldModelType::Type{UnfoldModel},
    design::Dict,
    tbl::DataFrame,
    data::AbstractArray;
    kwargs...,
    )
    detectedType = designToModeltype(design)

    fit(detectedType,design,tbl,data;kwargs...)
end


function StatsModels.fit(
    UnfoldModelType::Type{<:UnfoldModel},
    design::Dict,
    tbl::DataFrame,
    data::AbstractArray;
    kwargs...,
)
fit(UnfoldModelType(design),design,tbl,data;kwargs...)
end


function StatsModels.fit(
    uf::UnfoldModel,#Union{UnfoldLinearMixedModel,UnfoldLinearModel,UnfoldLinearMixedModelContinuousTime,UnfoldLinearModelContinuousTime},
    design::Dict,
    tbl::DataFrame,
    data::AbstractArray;
    kwargs...,
    )
    
    designmatrix!(uf, tbl; kwargs...)
    fit!(uf, data; kwargs...)

    return uf
end

function StatsModels.fit(
    UnfoldModelType::Type{T},
    X::DesignMatrix,
    data::AbstractArray;
    kwargs...,
) where {T<:Union{<:UnfoldModel}}
    if UnfoldModelType == UnfoldModel
        error(
            "Can't infer model automatically, specify with e.g. fit(UnfoldLinearModel...) instead of fit(UnfoldModel...)",
        )
    end
    uf = UnfoldModelType(Dict(), X)

    fit!(uf, data; kwargs...)

    return uf
end

isMixedModelFormula(f::FormulaTerm) = isMixedModelFormula(f.rhs)
isMixedModelFormula(f::Tuple) = any(isMixedModelFormula.(f))

isMixedModelFormula(f::InteractionTerm) = false
isMixedModelFormula(f::ConstantTerm) = false
isMixedModelFormula(f::Term) = false
#isMixedModelFormula(f::FunctionTerm) = false
function isMixedModelFormula(f::FunctionTerm) 
    try 
        isMixedModelFormula(f.f)
    catch
        isMixedModelFormula(f.forig) # StatsMoels  <0.7
    end
end
isMixedModelFormula(f::Function) = false  # catch all

isMixedModelFormula(f::typeof(|)) = true

function designToModeltype(design)
    # autoDetect
    tmp = collect(values(design))[1]
    f = tmp[1] # formula
    t = tmp[2] # Vector or BasisFunction

    isMixedModel = isMixedModelFormula(f)

    if typeof(t) <: BasisFunction
        if isMixedModel
            UnfoldModelType = UnfoldLinearMixedModelContinuousTime
        else
            UnfoldModelType = UnfoldLinearModelContinuousTime
        end
    else
        if isMixedModel
            UnfoldModelType = UnfoldLinearMixedModel
        else
            UnfoldModelType = UnfoldLinearModel
        end
    end
    return UnfoldModelType
end


# helper function for 1 channel data
function StatsModels.fit(
    ufmodel::T,
    design::Dict,
    tbl::DataFrame,
    data::AbstractVector,
    args...;
    kwargs...,
) where {T<:Union{<:UnfoldModel}}
    @debug("data array is size (X,), reshaping to (1,X)")
    data = reshape(data, 1, :)
    return fit(ufmodel, design, tbl, data, args...; kwargs...)
end

# helper to reshape a 
function StatsModels.fit(
    ufmodel::T,
    design::Dict,
    tbl::DataFrame,
    data::AbstractMatrix,
    args...;
    kwargs...)where {T<:Union{UnfoldLinearMixedModel,UnfoldLinearModel}}
    @debug("MassUnivariate data array is size (X,Y), reshaping to (1,X,Y)")
    data = reshape(data, 1, size(data)...)
    return fit(ufmodel, design, tbl, data, args...; kwargs...)
end




function StatsModels.fit!(
    uf::Union{UnfoldLinearModelContinuousTime,UnfoldLinearModel},
    data;
    solver = (x, y) -> solver_default(x, y),
    kwargs...,
)

    @assert ~isempty(designmatrix(uf))
    @assert typeof(first(values(design(uf)))[1]) <: FormulaTerm "InputError in design(uf) - :key=>(FORMULA,basis/times), formula not found. Maybe formula wasn't at the first place?"
    @assert (typeof(first(values(design(uf)))[2]) <: AbstractVector) ⊻ (typeof(uf) <: UnfoldLinearModelContinuousTime) "InputError: Either a basis function was declared, but a UnfoldLinearModel was built, or a times-vector (and no basis function) was given, but a UnfoldLinearModelContinuousTime was asked for."
    if isa(uf,UnfoldLinearModel)
        @assert length(first(values(design(uf)))[2]) == size(data,length(size(data))-1) "Times Vector does not match second last dimension of input data - forgot to epoch?"
    end
   
    X = modelmatrix(uf)

    @debug "UnfoldLinearModel(ContinuousTime), datasize: $(size(data))"
    
    if isa(uf,UnfoldLinearModel)
        d = designmatrix(uf)

        if isa(X,Vector)
        # mass univariate with multiple events fitted at the same time
        
        coefs = []
        for m = 1:length(X)
            # the main issue is, that the designmatrices are subsets of the event table - we have 
            # to do the same for the data, but data and designmatrix dont know much about each other.
            # Thus we use parentindices() to get the original indices of the @view events[...] from desigmatrix.jl
            push!(coefs,solver(X[m], @view data[:,:,parentindices(d.events[m])[1]]))
        end
        @debug @show [size(c.estimate) for c in coefs]
        uf.modelfit = LinearModelFit(
            cat([c.estimate for c in coefs]...,dims=3),
            [c.info for c in coefs],
            cat([c.standarderror for c in coefs]...,dims=3)
        )
        return # we are done here
   
        elseif isa(d.events,SubDataFrame)
            # in case the user specified an event to subset (and not any) we have to use the view from now on
            data = @view data[:,:,parentindices(d.events)[1]]
        end
    end


        # mass univariate, data = ch x times x epochs
        X, data = zeropad(X, data)

        uf.modelfit = solver(X, data)
        return

end


