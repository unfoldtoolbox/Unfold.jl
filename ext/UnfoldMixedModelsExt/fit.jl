
"""
fit!(uf::UnfoldModel,data::Union{<:AbstractArray{T,2},<:AbstractArray{T,3}}) where {T<:Union{Missing, <:Number}}

Fit a DesignMatrix against a 2D/3D Array data along its last dimension
Data is typically interpreted as channel x time (with basisfunctions) or channel x time x epoch (for mass univariate)

Returns an UnfoldModel object

# Examples
```julia-repl
```

"""
function StatsModels.fit!(
    uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
    data::AbstractArray;
    kwargs...,
)

    #@assert length(first(values(design(uf)))[2])
    if uf isa UnfoldLinearMixedModel
        if ~isempty(Unfold.design(uf))
        @assert length(Unfold.times(Unfold.design(uf))) == size(data,length(size(data))-1) "Times Vector does not match second last dimension of input data - forgot to epoch, or misspecified 'time' vector?"
        end
    end
    # function content partially taken from MixedModels.jl bootstrap.jl
    df = Array{NamedTuple,1}()
    dataDim = length(size(data)) # surely there is a nicer way to get this but I dont know it

    Xs = modelmatrix(uf)
    # If we have3 dimension, we have a massive univariate linear mixed model for each timepoint
    if dataDim == 3
        firstData = data[1, 1, :]
        ntime = size(data, 2)
    else
        # with only 2 dimension, we run a single time-expanded linear mixed model per channel/voxel
        firstData = data[1, :]
        ntime = 1
    end
    nchan = size(data, 1)

    Xs = (equalizeLengths(Xs[1]),Xs[2:end]...)
    _,data = zeropad(Xs[1],data)
    # get a un-fitted mixed model object
    
    Xs = disallowmissing.(Xs)

    mm = LinearMixedModel_wrapper(Unfold.formula(uf), firstData, Xs)
    # prepare some variables to be used
    βsc, θsc = similar(MixedModels.coef(mm)), similar(mm.θ) # pre allocate
    p, k = length(βsc), length(θsc)
    #β_names = (Symbol.(fixefnames(mm))..., )

    β_names = (Symbol.(vcat(fixefnames(mm)...))...,)
    β_names = (unique(β_names)...,)

    @assert(
        length(β_names) == length(βsc),
        "Beta-Names & coefficient length do not match. Did you provide two identical basis functions?"
    )

    @debug println("beta_names $β_names")
    @debug println("uniquelength: $(length(unique(β_names))) / $(length(β_names))")
    # for each channel
    prog = Progress(nchan * ntime, 0.1)
    #@showprogress .1 
    for ch in range(1, stop = nchan)
        # for each time
        for t in range(1, stop = ntime)

            #@debug "ch:$ch/$nchan, t:$t/$ntime"
            @debug "data-size: $(size(data))"
            #@debug println("mixedModel: $(mm.feterms)")
            if ndims(data) == 3
                MixedModels.refit!(mm, data[ch, t, :])
            else
                MixedModels.refit!(mm, data[ch, :])
            end
            #@debug println(MixedModels.fixef!(βsc,mm))

            β = NamedTuple{β_names}(MixedModels.fixef!(βsc, mm))

            out = (
                objective = mm.objective,
                σ = mm.σ,
                β = NamedTuple{β_names}(MixedModels.fixef!(βsc, mm)),
                se = SVector{p,Float64}(MixedModels.stderror!(βsc, mm)), #SVector not necessary afaik, took over from MixedModels.jl
                θ = SVector{k,Float64}(MixedModels.getθ!(θsc, mm)),
                channel = ch,
                timeIX = ifelse(dataDim == 2, NaN, t),
            )
            push!(df, out)
            ProgressMeter.next!(prog; showvalues = [(:channel, ch), (:time, t)])
        end
    end

    uf.modelfit = UnfoldMixedModelFitCollection(
        df,
        deepcopy(mm.λ),
        getfield.(mm.reterms, :inds),
        copy(mm.optsum.lowerbd),
        NamedTuple{Symbol.(fnames(mm))}(map(t -> (t.cnames...,), mm.reterms)),
    )


    return uf.modelfit
end

function StatsModels.coef(
    uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    beta = [x.β for x in MixedModels.tidyβ(modelfit(uf))]
    return reshape_lmm(uf, beta)
end

function MixedModels.ranef(
    uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    sigma = [x.σ for x in MixedModels.tidyσs(modelfit(uf))]
    return reshape_lmm(uf, sigma)
end

function reshape_lmm(uf::UnfoldLinearMixedModel, est)
    ntime = length(collect(values(design(uf)))[1][2])
    nchan = modelfit(uf).fits[end].channel
    return permutedims(reshape(est, :, ntime, nchan), [3 2 1])
end
function reshape_lmm(uf::UnfoldLinearMixedModelContinuousTime, est)
    nchan = modelfit(uf).fits[end].channel
    return reshape(est, :, nchan)'

end


LinearMixedModel_wrapper(form,data::Array{<:Union{TData},1},Xs;wts = []) where {TData<:Union{Missing,Number}}= @error("currently no support for missing values in MixedModels.jl")

"""
$(SIGNATURES)

Wrapper to generate a LinearMixedModel. Code taken from MixedModels.jl and slightly adapted.

"""
function LinearMixedModel_wrapper(
    form,
    data::Array{<:Union{TData},1},
    Xs;
    wts = [],
) where {TData<:Number}
    #    function LinearMixedModel_wrapper(form,data::Array{<:Union{Missing,TData},1},Xs;wts = []) where {TData<:Number}
    Xs = (equalizeLengths(Xs[1]),Xs[2:end]...)
    # XXX Push this to utilities zeropad
    # Make sure X & y are the same size
    m = size(Xs[1])[1]


    if m != size(data)[1]
        fe,data = zeropad(Xs[1],data)
        
        Xs = changeMatSize!(size(data)[1], fe, Xs[2:end])
    end

    y = (reshape(float(data), (:, 1)))
    
    MixedModels.LinearMixedModel(y, Xs, form, wts)
end

function MixedModels.LinearMixedModel(y, Xs, form::Array, wts)


    form_combined = form[1]
    for f in form[2:end]

        form_combined =
            form_combined.lhs ~
                MatrixTerm(form_combined.rhs[1] + f.rhs[1]) +
                form_combined.rhs[2:end] +
                f.rhs[2:end]
    end
    MixedModels.LinearMixedModel(y, Xs, form_combined, wts)
end


isMixedModelFormula(f::typeof(MixedModels.zerocorr)) = true