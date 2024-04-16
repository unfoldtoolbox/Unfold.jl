StatsModels.modelmatrix(
    uf::Union{UnfoldLinearMixedModelContinuousTime,<:UnfoldLinearMixedModel},
) = modelmatrix(designmatrix(uf))


function StatsModels.modelmatrix(
    Xs::Vector{
        <:Union{DesignMatrixLinearMixedModel,<:DesignMatrixLinearMixedModelContinuousTime},
    },
)
    @debug "modelmatrix" typeof(Xs)

    #X_vec = getfield.(designmatrix(uf), :modelmatrix)
    Xcomb = Xs[1]
    for k = 2:length(Xs)
        @debug typeof(Xcomb) typeof(Xs[k])
        modelmatrix1 = Unfold.modelmatrices(Xcomb)
        modelmatrix2 = Unfold.modelmatrices(Xs[k])

        @debug typeof(modelmatrix1), typeof(modelmatrix2)
        Xcomb_temp = Unfold.extend_to_larger(modelmatrix1, modelmatrix2)
        @debug "tmp" typeof(Xcomb_temp)
        Xcomb = lmm_combine_modelmatrices!(Xcomb_temp, Xcomb, Xs[k])
        @debug "Xcomb" typeof(Xcomb)
    end
    Xs = length(Xs) > 1 ? Xcomb : [Xs[1].modelmatrix]
    return Xs
end
"""
fit!(uf::UnfoldModel,data::Union{<:AbstractArray{T,2},<:AbstractArray{T,3}}) where {T<:Union{Missing, <:Number}}

Fit a DesignMatrix against a 2D/3D Array data along its last dimension
Data is typically interpreted as channel x time (with basisfunctions) or channel x time x epoch (for mass univariate)

- `show_progress` (default:true), deactivate the progressmeter

Returns an UnfoldModel object

# Examples
```julia-repl
```

"""
function StatsModels.fit!(
    uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
    data::AbstractArray{T};
    show_progress = true,
    kwargs...,
) where {T}

    #@assert length(first(values(design(uf)))[2])
    if uf isa UnfoldLinearMixedModel
        if ~isempty(Unfold.design(uf))
            @assert length(Unfold.times(Unfold.design(uf))[1]) ==
                    size(data, length(size(data)) - 1) "Times Vector does not match second last dimension of input data - forgot to epoch, or misspecified 'time' vector?"
        end
    end
    # function content partially taken from MixedModels.jl bootstrap.jl
    df = Array{NamedTuple,1}()
    dataDim = length(size(data)) # surely there is a nicer way to get this but I dont know it
    @debug typeof(uf)
    #Xs = modelmatrix(uf)
    Xs = modelmatrix(uf)
    if isa(Xs, Vector)
        Xs = Xs[1]
    end
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


    Xs = (Unfold.extend_to_larger(Xs[1]), Xs[2:end]...)#(Unfold.extend_to_larger(Xs[1]), Xs[2:end]...)
    _, data = Unfold.equalize_size(Xs[1], data)
    # get a un-fitted mixed model object

    Xs = (disallowmissing(Xs[1]), Xs[2:end]...)
    #Xs = (Matrix(Xs[1]),Xs[2:end]...)
    @debug "firstdata" size(firstData)
    mm = LinearMixedModel_wrapper(Unfold.formulas(uf), firstData, Xs)
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

    #@debug println("beta_names $β_names")
    #@debug println("uniquelength: $(length(unique(β_names))) / $(length(β_names))")
    # for each channel
    prog = Progress(nchan * ntime, 0.1)
    #@showprogress .1 
    for ch in range(1, stop = nchan)
        # for each time
        for t in range(1, stop = ntime)

            #@debug "ch:$ch/$nchan, t:$t/$ntime"
            #@debug "data-size: $(size(data))"
            #@debug println("mixedModel: $(mm.feterms)")
            if ndims(data) == 3
                MixedModels.refit!(mm, data[ch, t, :]; progress = false)
            else
                #@debug size(mm.y)
                MixedModels.refit!(mm, data[ch, :]; progress = false)
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
            if show_progress
                ProgressMeter.next!(prog; showvalues = [(:channel, ch), (:time, t)])
            end
        end
    end

    uf.modelfit = UnfoldLinearMixedModelFit{T,ndims(data)}(
        LinearMixedModelFitCollection{T}(
            df,
            deepcopy(mm.λ),
            getfield.(mm.reterms, :inds),
            copy(mm.optsum.lowerbd),
            NamedTuple{Symbol.(fnames(mm))}(map(t -> (t.cnames...,), mm.reterms)),
        ),
    )


    return uf
end

function StatsModels.coef(
    uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    beta = [x.β for x in MixedModels.tidyβ(uf)]
    return reshape_lmm(uf, beta)
end

function MixedModels.ranef(
    uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    sigma = [x.σ for x in MixedModels.tidyσs(uf)]
    return reshape_lmm(uf, sigma)
end

function reshape_lmm(uf::UnfoldLinearMixedModel, est)
    ntime = length(Unfold.times(uf)[1])
    @debug ntime
    nchan = modelfit(uf).fits[end].channel
    return permutedims(reshape(est, :, ntime, nchan), [3 2 1])
end
function reshape_lmm(uf::UnfoldLinearMixedModelContinuousTime, est)
    nchan = modelfit(uf).fits[end].channel
    return reshape(est, :, nchan)'

end


LinearMixedModel_wrapper(
    form,
    data::Array{<:Union{TData},1},
    Xs;
    wts = [],
) where {TData<:Union{Missing,Number}} =
    @error("currently no support for missing values in MixedModels.jl")

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
    @debug "LMM wrapper, $(typeof(Xs))"
    Xs = (Unfold.extend_to_larger(Xs[1]), Xs[2:end]...)
    # XXX Push this to utilities equalize_size
    # Make sure X & y are the same size
    @assert isa(Xs[1], AbstractMatrix) & isa(Xs[2], ReMat) "Xs[1] was a $(typeof(Xs[1])), should be a AbstractMatrix, and Xs[2] was a $(typeof(Xs[2])) should be a ReMat"
    m = size(Xs[1])[1]


    if m != size(data)[1]
        fe, data = Unfold.equalize_size(Xs[1], data)

        Xs = change_modelmatrix_size!(size(data)[1], fe, Xs[2:end])
    end

    #y = (reshape(float(data), (:, 1)))
    y = data

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
    #@debug typeof(form_combined)
    @debug typeof(y), typeof(Xs), typeof(wts)
    MixedModels.LinearMixedModel(y, Xs, form_combined, wts)
end


isa_lmm_formula(f::typeof(MixedModels.zerocorr)) = true
