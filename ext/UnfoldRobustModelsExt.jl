module UnfoldRobustModelsExt
using Unfold
using RobustModels
using ProgressMeter
using Missings


function solver_robust(
    X,
    data::AbstractArray{T,3};
    estimator = MEstimator{TukeyLoss}(),
    rlmOptions = (initial_scale = :mad,),
) where {T<:Union{Missing,<:Number}}
    #beta = zeros(Union{Missing,Number},size(data, 1), size(data, 2), size(X, 2))
    beta = zeros(T, size(data, 1), size(data, 2), size(X, 2))
    @showprogress 0.1 for ch = 1:size(data, 1)
        for t = 1:size(data, 2)
            #@debug("$(ndims(data,)),$t,$ch")

            dd = view(data, ch, t, :)
            ix = @. !ismissing(dd)
            # init beta
            #if ch == 1 && t==1
            #elseif ch > 1 && t == 1
            #	copyto!(view(beta, ch, 1,:), view(beta, ch-1, 1,:))
            #else
            # 	copyto!(view(beta, ch,t, :), view(beta, ch,t-1, :))
            #end

            X_local = disallowmissing((X[ix, :])) # view crashes robust model here. XXX follow up once
            # https://github.com/JuliaStats/GLM.jl/issues/470 received a satisfying result
            y_local = disallowmissing(@view(data[ch, t, ix]))

            # if not specified otherwise
            # this results in strange regularizing effects
            #if :initial_coef ∉ keys(rlmOptions)
            #    rlmOptions = merge(rlmOptions,(initial_coef=@view(beta[ch,t,:]),))
            #end
            m = rlm(X_local, y_local, estimator; rlmOptions...)
            beta[ch, t, :] .= coef(m)

        end
    end


    modelfit = Unfold.LinearModelFit{T,3}(beta, ["solver_robust"])

    return modelfit
end

end
