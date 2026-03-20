"""
    solver_cv(; k=5, shuffle=true, inner_solver=(X,y)->Unfold.solver_default(X,y), fit_test=false)

A cross-validation solver for Unfold. This function generates a solver that can be passed to `Unfold.fit` to perform k-fold cross-validation on the data.
# Arguments
- `n_folds::Int`: Number of folds for cross-validation (default: 5).
- `shuffle::Bool`: Whether to shuffle the data before splitting into folds (default: true).
- `inner_solver::Function`: A function that takes a design matrix `X` and data `y` and returns a model fit. This is the solver that will be used for each fold (default: `Unfold.solver_default`).
- `fit_test::Bool`: Whether to fit the model on the test set as well and return those estimates (default: false).
# Returns
A function that can be used as a solver in `Unfold.fit`, which performs k-fold cross-validation and returns a `LinearModelFitCV` object containing the estimates, standard errors, and fold indices.
"""
function solver_cv(;
    rng = MersenneTwister(),
    n_folds = 5,
    shuffle = true,
    inner_solver = (X, y) -> Unfold.solver_default(X, y),
    fit_test = false,
)

    # constructing the solver-function called by Unfold's fit()
    function cv_kernel(X, y)
        #@assert length(_X) == 1 "multi-event CV not supported right now"
        #X = _X[1]
        n_trials = size(y, 3)
        indices = 1:n_trials

        # Generate split indices
        folds_iterator = kfolds(shuffle ? shuffleobs(rng, indices) : indices, k = n_folds)

        train_results = []
        test_results = [] # Only used if fit_test=true
        fold_indices = []

        for (i, (train_idx, test_idx)) in enumerate(folds_iterator)
            # Split the data
            res_train = inner_solver(X[train_idx, :], y[:, :, train_idx])
            push!(train_results, res_train)

            if fit_test
                res_test = inner_solver(X[test_idx, :], y[:, :, test_idx])
                push!(test_results, res_test)
            end
            push!(fold_indices, (train = train_idx, test = test_idx))
        end

        # 3. Concatenate results along a new 4th dimension (the folds)
        # fold_results[i].estimate is (chan, time, coeff)
        all_estimates = cat([r.estimate for r in train_results]..., dims = 4)
        all_stderror = cat([r.standarderror for r in train_results]..., dims = 4)

        # Prepare Info field
        info_dict = (
            inner_infos = [r.info for r in train_results],
            test_estimates = fit_test ?
                             cat([r.estimate for r in test_results]..., dims = 4) :
                             nothing,
        )
        return LinearModelFitCV(
            estimate = all_estimates,
            info = info_dict,
            standarderror = all_stderror,
            folds = fold_indices,
        )
    end

    return cv_kernel
end

coef(m::LinearModelFitCV) = dropdims(mean(m.estimate, dims = 4), dims = 4)
#stderror(m::LinearModelFitCV) = dropdims(mean(m.standarderror, dims = 4), dims = 4)
function stderror(m::LinearModelFitCV{T,N}) where {T,N}
    if isempty(m.standarderror)
        stderror = fill(nothing, size(m.estimate)[1:3])
    else
        stderror = dropdims(mean(T.(m.standarderror), dims = 4), dims = 4)
    end
end
