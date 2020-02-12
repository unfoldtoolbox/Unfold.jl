# the (sparse) lm fit & mixed model fit functionality is populated here
function fit(ufdesign::unfold.DesignMatrix,y)
if  0#formula has random effects
    out = fit_lmm(ufdesign,y)
else
    out = fit_lm(ufdesign.Xdc,y)
end
return out
end

function fit_lm(X,y)
    # due to our lazy evaluation of X, it could be that it is smaller than y.
    # In that case we fill it with nans

    if size(X,1)<=length(y)
        y = y[1:size(X,1)]
    end

    b,history = lsqr(X,y,log=true)
    return(b,history)
end
