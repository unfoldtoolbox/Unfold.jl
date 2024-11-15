calc_Rxy!(R_xy, Xt, data) = mul!(R_xy, Xt, data)
#----- PINV -----
"""
    $(SIGNATURES)
calculates pinv of the designmatrix for later use in the solver-step. This is helpful in case you have many chanels
"""
prepare_pinv(all::Tuple) = prepare_pinv(all...)
prepare_pinv(a, b::AbstractArray, all::Tuple) = prepare_pinv(a, b, all...)
function prepare_pinv(Ĥ, data, Xt, R_xx, R_xy)
    R_xx_pinv = pinv(R_xx)
    return Ĥ, data, (Xt, R_xx_pinv, R_xy)
end

function prepare_pinv(beta, data::AbstractArray{T,3}, X::AbstractArray) where {T}
    return beta, data, (pinv(X),)
end


"""
    $(SIGNATURES)
b .= pinv(X)*data

or

b .= pinv(X'X)*X'y
"""

function solver_pinv!(beta, data, Xt, R_xx_pinv, R_xy)
    calc_Rxy!(R_xy, Xt, data)
    beta .= R_xx_pinv * R_xy
end
function solver_pinv!(beta, data, X_pinv)
    beta .= X_pinv * data
end

function solver_qr!(beta, data, Xt, R_xx_qr, R_xy)
    calc_Rxy!(R_xy, Xt, data)
    beta .= R_xx_qr \ R_xy
end
function solver_qr!(beta, data, X_qr)
    beta .= X_qr \ data
end


#----- Cholesky ---
prepare_cholesky(all::Tuple) = prepare_cholesky(all...)
prepare_cholesky(Ĥ, data, all::Tuple) = prepare_cholesky(Ĥ, data, all...)
prepare_cholesky(Ĥ, data, Xt, R_xx, R_xy) = (Ĥ, data, (Xt, cholesky(R_xx), R_xy))

function solver_cholesky!(beta, data, Xt, XtX_cholesky, R_xy)
    Unfold.calc_Rxy!(R_xy, Xt, data)
    beta .= XtX_cholesky \ R_xy
end

#----- QR -----
prepare_qr(all::Tuple) = prepare_qr(all...)
prepare_qr(a, b::AbstractArray, all::Tuple) = prepare_qr(a, b, all...)
function prepare_qr(Ĥ, data, Xt, R_xx, R_xy)
    #@info typeof(R_xx) typeof(R_xy) typeof(Xt)
    R_xx_qr = qr(R_xx)
    return Ĥ, data, (Xt, R_xx_qr, R_xy)
end

#----- LSMR -----
solver_lsmr!(beta, data, X) = IterativeSolvers.lsmr!(beta, X, data; log = true)[2]


#----- CG -----
function solver_cg!(beta, data, Xt, R_xx, R_xy)
    calc_Rxy!(R_xy, Xt, data)
    _, h = Unfold.IterativeSolvers.cg!(beta, R_xx, R_xy, log = true)
end

#----- Intern -----
function solver_intern!(beta, data, Xt, R_xx, R_xy)
    calc_Rxy!(R_xy, Xt, data)
    beta .= R_xx \ R_xy
end
