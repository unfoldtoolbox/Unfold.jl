# Generating the designmatrix happens here

struct DesignMatrix{T<:AbstractArray}
        formula::FormulaTerm
        events::DataFrame
        basisfunction::BasisFunction
        X::Any
        Xdc::Any
end

function DesignMatrix(formula::FormulaTerm,events::DataFrame,basisfunction::BasisFunction)
        events[:,:y] .=0
        resp,X = modelcols(apply_schema(formula,schema(formula,events)),events)
        ncolsX = size(X)[2]
        nrowsX = size(X)[1]
        # for now we lazyly define the size of Xdc, but I'm sure this will bite me in the ass later
        npos = sum(basisfunction.times.>=0)
        nneg = sum(basisfunction.times.<0)
        ntimes = length(basisfunction.times)
        Xdc = spzeros(ceil(maximum(events.latency))+npos+1,ntimes*ncolsX)
        for col in 1:ncolsX
                for row in 1:nrowsX
                        content = X[row,col]
                        # move it by the event latency
                        timing = basisfunction.times
                        onset = events[row,:latency]
                        basis = basisfunction.kernel(onset)
                        Gc = basis .*content
                        fromRowIx = floor(onset)-nneg
                        toRowIx = floor(onset)+npos
                        # border case of very early event
                        if fromRowIx<1
                                fromRowIx = 1
                                Gc = Gc[sum(fromRowIx.<0):end,:]
                        end
                        fromColIx = 1+(col-1)*ntimes
                        toColIx = fromColIx + ntimes

                        Xdc[fromRowIx:toRowIx,fromColIx:toColIx-1] = sparse(Gc)
                end
        end

        return(formula,events,basisfunction,X,Xdc)

end
