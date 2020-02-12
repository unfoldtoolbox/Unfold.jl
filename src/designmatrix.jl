# Generating the designmatrix happens here

struct DesignMatrix{T1<:AbstractMatrix,T2<:AbstractMatrix,T3<:AbstractMatrix}#TODO actually the last two are sparse matrices
        formula::FormulaTerm
        tbl::DataFrame
        basisfunction::BasisFunction
        X::T1
        Xdc::T2
        Zdc::T3
        DesignMatrix(formula::FormulaTerm,tbl::DataFrame,basisfunction::BasisFunction,X::T1,Xdc::T2,Zdc::T3) where {T1,T2,T3} = new{T1,T2,T3}(formula,tbl,basisfunction,X,Xdc,Zdc)
end

function Base.show(io::IO, obj::DesignMatrix)
        println(io, "DesignMatrix object")
        println(io, "formula: $(obj.formula)")
        println(io)
        println(io, "tbl: $(size(obj.tbl))")
        println(io)
        println(io, "basisfunction:")
        println(io, "$(obj.basisfunction)")
        println(io)
        println(io, "X: $(size(obj.X))")
        println(io, "Xdc: $(size(obj.Xdc))")
end

function DesignMatrix(formula::FormulaTerm,tbl::DataFrame,basisfunction::BasisFunction)
        tbl[:,:y] .=0
        # TODO Ask if formulas without lhs can exist

        #tbl, _ = StatsModels.missing_omit(tbl, formula)
        #formula = apply_schema(formula,schema(formula,tbl),LinearMixedModel)
        formula = apply_schema(formula,schema(formula,tbl))
        # TODO What does the LinearMixedMOdel do?

        resp,X = modelcols(formula,tbl)
        ncolsX = size(X)[2]
        nrowsX = size(X)[1]

        npos = sum(basisfunction.times.>=0)
        nneg = sum(basisfunction.times.<0)
        ntimes = length(basisfunction.times)
        fixefCols = zeros(ncolsX)
        A = spzeros(ceil(maximum(tbl.latency))+npos+1,ntimes*ncolsX)
        for col in 1:ncolsX
                # TODO:This has big potential for speedup if necesary, as columns share the events
                for row in 1:nrowsX
                        content = X[row,col]
                        # move it by the event latency
                        timing = basisfunction.times
                        onset = tbl[row,:latency]
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

                        A[fromRowIx:toRowIx,fromColIx:toColIx-1] = sparse(Gc)
                end

                fixefCols(col)= !isa(X,ReTerm)

        end
        print(fixefCols)
        Xdc = A
#        Xdc = A[:,fixefCols]
        #Zdc = A[:,!fixefCols]
        Zdc = copy(Xdc)
        return(DesignMatrix(formula,tbl,basisfunction,X,Xdc,Zdc))
end
