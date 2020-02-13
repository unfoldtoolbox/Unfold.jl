struct DesignMatrix{T1<:AbstractMatrix,T2<:AbstractMatrix,T3}#TODO actually the last two are sparse matrices
        formula::FormulaTerm
        X::T1
        Xdc::T2
        Zdc::T3
        #DesignMatrix(formula::FormulaTerm,X::T1,Xdc::T2,Zdc::T3) where {T1,T2,T3} = new{T1,T2,T3}(formula,X,Xdc,Zdc)
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



# struct for behavior
struct TimeExpandedTerm{T} <: AbstractTerm
    term::T
    basisfunction::BasisFunction
    eventtime::Symbol
end

function TimeExpandedTerm(term,basisfunction;eventtime=:latency)
        TimeExpandedTerm(term, basisfunction,eventtime)
end


function Base.show(io::IO, p::TimeExpandedTerm)
        #print(io, "timeexpand($(p.term), $(p.basisfunction.type),$(p.basisfunction.times))")
        println(io,"")
 end


function StatsModels.modelcols(term::TimeExpandedTerm,tbl)

        X = modelcols(term.term,tbl)

        npos = sum(term.basisfunction.times.>=0)
        nneg = sum(term.basisfunction.times.<0)
        X = reshape(X,size(X,1),:)
        ncolsX = size(X)[2]
        nrowsX = size(X)[1]

        ntimes = length(term.basisfunction.times)
        A = spzeros(ceil(maximum(tbl.latency))+npos+1,ntimes*ncolsX)
        for row in 1:nrowsX
                onset = tbl[term.eventtime][row]

                basis = term.basisfunction.kernel(onset)
                for col in 1:ncolsX
                        content = X[row,col]
                        # move it by the event latency

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
        end
        return(A)

end
function DesignMatrix(formula::FormulaTerm,tbl::DataFrame,basisfunction::BasisFunction)

        formula = apply_schema(formula,schema(formula,tbl),LinearMixedModel)


        resp,X = modelcols(formula,tbl)

        A = spzeros(ceil(maximum(tbl.latency))+npos+1,ntimes*ncolsX)


        return(DesignMatrix(formula,tbl,basisfunction,X,Xdc,Zdc))
end
