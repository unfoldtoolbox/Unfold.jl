# struct for behavior
struct TimeExpandedTerm{T} <: AbstractTerm
    term::T
    basisfunction::BasisFunction
    eventtime::Symbol
end

function StatsModels.coefnames(term::TimeExpandedTerm)
        names = coefnames(term.term)
        times = term.basisfunction.times
        if typeof(names) == String
                names = [names]
        end
        return kron(names.*" : ",string.(times))
end
function StatsModels.coefnames(term::MixedModels.ZeroCorr)
        coefnames(term.term)
end
function StatsModels.coefnames(term::RandomEffectsTerm)
        coefnames(term.lhs)
end
function TimeExpandedTerm(term,basisfunction;eventtime=:latency)
        TimeExpandedTerm(term, basisfunction,eventtime)
end



function Base.show(io::IO, p::TimeExpandedTerm)
        #print(io, "timeexpand($(p.term), $(p.basisfunction.type),$(p.basisfunction.times))")
        println(io,"$(coefnames(p))")
 end



function StatsModels.modelcols(term::TimeExpandedTerm,tbl)
        X = modelcols(term.term,tbl)
        time_expand(X,term,tbl)


end
function StatsModels.modelcols(term::TimeExpandedTerm{<:RandomEffectsTerm},tbl)
        # to start with, timeexpand the random effect lhs
        reMat = modelcols(term.term,tbl)

        z = transpose(time_expand(transpose(reMat.z),term,tbl))
        ntimes = length(term.basisfunction.times)

        #check for overlap
        group = tbl[term.term.rhs.sym]
        time = tbl[term.eventtime]
        onsets = time_expand_getRandomGrouping(group,time,term.basisfunction)
        #print(size(reMat.z))
        refs = zeros(size(z)[2]).+1
        for (i,o) in enumerate(onsets[2:end])
                if (minimum(o) <= maximum(onsets[i+1])) & (maximum(o) <= minimum(onsets[i+1]))
                        error("overlap in random effects structure detected, not currently supported")
                end
                refs[o] .= tbl[term.term.rhs.sym][i]
        end

        # because the above method with refs did not work, because he complains about not increasing values, I do it hacky:
        # we can assume no overlap by here!
        refs = zeros(size(z)[2]).+1
        uGroup = unique(group)
        for (i,g) = enumerate(uGroup[1:end-1])
                ix_start = findfirst(g.==group)
                if i == 1
                        ix_start = 1
                end
                ix_end = findfirst(uGroup[i+1].==group)

                refs[Int64(time[ix_start]):Int64(time[ix_end])] .= g
        end
        # due to local scope
        ix_end = findlast(uGroup[end-1].==group)
        refs[Int64(time[ix_end]):end] .= uGroup[end]

        wtz = z
        trm = term

        S = size(z, 1)
        T = eltype(z)
        λ  = LowerTriangular(Matrix{T}(I, S, S))

        inds = MixedModels.sizehint!(Int[], (S * (S + 1)) >> 1)
        m = reshape(1:abs2(S), (S, S))
        inds = sizehint!(Int[], (S * (S + 1)) >> 1)
        for j = 1:S, i = j:S
                if i == j # for diagonal
                        push!(inds, m[i, j])
                end
        end



        levels = reMat.levels#vcat(transpose(hcat(repeat([reMat.levels],1,ntimes)...))...)
        refs =refs


        # reMat.levels doesnt change
        cnames = coefnames(term)
        #print(refs)
        adjA = MixedModels.adjA(refs, z)
        scratch = Matrix{T}(undef, (S, length(uGroup)))

        ReMat{T,S}(term.term.rhs,
        refs,
        levels,
        cnames,
        z,
        wtz,
        λ,
        inds,
        adjA,
        scratch)

end


function time_expand_getRandomGrouping(tblGroup,tblLatencies,basisfunction)
        ranges = time_expand_getTimeRange.(tblLatencies,Ref(basisfunction))


end
function time_expand_getTimeRange(onset,basisfunction)
        npos = sum(basisfunction.times.>=0)
        nneg = sum(basisfunction.times.<0)

        basis = basisfunction.kernel(onset)

        fromRowIx = floor(onset)-nneg
        toRowIx = floor(onset)+npos

        range(fromRowIx,stop=toRowIx)

end
function time_expand(X,term,tbl)

        npos = sum(term.basisfunction.times.>=0)
        nneg = sum(term.basisfunction.times.<0)
        # make sure that this is a 2d matrix
        X = reshape(X,size(X,1),:)
        ncolsX = size(X)[2]
        nrowsX = size(X)[1]

        ntimes = length(term.basisfunction.times)
        A = spzeros(ceil(maximum(tbl.latency))+npos+1,ntimes*ncolsX)
        for row in 1:nrowsX
                onset = tbl[term.eventtime][row]

                basis = term.basisfunction.kernel(onset)
                for col in 1:ncolsX
                        fromRowIx = floor(onset)-nneg
                        toRowIx = floor(onset)+npos

                        content = X[row,col]
                        # move it by the event latency

                        Gc = basis .*content

                        # border case of very early event
                        if fromRowIx<1
                                tmp = (abs(fromRowIx)+2)
                                Gc = Gc[tmp:end,:]
                                fromRowIx = 1
                        end
                        fromColIx = 1+(col-1)*ntimes
                        toColIx = fromColIx + ntimes

                        A[fromRowIx:toRowIx,fromColIx:toColIx-1] = A[fromRowIx:toRowIx,fromColIx:toColIx-1]+sparse(Gc)
                end
        end
        return(A)

end
