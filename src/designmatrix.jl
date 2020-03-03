using StatsModels

struct UnfoldDesignmatrix
        formulas
        Xs
end

# with basis expansion
function unfoldDesignmatrix(type,f,tbl,basisfunction;kwargs...)
        println(kwargs)
        Xs,form = generateDesignmatrix(type,f,tbl,basisfunction; kwargs...)
        UnfoldDesignmatrix(form,Xs)
end

#without basis expansion
function unfoldDesignmatrix(type,f,tbl;kwargs...)
        Xs,form = generateDesignmatrix(type,f,tbl,nothing; kwargs...)
        UnfoldDesignmatrix(form,Xs)

end
function generateDesignmatrix(type,f,tbl,basisfunction;contrasts= Dict{Symbol,Any}())
        form = apply_schema(f, schema(f, tbl, contrasts), LinearMixedModel)
        println("generateDesignmatrix")
        if !isnothing(basisfunction)

                println(typeof(form.rhs))
                if type <: UnfoldLinearMixedModel
                        println("Mixed Model")
                        form = FormulaTerm(form.lhs, TimeExpandedTerm.(form.rhs,Ref(basisfunction)))
                else
                        println("not mixed model, $type")
                        form = FormulaTerm(form.lhs, TimeExpandedTerm(form.rhs,basisfunction))
                end
        end
        X = modelcols(form.rhs, tbl)
        return X,form
end


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
        println(io,"$(coefnames(p))")
end




# Timeexpand the fixed effect part
function StatsModels.modelcols(term::TimeExpandedTerm,tbl)
        X = modelcols(term.term,tbl)
        time_expand(X,term,tbl)
end

# This function timeexpands the random effects and generates a ReMat object
function StatsModels.modelcols(term::TimeExpandedTerm{<:RandomEffectsTerm},tbl)
        ntimes = length(term.basisfunction.times)

        # get the non-timeexpanded reMat
        reMat = modelcols(term.term,tbl)

        # Timeexpand the designmatrix
        z = transpose(time_expand(transpose(reMat.z),term,tbl))



        # First we check if there is overlap in the timeexpanded term. If so, we cannot continue. Later implementations will remedy that
        group = tbl[term.term.rhs.sym]
        time = tbl[term.eventtime]

        # get the from-to onsets of the grouping varibales
        onsets = time_expand_getRandomGrouping(group,time,term.basisfunction)
        #print(size(reMat.z))
        refs = zeros(size(z)[2]).+1
        for (i,o) in enumerate(onsets[2:end])
                # check for overlap
                if (minimum(o) <= maximum(onsets[i+1])) & (maximum(o) <= minimum(onsets[i+1]))
                        error("overlap in random effects structure detected, not currently supported")
                end
        end

        # From now we can assume no overlap
        # We want to fnd which subject is active when
        refs = zeros(size(z)[2]).+1
        uGroup = unique(group)

        for (i,g) = enumerate(uGroup[1:end])

                ix_start = findfirst(g.==group)
                ix_end = findlast(g.==group)
                if i == 1
                        time_start = 1
                else
                        time_start = time[ix_start]
                        time_start = time_start - sum(term.basisfunction.times.<=0)
                end
                if i == length(uGroup)
                        time_stop = size(refs,1)
                else
                        time_stop = time[ix_end]
                        time_stop = time_stop + sum(term.basisfunction.times.>0)
                end
                if time_start < 0
                        time_start = 1
                end

                if time_stop > size(refs,1)
                        time_stop = size(refs,1)
                end


                #println("$g,$time_start,$time_stop")
                refs[Int64(time_start):Int64(time_stop)] .= g
        end

        # Other variables with implementaions taken from the LinerMixedModel function
        wtz = z
        trm = term

        S = size(z, 1)
        T = eltype(z)
        λ  = LowerTriangular(Matrix{T}(I, S, S))

        inds = MixedModels.sizehint!(Int[], (S * (S + 1)) >> 1)
        m = reshape(1:abs2(S), (S, S))
        inds = sizehint!(Int[], (S * (S + 1)) >> 1)
        for j = 1:S, i = j:S
                # We currently restrict to diagonal entries
                if i == j # for diagonal
                        push!(inds, m[i, j])
                end
        end

        levels = reMat.levels
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

# Get the timeranges where the random grouping variable was applied
function time_expand_getRandomGrouping(tblGroup,tblLatencies,basisfunction)
        ranges = time_expand_getTimeRange.(tblLatencies,Ref(basisfunction))
end

# helper function to get the ranges from where to where the basisfunction is added
function time_expand_getTimeRange(onset,basisfunction)
        npos = sum(basisfunction.times.>=0)
        nneg = sum(basisfunction.times.<0)

        basis = basisfunction.kernel(onset)

        fromRowIx = floor(onset)-nneg
        toRowIx = floor(onset)+npos

        range(fromRowIx,stop=toRowIx)
end


# Applies the timebasis kernel saved in the "term"
function time_expand(X,term,tbl)
        npos = sum(term.basisfunction.times.>=0)
        nneg = sum(term.basisfunction.times.<0)
        # make sure that this is a 2d matrix
        #println("timeexpand")
        #println(typeof(X))
        #if typeof(X) == NTuple{N,Union{<:Array{Float64,2}, <:ReMat}} where {N}
#                println("selecting first")
#                X = X[1]#choose the MatrixTerm
#        end
#        println(typeof(X))
        X = reshape(X,size(X,1),:)
        ncolsX = size(X)[2]
        nrowsX = size(X)[1]

        ntimes = length(term.basisfunction.times)
        A = spzeros(ceil(maximum(tbl.latency))+npos+1,ntimes*ncolsX)

        onsets = tbl[term.eventtime]
        bases = term.basisfunction.kernel.(onsets)

        for row in 1:nrowsX
                onset = tbl[term.eventtime][row]

                basis = term.basisfunction.kernel(onset)
                for col in 1:ncolsX
                        fromRowIx = floor(onsets[row])-nneg
                        toRowIx   = floor(onsets[row])+npos


                        content = X[row,col]

                        Gc = bases[row] .*content
                        # border case of very early event
                        if fromRowIx<1
                                tmp = (abs(fromRowIx)+2)
                                Gc = Gc[tmp:end,:]
                                fromRowIx = 1
                        end
                        fromColIx = 1+(col-1)*ntimes
                        toColIx = fromColIx + ntimes
                        tmp = A[fromRowIx:toRowIx,fromColIx:toColIx-1]+Gc
                        A[fromRowIx:toRowIx,fromColIx:toColIx-1] = tmp
                end
        end
        return(A)
end
# Applies the timebasis kernel saved in the "term"
function time_expand_genSparse(X,term,tbl)
        # implementation that generates a fresh sparse matrix instead of adding to a large one
        error('this implementation is not working / tested properly, it is also not faster ')

        npos = sum(term.basisfunction.times.>=0)
        nneg = sum(term.basisfunction.times.<0)
        # make sure that this is a 2d matrix

        X = reshape(X,size(X,1),:)
        ncolsX = size(X)[2]
        nrowsX = size(X)[1]
        to = TimerOutput()

        ntimes = length(term.basisfunction.times)
        A = spzeros(ceil(maximum(tbl.latency))+npos+1,ntimes*ncolsX)


        # calculate the basisfunction for each onset
        onsets = tbl[term.eventtime]
        @timeit to "getbasis" Gc = term.basisfunction.kernel.(onsets)

        fromRowIx = floor.(onsets).-nneg
        toRowIx = floor.(onsets).+npos
        println(size(Gc))
        println(size(X))



        @timeit to "colLoop" for col in 1:ncolsX # predictors
                Gc = Gc .* X[:,col]
                println(size(Gc))

                @timeit to "rowLimit" for row in 1:nrowsX # events
                        if  fromRowIx[row] .< 1
                                tmp = (abs(fromRowIx[row])+2)
                                Gc[row] = Gc[row][tmp:end,:]
                                fromRowIx[row] = 1
                        end
                end


                fromColIx = 1+(col-1)*ntimes
                toColIx = fromColIx + ntimes

                Gc_flat = vcat(Gc...)

                rowix = []
                @timeit to "rowappend" for row in 1:nrowsX
                        r = range(fromRowIx[row],stop=toRowIx[row])
                        append!(rowix,r)
                end

                ix_un = unique(rowix)
                rowix_un = []
                colix_un = []
                content_un = fill(0.,(size(ix_un,1),ntimes))

                @timeit to "rowmatch" for r in 1:length(ix_un)
                        ix_tmp = ix_un[r] .== ix_un
                        #println(dump(ix_tmp))
                        tmp = sum(Gc_flat[ix_tmp,:],dims=1)
                        #println(dump(tmp))

                        content_un[r,:] = tmp
                end
                @timeit to "columnExpand" for cTimeexp in 1:ntimes
                        append!(colix_un, repeat([cTimeexp+(ncolsX-1)*ntimes],size(ix_un,1)))
                        append!(rowix_un,ix_un)

                end

                println("row $(size(rowix_un)), col $(size(colix_un)),content $(size(content_un))")
                tmp = content_un[:]
                @timeit to "genSparse" A = sparse(rowix_un,colix_un,tmp)
                #@timeit to "colAssign" A[ix_un,fromColIx:toColIx-1] = content_un


        end
        println(to)
        return(A)
end

## Coefnames
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
