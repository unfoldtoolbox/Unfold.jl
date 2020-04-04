struct TimeExpandedTerm{T<:AbstractTerm} <: AbstractTerm
        term::T
        basisfunction::BasisFunction
        eventfields::Array{Symbol}
end


struct ZeroCorr2{T<:RandomEffectsTerm} <: AbstractTerm
    term::T
end

struct UnfoldDesignmatrix
        formulas
        Xs
end

function combineDesignmatrices(X1::UnfoldDesignmatrix,X2::UnfoldDesignmatrix)
        Xs1 = X1.Xs
        Xs2 = X2.Xs
        if typeof(X1.Xs) <: SparseMatrixCSC
                # easy case

                sX1 = size(Xs1,1)
                sX2 = size(Xs2,1)

                # append 0 to the shorter designmat
                if sX1 < sX2
                        Xs1 = SparseMatrixCSC(sX2, Xs1.n, Xs1.colptr, Xs1.rowval, Xs1.nzval)
                elseif sX2 < sX1
                        Xs2 = SparseMatrixCSC(sX1, Xs2.n, Xs2.colptr, Xs2.rowval, Xs2.nzval)

                end
                Xcomb = hcat(Xs1,Xs2)
        else

        end
        if X1.formulas.rhs.basisfunction.times.step != X2.formulas.rhs.basisfunction.times.step
                @warn("Concatenating formulas with different sampling rates. Be sure that this is what you want.")
        end
        UnfoldDesignmatrix([X1.formulas X2.formulas],Xcomb)
end

Base.:+(X1::UnfoldDesignmatrix, X2::UnfoldDesignmatrix) = combineDesignmatrices(X1,X2)
# with basis expansion
function unfoldDesignmatrix(type,f,tbl,basisfunction;kwargs...)
        Xs,form = generateDesignmatrix(type,f,tbl,basisfunction; kwargs...)
        UnfoldDesignmatrix(form,Xs)
end

#without basis expansion
function unfoldDesignmatrix(type,f,tbl;kwargs...)
        Xs,form = generateDesignmatrix(type,f,tbl,nothing; kwargs...)
        UnfoldDesignmatrix(form,Xs)

end
function generateDesignmatrix(type,f,tbl,basisfunction;contrasts= Dict{Symbol,Any}(), kwargs...)
        form = apply_schema(f, schema(f, tbl, contrasts), LinearMixedModel)
        @debug("generateDesignmatrix")
        if !isnothing(basisfunction)


                if type <: UnfoldLinearMixedModel
                        @debug("Mixed Model")
                        form = FormulaTerm(form.lhs, TimeExpandedTerm.(form.rhs,Ref(basisfunction),get(kwargs,:eventfields,nothing)))
                else
                        @debug("not mixed model, $type")

                        form = FormulaTerm(form.lhs, TimeExpandedTerm(form.rhs,basisfunction,get(kwargs,:eventfields,nothing)))
                end
        end
        X = modelcols(form.rhs, tbl)
        return X,form
end

function TimeExpandedTerm(term,basisfunction,eventfields::Nothing)
        TimeExpandedTerm(term, basisfunction)
end
function TimeExpandedTerm(term,basisfunction,eventfields::Symbol)
        TimeExpandedTerm(term, basisfunction,[eventfields])
end

function TimeExpandedTerm(term,basisfunction;eventfields=:latency)
        TimeExpandedTerm(term, basisfunction,eventfields)
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
#function StatsModels.modelcols(term::TimeExpandedTerm{<:Union{<:RandomEffectsTerm,<:AbstractTerm{<:RandomEffectsTerm}}},tbl)
# exchange this to get ZeroCorr to work



        # get the non-timeexpanded reMat
        reMat = modelcols(term.term,tbl)

        # Timeexpand the designmatrix
        z = transpose(time_expand(transpose(reMat.z),term,tbl))



        # First we check if there is overlap in the timeexpanded term. If so, we cannot continue. Later implementations will remedy that
        #println(dump(term,))
        if hasfield(typeof(term.term),:rhs)
                rhs = term.term.rhs


        elseif hasfield(typeof(term.term.term),:rhs)
                # we probably have something like zerocorr, which does not need to show a .rhs necessarily
                rhs = term.term.term.rhs
        else

                printn("term.term: $(dump(term.term))")
                error("unknown RE structure, has no field .rhs:$(typeof(term.term))")

        end

        group = tbl[rhs.sym]
        time = tbl[term.eventfields[1]]

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
        for j in 1:S
                for i in j:S
                        # We currently restrict to diagonal entries
                        # Once mixedmodels#293 is pushed, we can relax this and use zerocorr()
                        if i == j # for diagonal
                                push!(inds, m[i, j])
                        end
                end
        end

        levels = reMat.levels
        refs =refs


        # reMat.levels doesnt change
        cnames = coefnames(term)
        #print(refs)
        adjA = MixedModels.adjA(refs, z)
        scratch = Matrix{T}(undef, (S, length(uGroup)))

        ReMat{T,S}(rhs,
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

function time_expand(X,term,tbl)

        #to = TimerOutput()
        nColsBasis = size(term.basisfunction.kernel(0),2)



        # We need the


        X = reshape(X,size(X,1),:)

        ncolsX = size(X)[2]
        nrowsX = size(X)[1]
        ncolsXdc = nColsBasis*ncolsX

        onsets = tbl[term.eventfields[1]]


        if typeof(term.eventfields) <:Array && length(term.eventfields) == 1
                bases = term.basisfunction.kernel.(tbl[term.eventfields[1]])
        else

                bases = term.basisfunction.kernel.(eachrow(tbl[term.eventfields]))

        end

        # generate rowindices
        rows =  copy(rowvals.(bases))


        for r in 1:length(rows)
                rows[r] .+= floor(onsets[r]-1)+term.basisfunction.shiftOnset
        end

        rows = vcat(rows...)
        rows = repeat(rows,ncolsX)

        # generate column indices
        cols = []
        #@timeit to "Col"

        for Xcol in 1:ncolsX
                for b in 1:length(bases)
                        for c in 1:nColsBasis
                                push!(cols,repeat([c+(Xcol-1)*nColsBasis],length(nzrange(bases[b],c))))
                        end
                end
        end

        cols = vcat(cols...)

        # generate values
        vals = []
        for Xcol in 1:ncolsX
                push!(vals,vcat(nonzeros.(bases).*X[:,Xcol]...))
        end

        vals = vcat(vals...)
        ix = rows.>0
        #@timeit to "generate"
        A = sparse(rows[ix],cols[ix],vals[ix])
        #println(to)
        return A
end
## Coefnames
function StatsModels.coefnames(term::TimeExpandedTerm)
        names = coefnames(term.term)
        colnames = term.basisfunction.colnames
        if typeof(names) == String
                names = [names]
        end
        return kron(names.*" : ",string.(colnames))
end

function StatsModels.coefnames(term::MixedModels.ZeroCorr)
        coefnames(term.term)
end

function StatsModels.coefnames(term::RandomEffectsTerm)
        coefnames(term.lhs)
end
