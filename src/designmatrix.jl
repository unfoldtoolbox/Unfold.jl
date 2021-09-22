"""
Object with a *term* and an applicable *BasisFunction* and a eventfield that are later passed to the basisfunction.

$(TYPEDEF)
$(TYPEDSIGNATURES)
$(FIELDS)
# Examples
```julia-repl
julia>  b = TimeExpandedTerm(term,kernel,[:latencyTR,:durationTR])
```
"""
struct TimeExpandedTerm{T<:AbstractTerm} <: AbstractTerm
        "Term that the basis function is applied to. This is regularly called in other functions to get e.g. term-coefnames and timeexpand those"
        term::T
        "Kernel that determines what should happen to the designmatrix of the term"
        basisfunction::BasisFunction
        "Which fields of the event-table should be passed to the basisfunction.Important: The first entry has to be the event-latency in samples!"
        eventfields::Array{Symbol}
end


function TimeExpandedTerm(term::NTuple{N,Union{<:AbstractTerm,<:MixedModels.RandomEffectsTerm}},basisfunction,eventfields::Array{Symbol,1}) where {N}
        # for mixed models, apply it to each Term
        TimeExpandedTerm.(term,Ref(basisfunction),Ref(eventfields))
end
function TimeExpandedTerm(term,basisfunction,eventfields::Nothing)
        # if the eventfield is nothing, call default
        TimeExpandedTerm(term, basisfunction)
end

function TimeExpandedTerm(term,basisfunction,eventfields::Symbol)
        # call with symbol, convert to array of symbol
        TimeExpandedTerm(term, basisfunction,[eventfields])
end

function TimeExpandedTerm(term,basisfunction;eventfields=[:latency])
        # default call, use latency
        TimeExpandedTerm(term, basisfunction,eventfields)
end

function Base.show(io::IO, p::TimeExpandedTerm)
        #print(io, "timeexpand($(p.term), $(p.basisfunction.type),$(p.basisfunction.times))")
        println(io,"$(coefnames(p))")
end



"""
To be explored abstractTerm structure to automatically detect randomeffects. Not in use
"""
struct ZeroCorr2{T<:RandomEffectsTerm} <: AbstractTerm
    term::T
end

"""
Object with a *term* and an applicable *BasisFunction* and a eventfield that are later passed to the basisfunction.

$(TYPEDEF)
$(TYPEDSIGNATURES)
$(FIELDS)
# Examples
```julia-repl
julia>  b = TimeExpandedTerm(term,kernel,[:latencyTR,:durationTR])
```
"""



"""
$(SIGNATURES)
Combine two UnfoldDesignmatrices. This allows combination of multiple events.

This also allows to define events with different lengths.

Not supported for models without timebasis, as it is not needed there (one can simply run multiple models)
# Examples
```julia-repl
julia>  basisfunction1 = firbasis(τ=(0,1),sfreq = 10,name="basis1")
julia>  basisfunction2 = firbasis(τ=(0,0.5),sfreq = 10,name="basis2")
julia>  Xdc1          = designmatrix(UnfoldLinearModel,@formula 0~1,tbl_1,basisfunction1)
julia>  Xdc2          = designmatrix(UnfoldLinearModel,@formula 0~1,tbl_2,basisfunction2)
julia>  combineDesignmatrices(Xdc1,Xdc2)
julia>  Xdc = Xdc1+Xdc2 # equivalently
```

"""
function combineDesignmatrices(X1::DesignMatrix,X2::DesignMatrix)
        
        X1 = deepcopy(X1)
        X2 = deepcopy(X2)
        Xs1 = get_Xs(X1.Xs)
        Xs2 = get_Xs(X2.Xs)

        sX1 = size(Xs1,1)
        sX2 = size(Xs2,1)
        
        # append 0 to the shorter designmat
        if sX1 < sX2
                Xs1 = SparseMatrixCSC(sX2, Xs1.n, Xs1.colptr, Xs1.rowval, Xs1.nzval)
        elseif sX2 < sX1
                Xs2 = SparseMatrixCSC(sX1, Xs2.n, Xs2.colptr, Xs2.rowval, Xs2.nzval)
        end
        
        Xcomb = hcat(Xs1,Xs2)
     
        #if !(Float64(X1.formulas.rhs.basisfunction.times.step) ≈ Float64(X2.formulas.rhs.basisfunction.times.step))
        #        @warn("Concatenating formulas with different sampling rates. Be sure that this is what you want.")
        #end
        if typeof(X1.Xs) <: Tuple
                # we have random effects                
                # combine REMats in single-eventtpe formulas ala y ~ (1|x) + (a|x)
                Xs1 = MixedModels._amalgamate([X1.Xs[2:end]...],Float64)
                Xs2 = MixedModels._amalgamate([X2.Xs[2:end]...],Float64)
                
                Xcomb = (Xcomb,Xs1...,Xs2...)
                
                # Next we make the ranefs all equal size
                equalizeReMatLengths!(Xcomb[2:end])
                
                # check if ranefs can be amalgamated. If this fails, then MixedModels tried to amalgamate over different eventtypes and we should throw the warning
                # if it success, we have to check if the size before and after is identical. If it is not, it tried to amalgamize over different eventtypes which were of the same length
                
                try
                reterms = MixedModels._amalgamate([Xcomb[2:end]...],Float64)
                if length(reterms) != length(Xcomb[2:end])
                        throw("IncompatibleRandomGroupings")
                end
                catch e
                        @error "Error, you seem to have two different eventtypes with the same random-effect grouping variable. \n
                        This is not allowed, you have to rename one. Example:\n
                        eventA: y~1+(1|item) \n
                        eventB: y~1+(1|item)  \n
                        This leads to this error. Rename the later one\n
                        eventB: y~1+(1|itemB) "
                        throw("IncompatibleRandomGroupings")
                end
                


        end
        
        if typeof(X1.formulas) <: FormulaTerm
                DesignMatrix([X1.formulas X2.formulas],Xcomb,[X1.events, X2.events])
        else
                DesignMatrix([X1.formulas... X2.formulas],Xcomb,[X1.events, X2.events])
        end
end

function changeMatSize!(m,fe,remats)
        changeReMatSize!.(remats,Ref(m))
        fe = SparseMatrixCSC(m, fe.n, fe.colptr,fe.rowval, fe.nzval)
        return (fe,remats...,)
end
function changeReMatSize!(remat::MixedModels.AbstractReMat,m ::Integer)
        
        n = m-length(remat.refs) # missing elements
        if n<0
                deleteat!(remat.refs,range(m+1,stop=length(remat.refs)))
                remat.adjA = remat.adjA[:,1:(m+1)]
                remat.wtz = remat.wtz[:,1:(m+1)]
                remat.z = remat.z[:,1:(m+1)]
        elseif n>0
                append!(remat.refs,repeat([remat.refs[end]],n))
                remat.adjA = [remat.adjA sparse(repeat([0],size(remat.adjA)[1],n))]
                remat.wtz = [remat.wtz zeros((size(remat.wtz)[1],n))]
                remat.z = [remat.z zeros((size(remat.z)[1],n))]
        end
        
        #ReMat{T,S}(remat.trm, refs, levels, cnames, z, wtz, λ, inds, adjA, scratch)
end
function equalizeReMatLengths!(remats::NTuple{A,MixedModels.AbstractReMat}) where {A}
        # find max length
        m = maximum([x[1] for x in size.(remats)])
        @debug print("combining lengths: $m")
        # for each reMat
        for k = range(1,length=length(remats))
                remat = remats[k]
                if size(remat)[1] == m
                        continue
                end
                # prolong if necessary
                changeReMatSize!(remat,m)
                
        end
end


Base.:+(X1::DesignMatrix, X2::DesignMatrix) = combineDesignmatrices(X1,X2)

Base.:+(X1::Nothing, X2::DesignMatrix) = X2

function get_Xs(Xs::Tuple)
        return Xs[1]
end
function get_Xs(Xs::SparseMatrixCSC)
        return Xs
end

"""
$(SIGNATURES)
designmatrix(type, f, tbl; kwargs...)
Return a *DesignMatrix* used to fit the models.
# Arguments
- type::Union{UnfoldLinearMixedModel,UnfoldLinearModel}
- f::FormulaTerm: Formula to be used in this designmatrix
- tbl: Events (usually a data frame) to be modelled
- basisfunction::BasisFunction: basisfunction to be used in modeling (if specified)
- contrasts::Dict: (optional) contrast to be applied to formula
- eventfields::Array: (optional) Array of symbols which are passed to basisfunction event-wise. 
First field of array always defines eventonset in samples. Default is [:latency]

# Examples
```julia-repl
julia>  designmatrix(UnfoldLinearModel,f,tbl,basisfunction1)
```

"""

function designmatrix(type,f::Union{Tuple,FormulaTerm},tbl,basisfunction;contrasts= Dict{Symbol,Any}(), kwargs...)
        @debug("generating DesignMatrix")
        form = apply_schema(f, schema(f, tbl, contrasts), MixedModels.LinearMixedModel)
        form = apply_basisfunction(form,basisfunction,get(Dict(kwargs),:eventfields,nothing))

        # Evaluate the designmatrix
        if (!isnothing(basisfunction)) & (type<:UnfoldLinearMixedModel)
                X = modelcols.(form.rhs, Ref(tbl))
        else
                X = modelcols(form.rhs, tbl)
        end
        return DesignMatrix(form,X,tbl)
end
function designmatrix(type,f,tbl;kwargs...)
        return designmatrix(type,f,tbl,nothing;kwargs...)
end

# specify for abstract interface
designmatrix(uf::UnfoldModel) = uf.designmatrix

function designmatrix(uf::UnfoldModel,tbl;eventcolumn = :event,contrasts= Dict{Symbol,Any}(), kwargs...)
        X = nothing
        fDict = design(uf)
        for (eventname,f) in pairs(fDict)        
                if eventname == Any
                        eventTbl = tbl
                else
                        if !(eventcolumn ∈ names(tbl))
                                error("Couldnt find columnName: "*string(eventcolumn)*" in event-table.  Maybe need to specify eventcolumn=:correctColumnName (default is ':event') \n names(tbl) = "*join(names(tbl),","))
                        end
                        eventTbl = tbl[tbl[:,eventcolumn].==eventname,:]
                end
                if isempty(eventTbl)
                        error("eventTable empty after subsetting. Couldnt find event '"*string(eventname)*"'{"*string(typeof(eventname))*"}, in field tbl[:,:"*string(eventcolumn)*"].?")
                end

                fIx = collect(typeof.(f) .<:FormulaTerm)
                bIx = collect(typeof.(f) .<:BasisFunction)
        


                if any(bIx)
                        # timeContinuos way
                        # TODO there should be a julian way to do this distinction

                        X = X+designmatrix(typeof(uf),f[fIx],eventTbl,collect(f[bIx])[1];contrasts= contrasts, kwargs...)
                else
                        # normal way
                        X = X+designmatrix(typeof(uf),f[fIx],eventTbl;contrasts= contrasts, kwargs...)
                end
        end
        return X
end

import Base.isempty
Base.isempty(d::DesignMatrix) = isempty(d.Xs)

"""
$(SIGNATURES)
timeexpand the rhs-term of the formula with the basisfunction

"""

function apply_basisfunction(form,basisfunction::BasisFunction,eventfields)
        @debug("apply_basisfunction")
        return FormulaTerm(form.lhs, TimeExpandedTerm(form.rhs,basisfunction,eventfields))
end

function apply_basisfunction(form,basisfunction::Nothing,eventfields)
        # in case of no basisfunctin, do nothing
        return form
end


function designmatrix!(uf::UnfoldModel,evts;kwargs...)
        X = designmatrix(uf,evts;kwargs...)
        uf.designmatrix = X
    end
    

function StatsModels.modelmatrix(uf::UnfoldLinearModelContinuousTime;basisfunction=true)
        if basisfunction
                return modelmatrix(designmatrix(uf))
        else
                return modelmatrix(design(uf),events(uf))
        end
end


StatsModels.modelmatrix(uf::UnfoldModel) =  modelmatrix(designmatrix(uf))#modelmatrix(uf.design,uf.designmatrix.events)
StatsModels.modelmatrix(d::DesignMatrix) =  d.Xs
StatsModels.modelmatrix(d::Dict,events) = modelcols(formula(d).rhs,events)

formula(uf::UnfoldModel) = formula(designmatrix(uf))
formula(d::DesignMatrix) = d.formulas

events(uf::UnfoldModel) = events(designmatrix(uf))
events(d::DesignMatrix) = d.events

design(uf::UnfoldModel) = uf.design

function formula(d::Dict) #TODO Specify Dict better
        return collect(values(d))[1][1] # give back first formula for now
end

"""
$(SIGNATURES)
calculates the actual designmatrix for a timeexpandedterm. Multiple dispatch on StatsModels.modelcols
"""
# Timeexpand the fixed effect part
function StatsModels.modelcols(term::TimeExpandedTerm,tbl)

        X = modelcols(term.term,tbl)

        time_expand(X,term,tbl)
end


mutable struct SparseReMat{T,S} <: MixedModels.AbstractReMat{T}
        trm
        refs::Vector{Int32}
        levels
        cnames::Vector{String}
        z::SparseMatrixCSC{T}
        wtz::SparseMatrixCSC{T}
        λ::LowerTriangular{T,Matrix{T}}
        inds::Vector{Int}
        adjA::SparseMatrixCSC{T,Int32}
        scratch::Matrix{T}
    end
# This function timeexpands the random effects and generates a ReMat object
#function StatsModels.modelcols(term::TimeExpandedTerm{<:RandomEffectsTerm},tbl)
function StatsModels.modelcols(term::TimeExpandedTerm{<:Union{<:RandomEffectsTerm,<:MixedModels.AbstractReTerm}},tbl)
# exchange this to get ZeroCorr to work


        tbl = DataFrame(tbl)
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

                println("term.term: $(dump(term.term))")
                error("unknown RE structure, has no field .rhs:$(typeof(term.term))")

        end

        group = tbl[!,rhs.sym]
        time = tbl[!,term.eventfields[1]]

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
                refs[Int64(time_start):Int64(time_stop)] .= i
        end

        # Other variables with implementaions taken from the LinerMixedModel function
        wtz = z
        trm = term

        S = size(z, 1)
        T = eltype(z)
        λ  = LowerTriangular(Matrix{T}(I, S, S))

        inds = MixedModels.sizehint!(Int[], (S * (S + 1)) >> 1) # double trouble? can this line go?
        m = reshape(1:abs2(S), (S, S))
        inds = sizehint!(Int[], (S * (S + 1)) >> 1) 
        for j in 1:S
                for i in j:S
                        # We currently restrict to diagonal entries
                        # Once mixedmodels#293 is pushed, we can relax this and use zerocorr()
                        if !(typeof(term.term)<:MixedModels.ZeroCorr) || (i == j) # for diagonal
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

        # Once MixedModels.jl supports it, can be replaced with: SparseReMat
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

"""
$(SIGNATURES)
Get the timeranges where the random grouping variable was applied
"""
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


"""
$(SIGNATURES)
performs the actual time-expansion in a sparse way.

 - Get the non-timeexpanded designmatrix X from StatsModels.
 - evaluate the basisfunction kernel at each event
 - calculate the necessary rows, cols and values for the sparse matrix
 Returns SparseMatrixCSC 
"""
function time_expand(X,term,tbl)
        ncolsBasis = size(term.basisfunction.kernel(0),2)
        X = reshape(X,size(X,1),:)

        ncolsX = size(X)[2]
        nrowsX = size(X)[1]
        ncolsXdc = ncolsBasis*ncolsX

        # this is the predefined eventfield, usually "latency"
        #println(term.eventfields)
        #println(tbl[term.eventfields[1]][1:10])
        tbl = DataFrame(tbl)
        onsets = tbl[!,term.eventfields[1]]

        if typeof(term.eventfields) <:Array && length(term.eventfields) == 1
                bases = term.basisfunction.kernel.(tbl[!,term.eventfields[1]])
        else
                bases = term.basisfunction.kernel.(eachrow(tbl[!,term.eventfields]))
        end

        # generate rowindices
        rows =  copy(rowvals.(bases))
        # this shift is necessary as some basisfunction time-points can be negative. But a matrix is always from 1:τ. Thus we have to shift it backwards in time.
        # The onsets are onsets-1 XXX not sure why.
        for r in 1:length(rows)
                rows[r] .+= floor(onsets[r]-1)+term.basisfunction.shiftOnset
        end

        rows = vcat(rows...)
        rows = repeat(rows,ncolsX)

        # generate column indices
        cols = []
        for Xcol in 1:ncolsX
                for b in 1:length(bases)
                        for c in 1:ncolsBasis
                                push!(cols,repeat([c+(Xcol-1)*ncolsBasis],length(nzrange(bases[b],c))))
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
        A = sparse(rows[ix],cols[ix],vals[ix])

        return A
end

"""
$(SIGNATURES)
coefnames of a TimeExpandedTerm concatenates the basis-function name with the kronecker product of the term name and the basis-function colnames. Separator is ' : '
Some examples for a firbasis:
        basis_313 : (Intercept) : 0.1
        basis_313 : (Intercept) : 0.2
        basis_313 : (Intercept) : 0.3
        ...
"""
function StatsModels.coefnames(term::TimeExpandedTerm)
        terms = coefnames(term.term)
        colnames = term.basisfunction.colnames
        name = term.basisfunction.name
        if typeof(terms) == String
                terms = [terms]
        end
        return name.*" : ".*kron(terms.*" : ",string.(colnames))
end

function termnames(term::TimeExpandedTerm)
        terms = coefnames(term.term)
        colnames = term.basisfunction.colnames
        if typeof(terms) == String
                terms = [terms]
        end
        return vcat(repeat.([[t] for t in terms],length(colnames))...)
end


function colname_basis(term::TimeExpandedTerm)
        terms = coefnames(term.term)
        colnames = term.basisfunction.colnames
        if typeof(terms) == String
                terms = [terms]
        end
        return repeat(colnames,length(terms))
end


function StatsModels.coefnames(terms::AbstractArray{<:FormulaTerm})
        return coefnames.(Base.getproperty.(terms,:rhs))
end
function StatsModels.coefnames(term::MixedModels.ZeroCorr)
        coefnames(term.term)
end

function StatsModels.coefnames(term::RandomEffectsTerm)
        coefnames(term.lhs)
end
