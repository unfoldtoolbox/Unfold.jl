"""
$(SIGNATURES)
using DataFrames: AbstractAggregate
Combine two UnfoldDesignmatrices. This allows combination of multiple events.

This also allows to define events with different lengths.

Not supported for models without timebasis, as it is not needed there (one can simply run multiple models)
# Examples
```julia-repl
julia>  basisfunction1 = firbasis(τ=(0,1),sfreq = 10,name="basis1")
julia>  basisfunction2 = firbasis(τ=(0,0.5),sfreq = 10,name="basis2")
julia>  Xdc1          = designmatrix(UnfoldLinearModelContinuousTime(Dict(Any=>(@formula 0~1,basisfunction1)),tbl_1)
julia>  Xdc2          = designmatrix(UnfoldLinearModelContinuousTime(Dict(Any=>(@formula 0~1,basisfunction2)),tbl_2)
julia>  Xdc = Xdc1+Xdc2 
```

"""
function combineDesignmatrices(X1::DesignMatrix, X2::DesignMatrix)

    # the reason for the assertion is simply that I found it too difficult to concatenate the formulas down below ;) should be easy to implement hough
    @assert !(isa(X1.formulas,AbstractArray) && isa(X2.formulas,AbstractArray)) "it is currently not possible to combine desigmatrices from two already concatenated designs - please concatenate one after the other"

    X1 = deepcopy(X1)
    X2 = deepcopy(X2)
    Xs1 = get_Xs(X1)
    Xs2 = get_Xs(X2)
    if typeof(Xs1) <:Vector
        Xcomb = [Xs1...,Xs2]
    else
        Xcomb = [Xs1,Xs2]
    end

    if typeof(X1.Xs) <: Tuple
      Xcomb = lmm_combineMats!(Xcomb,X1,X2)
    end

    if X1.formulas isa FormulaTerm
        # due to the assertion above, we can assume we have only 2 formulas here
        if X1.formulas.rhs isa Unfold.TimeExpandedTerm
            fcomb = Vector{FormulaTerm{<:InterceptTerm, <:TimeExpandedTerm}}(undef,2)
        else
            fcomb = Vector{FormulaTerm}(undef,2) # mass univariate case
        end
        fcomb[1] = X1.formulas
        fcomb[2] = X2.formulas
        return DesignMatrix(fcomb, Xcomb, [X1.events, X2.events])
    else
        if X1.formulas[1].rhs isa Unfold.TimeExpandedTerm
            # we can ignore length of X2, as it has to be a single formula due to the assertion above
            fcomb = Vector{FormulaTerm{<:InterceptTerm, <:TimeExpandedTerm}}(undef,length(X1.formulas)+1)
        else
            fcomb = Vector{FormulaTerm}(undef,length(X1.formulas)+1) # mass univariate case
        end
        fcomb[1:end-1] = X1.formulas
        fcomb[end] = X2.formulas
        return DesignMatrix(fcomb, Xcomb, [X1.events..., X2.events])
    end
end



Base.:+(X1::DesignMatrix, X2::DesignMatrix) = combineDesignmatrices(X1, X2)

Base.:+(X1::Nothing, X2::DesignMatrix) = X2

function get_Xs(Xs::Tuple)
    return Xs[1]
end
function get_Xs(Xs::SparseMatrixCSC)
    return Xs
end
function get_Xs(Xs::AbstractArray)
    # mass univariate case
    # mass univariate case with multiple events
    return Xs
end

"""
Typically returns the field X.Xs of the designmatrix

Compare to `modelmatrix` which further concatenates the designmatrices (in the UnfoldLinearModelContinuousTime) as needed
"""
function get_Xs(X::DesignMatrix)
    return get_Xs(X.Xs)
end

"""
$(SIGNATURES)
designmatrix(type, f, tbl; kwargs...)
Return a *DesignMatrix* used to fit the models.
# Arguments
- type::UnfoldModel
- f::FormulaTerm: Formula to be used in this designmatrix
- tbl: Events (usually a data frame) to be modelled
- basisfunction::BasisFunction: basisfunction to be used in modeling (if specified)
- contrasts::Dict: (optional) contrast to be applied to formula
- eventfields::Array: (optional) Array of symbols which are passed to basisfunction event-wise. 
First field of array always defines eventonset in samples. Default is [:latency]

# Examples
```julia-repl
julia>  designmatrix(UnfoldLinearModelContinuousTime,Dict(Any=>(f,basisfunction1),tbl)
```

"""
function designmatrix(
    type,
    f::Union{Tuple,FormulaTerm},
    tbl,
    basisfunction;
    contrasts = Dict{Symbol,Any}(),
    kwargs...,
)
    @debug("generating DesignMatrix")

    # check for missings-columns - currently missings not really supported in StatsModels
    s_tmp = schema(f,tbl,contrasts) # temporary scheme to get necessary terms
    neededCols = [v.sym for v in values(s_tmp.schema)]

    tbl_nomissing = DataFrame(tbl) # tbl might be a SubDataFrame due to a view - but we can't modify a subdataframe, can we?
    try 
        disallowmissing!(tbl_nomissing,neededCols) # if this fails, it means we have a missing value in a column we need. We do not support this
    catch e
        if e isa ArgumentError
            error(e.msg*"\n we tried to get rid of a event-column declared as type Union{Missing,T}. But there seems to be some actual missing values in there. 
            You have to replace them yourself (e.g. replace(tbl.colWithMissingValue,missing=>0)) or impute them otherwise.")
        else
           rethrow()
        end
    end

        


    form = unfold_apply_schema(type,f, schema(f, tbl_nomissing, contrasts))
    
    @debug "type: $type"
    if (type== UnfoldLinearMixedModel) || (type == UnfoldLinearMixedModelContinuousTime)
        # get all random effects

        check_groupsorting(form.rhs)
    end
    form =
        apply_basisfunction(form, basisfunction, get(Dict(kwargs), :eventfields, nothing))

    # Evaluate the designmatrix
    
    #note that we use tbl again, not tbl_nomissing.
    @debug typeof(form)
    if (!isnothing(basisfunction)) & (type <: UnfoldLinearMixedModel)
        X = modelcols.(form.rhs, Ref(tbl))
    else
        X = modelcols(form.rhs, tbl)
    end
    @debug typeof(X)
    return DesignMatrix(form, X, tbl)
end
function designmatrix(type, f, tbl; kwargs...)
    return designmatrix(type, f, tbl, nothing; kwargs...)
end

"""
wrapper to make apply_schema mixed models as extension possible
"""
unfold_apply_schema(type::Any,f,schema) = apply_schema(f,schema,UnfoldModel)


# specify for abstract interface
designmatrix(uf::UnfoldModel) = uf.designmatrix

function designmatrix(
    uf::UnfoldModel,
    tbl;
    eventcolumn = :event,
    contrasts = Dict{Symbol,Any}(),
    kwargs...,
)
    X = nothing
    fDict = design(uf)
    for (eventname, f) in pairs(fDict)
        
        @debug "Eventname, X:",eventname,X
        if eventname == Any
            eventTbl = tbl
        else
            if !((eventcolumn ∈ names(tbl)) | (eventcolumn ∈ propertynames(tbl)))
                error(
                    "Couldnt find columnName: " *
                    string(eventcolumn) *
                    " in event-table.  Maybe need to specify eventcolumn=:correctColumnName (default is ':event') \n names(tbl) = " *
                    join(names(tbl), ","),
                )
            end
            eventTbl = @view tbl[tbl[:, eventcolumn].==eventname, :] # we need a view so we can remap later if needed
        end
        if isempty(eventTbl)
            error(
                "eventTable empty after subsetting. Couldnt find event '" *
                string(eventname) *
                "'{" *
                string(typeof(eventname)) *
                "}, in field tbl[:,:" *
                string(eventcolumn) *
                "].? - maybe you need to specify it as a string instead of a symbol?",
            )
        end

        fIx = collect(typeof.(f) .<: FormulaTerm)
        bIx = collect(typeof.(f) .<: BasisFunction)



        if any(bIx)
            # timeContinuos way
            # TODO there should be a julian way to do this distinction

            X =
                X + designmatrix(
                    typeof(uf),
                    f[fIx],
                    eventTbl,
                    collect(f[bIx])[1];
                    contrasts = contrasts,
                    kwargs...,
                )
        else
            # normal way
            X =
                X +
                designmatrix(typeof(uf), f[fIx], eventTbl; contrasts = contrasts, kwargs...)
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
function apply_basisfunction(form, basisfunction::BasisFunction, eventfields)
    @debug("apply_basisfunction")
    return FormulaTerm(form.lhs, TimeExpandedTerm(form.rhs, basisfunction, eventfields))
end

function apply_basisfunction(form, basisfunction::Nothing, eventfields)
    # in case of no basisfunctin, do nothing
    return form
end


function designmatrix!(uf::UnfoldModel, evts; kwargs...)
    X = designmatrix(uf, evts; kwargs...)
    uf.designmatrix = X
end

function StatsModels.modelmatrix(uf::UnfoldLinearModel,basisfunction)
    if basisfunction
        @warn("basisfunction not defined for this kind of model")
    else
        return modelmatrix(uf)
    end
end

# catch all case
equalizeLengths(Xs::AbstractMatrix) = Xs

# UnfoldLinearMixedModelContinuousTime case
equalizeLengths(Xs::Tuple) = (equalizeLengths(Xs[1]),Xs[2:end]...)

# UnfoldLinearModel - they have to be equal already
equalizeLengths(Xs::Vector{<:AbstractMatrix}) = Xs 

#UnfoldLinearModelContinuousTime
equalizeLengths(Xs::Vector{<:SparseMatrixCSC}) = equalizeLengths(Xs...) 
equalizeLengths(Xs1::SparseMatrixCSC,Xs2::SparseMatrixCSC,args...) = equalizeLengths(equalizeLengths(Xs1,Xs2),args...)
function equalizeLengths(Xs1::SparseMatrixCSC,Xs2::SparseMatrixCSC)
    sX1 = size(Xs1, 1)
    sX2 = size(Xs2, 1)

    # append 0 to the shorter designmat
    if sX1 < sX2
        Xs1 = SparseMatrixCSC(sX2, Xs1.n, Xs1.colptr, Xs1.rowval, Xs1.nzval)
    elseif sX2 < sX1
        Xs2 = SparseMatrixCSC(sX1, Xs2.n, Xs2.colptr, Xs2.rowval, Xs2.nzval)
    end
    return hcat(Xs1,Xs2)
end
function StatsModels.modelmatrix(uf::UnfoldLinearModelContinuousTime,basisfunction = true)
    if basisfunction
        return modelmatrix(designmatrix(uf))
        #return hcat(modelmatrix(designmatrix(uf))...)
    else
        # replace basisfunction with non-timeexpanded one
        f = formula(uf)
        
        # probably a more julian way to do this...
        if isa(f,AbstractArray)
            return modelcols_nobasis.(f,events(uf))
        else
            return modelcols_nobasis(f,events(uf))
        end
        
    end
end

modelcols_nobasis(f::FormulaTerm,tbl::AbstractDataFrame) = modelcols(f.rhs.term,tbl)
StatsModels.modelmatrix(uf::UnfoldModel) = modelmatrix(designmatrix(uf))#modelmatrix(uf.design,uf.designmatrix.events)
StatsModels.modelmatrix(d::DesignMatrix) = equalizeLengths(d.Xs)


StatsModels.modelmatrix(d::Dict, events) = modelcols(formula(d).rhs, events)

formula(uf::UnfoldModel) = formula(designmatrix(uf))
formula(d::DesignMatrix) = d.formulas


events(uf::UnfoldModel) = events(designmatrix(uf))
events(d::DesignMatrix) = d.events

design(uf::UnfoldModel) = uf.design

function formula(d::Dict) #TODO Specify Dict better
    if length(values(d)) == 1
        return [c[1] for c in collect(values(d))][1]
    else
        return [c[1] for c in collect(values(d))]
    end
    
end

"""
$(SIGNATURES)
calculates the actual designmatrix for a timeexpandedterm. Multiple dispatch on StatsModels.modelcols
"""
function StatsModels.modelcols(term::TimeExpandedTerm, tbl)
    @debug term.term , first(tbl)
    X = modelcols(term.term, tbl)

    time_expand(X, term, tbl)
end




# helper function to get the ranges from where to where the basisfunction is added
function time_expand_getTimeRange(onset, basisfunction)
    npos = sum(times(basisfunction) .>= 0)
    nneg = sum(times(basisfunction) .< 0)

    #basis = kernel(basisfunction)(onset)

    fromRowIx = floor(onset) - nneg
    toRowIx = floor(onset) + npos

    range(fromRowIx, stop = toRowIx)
end


"""
$(SIGNATURES)
performs the actual time-expansion in a sparse way.

 - Get the non-timeexpanded designmatrix X from StatsModels.
 - evaluate the basisfunction kernel at each event
 - calculate the necessary rows, cols and values for the sparse matrix
 Returns SparseMatrixCSC 
"""
function time_expand(X, term, tbl)
    ncolsBasis = size(kernel(term.basisfunction)(0), 2)
    X = reshape(X, size(X, 1), :)

    ncolsX = size(X)[2]
    nrowsX = size(X)[1]
    ncolsXdc = ncolsBasis * ncolsX

    # this is the predefined eventfield, usually "latency"
    #println(term.eventfields)
    #println(tbl[term.eventfields[1]][1:10])
    tbl = DataFrame(tbl)
    onsets = tbl[!, term.eventfields[1]]

    if typeof(term.eventfields) <: Array && length(term.eventfields) == 1
        bases = kernel(term.basisfunction).(tbl[!, term.eventfields[1]])
    else
        bases = kernel(term.basisfunction).(eachrow(tbl[!, term.eventfields]))
    end

    # generate rowindices
    rows = copy(rowvals.(bases))
    # this shift is necessary as some basisfunction time-points can be negative. But a matrix is always from 1:τ. Thus we have to shift it backwards in time.
    # The onsets are onsets-1 XXX not sure why.
    for r = 1:length(rows)
        rows[r] .+= floor(onsets[r] - 1) + shiftOnset(term.basisfunction)
    end

    rows = vcat(rows...)
    rows = repeat(rows, ncolsX)

    # generate column indices
    cols = []
    for Xcol = 1:ncolsX
        for b = 1:length(bases)
            for c = 1:ncolsBasis
                push!(
                    cols,
                    repeat([c + (Xcol - 1) * ncolsBasis], length(nzrange(bases[b], c))),
                )
            end
        end
    end
    cols = vcat(cols...)

    # generate values
    #vals = []
    vals = Array{Union{Missing,Float64}}(undef,size(cols))
    ix = 1
    
    for Xcol = 1:ncolsX
        for (i,b) = enumerate(bases)
            b_nz = nonzeros(b)
            l = length(b_nz)
            
            vals[ix:ix+l-1] .= b_nz .* @view X[i, Xcol]
            ix = ix+l
        #push!(vals, )
        end
    end

    #vals = vcat(vals...)
    ix = rows .> 0 .&& vals .!= 0.
    A = sparse(rows[ix], cols[ix], vals[ix])
    
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
    colnames = Unfold.colnames(term.basisfunction)
    name = Unfold.name(term.basisfunction)
    if typeof(terms) == String
        terms = [terms]
    end
    return name .* " : " .* kron(terms .* " : ", string.(colnames))
end

function termnames(term::TimeExpandedTerm)
    terms = coefnames(term.term)
    colnames = colnames(term.basisfunction)
    if typeof(terms) == String
        terms = [terms]
    end
    return vcat(repeat.([[t] for t in terms], length(colnames))...)
end


function colname_basis(term::TimeExpandedTerm)
    terms = coefnames(term.term)
    colnames = colnames(term.basisfunction)
    if typeof(terms) == String
        terms = [terms]
    end
    return repeat(colnames, length(terms))
end


function StatsModels.coefnames(terms::AbstractArray{<:FormulaTerm})
    return coefnames.(Base.getproperty.(terms, :rhs))
end




function Base.show(io::IO,d::DesignMatrix)
    println(io,"Unfold.DesignMatrix")
    println(io,"Formulas: $(d.formulas)")
    #display(io,d.formulas)
    sz_evts =  isa(d.events,Vector) ? size.(d.events) : size(d.events)
    sz_Xs =  (isa(d.Xs,Vector) | isa(d.Xs,Tuple)) ? size.(d.Xs) : size(d.Xs)
    
    println(io,"\nSizes: Xs: $sz_Xs, events: $sz_evts")
    println(io,"\nuseful functions: formula(d),modelmatrix(d)")
    println(io,"Fields: .formulas, .Xs, .events")
end