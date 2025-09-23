"""
$(SIGNATURES)

Combine two UnfoldDesignmatrices. This allows combination of multiple events.

This also allows to define events with different lengths.

Not supported for models without timebasis, as it is not needed there (one can simply run multiple models)
# Examples
```julia-repl
julia>  basisfunction1 = firbasis(τ=(0,1),sfreq = 10,name="basis1")
julia>  basisfunction2 = firbasis(τ=(0,0.5),sfreq = 10,name="basis2")
julia>  Xdc1          = designmatrix(UnfoldLinearModelContinuousTime([Any=>(@formula 0~1,basisfunction1)],tbl_1)
julia>  Xdc2          = designmatrix(UnfoldLinearModelContinuousTime([Any=>(@formula 0~1,basisfunction2)],tbl_2)
julia>  Xdc = Xdc1+Xdc2
```

"""


Base.:+(X1::Vector{T}, X2::T) where {T<:AbstractDesignMatrix} = [X1..., X2]
Base.:+(X1::T, X2::T) where {T<:AbstractDesignMatrix} = [X1, X2]
Base.:+(X1::Nothing, X2::AbstractDesignMatrix) = [X2]


# helper to get the fixef of lmm but the normal matrix elsewhere

modelmatrices(modelmatrix::AbstractMatrix) = modelmatrix



"""
    modelmatrices(X::AbstractDesignMatrix)
    modelmatrices(X::Vector{<:AbstractDesignMatrix})
    modelmatrices(modelmatrix::AbstractMatrix)

Returns the modelmatrices (also called designmatrices) separately for the events. This is similar to `StatsModels.modelcols`, but merely access the precomputed designmatrix. If the designmatrix needs to be computed, please use `modelcols`

Compare to `modelmatrix` which further concatenates the designmatrices (in the ContinuousTime case).
"""
modelmatrices(X::AbstractDesignMatrix) = modelmatrices(X.modelmatrix)
modelmatrices(X::Vector{<:AbstractDesignMatrix}) = modelmatrices.(X)


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
    unfoldmodeltype::Type{<:UnfoldModel},
    f::Union{Tuple,FormulaTerm},
    tbl,
    basisfunction;
    contrasts = Dict{Symbol,Any}(),
    eventname = Any,
    show_warnings = true,
    kwargs...,
)
    @debug("generating DesignMatrix")

    # check for missings-columns - currently missings not really supported in StatsModels
    s_tmp = schema(f, tbl, contrasts) # temporary scheme to get necessary terms
    neededCols = [v.sym for v in values(s_tmp.schema)]

    tbl_nomissing = DataFrame(tbl) # tbl might be a SubDataFrame due to a view - but we can't modify a subdataframe, can we?
    try
        disallowmissing!(tbl_nomissing, neededCols) # if this fails, it means we have a missing value in a column we need. We do not support this
    catch e
        if e isa ArgumentError
            error(
                e.msg *
                "\n we tried to get rid of a event-column declared as type Union{Missing,T}. But there seems to be some actual missing values in there.
    You have to replace them yourself (e.g. replace(tbl.colWithMissingValue,missing=>0)) or impute them otherwise.",
            )
        else
            rethrow()
        end
    end

    # check datatypes for an Any
    any_ix = eltype.(eachcol(tbl[:, neededCols])) .== Any
    if show_warnings && any(any_ix)
        @warn "Careful, in your `evts` DataFrame you have a column with eltype `Any`: $(neededCols[any_ix]). This will automatically be cast to a Categorical Term -
        if it is supposed to be a continuous term/column, be sure to explicitly cast it to e.g. `Float64` via `data.$(neededCols[any_ix][1]) = Float64.(data.$(neededCols[any_ix][1]))), or to a string if categorical."
    end

    @debug "applying schema $unfoldmodeltype"
    form = unfold_apply_schema(unfoldmodeltype, f, schema(f, tbl_nomissing, contrasts))

    form = apply_basisfunction(
        form,
        basisfunction,
        get(Dict(kwargs), :eventfields, nothing),
        eventname,
    )

    # Evaluate the designmatrix
    X = modelcols(form.rhs, tbl_nomissing)


    designmatrixtype = typeof(designmatrix(unfoldmodeltype())[1])

    #return X
    return designmatrixtype(form, X, tbl)
end

"""
    designmatrix(type, f, tbl; kwargs...)
call without basis function, continue with basisfunction = `nothing`
"""
function designmatrix(type, f, tbl; kwargs...)
    return designmatrix(type, f, tbl, nothing; kwargs...)
end




"""
wrapper to make apply_schema mixed models as extension possible

Note: type is not necessary here, but for LMM it is for multiple dispatch reasons!
"""
unfold_apply_schema(type, f, schema) = apply_schema(f, schema, UnfoldModel)


# specify for abstract interface
designmatrix(uf::UnfoldModel) = uf.designmatrix



"""
    designmatrix(
        uf::UnfoldModel,
        tbl;
        eventcolumn = :event,
        contrasts = Dict{Symbol,Any}(),
        kwargs...,

Main function called from `fit(UnfoldModel...)`, generates the designmatrix, returns a list of `<:AbstractDesignMatrix`

"""
designmatrix(uf::UnfoldModel, tbl; kwargs...) =
    designmatrix(typeof(uf), design(uf), tbl; kwargs...)


"""
    designmatrix(
        T::Type{<:UnfoldModel},
        design_array::Vector{<:Pair},
        tbl;
        eventcolumn = :event,
        contrasts = Dict{Symbol,Any}(),
        kwargs...,

iteratively calls `designmatrix` for each event in the design_array, and returns a list of `<:AbstractDesignMatrix`

"""
function designmatrix(
    T::Type{<:UnfoldModel},
    design_array::Vector{<:Pair},
    tbl;
    eventcolumn = :event,
    contrasts = Dict{Symbol,Any}(),
    kwargs...,
)

    X = nothing
    for (eventname, f) in design_array

        @debug "Eventname, X:", eventname, X
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
            eventTbl = @view tbl[tbl[:, eventcolumn] .== eventname, :] # we need a view so we can remap later if needed
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
                    T,
                    f[fIx],
                    eventTbl,
                    collect(f[bIx])[1];
                    contrasts = contrasts,
                    eventname = eventname,
                    kwargs...,
                )
        else
            # normal way
            @debug f
            X =
                X + designmatrix(
                    T,
                    f[fIx],
                    eventTbl;
                    contrasts = contrasts,
                    eventname = eventname,
                    kwargs...,
                )
        end
    end
    return X
end

import Base.isempty
Base.isempty(d::AbstractDesignMatrix) = isempty(modelmatrices(d))

"""
$(SIGNATURES)
timeexpand the rhs-term of the formula with the basisfunction

"""
function apply_basisfunction(form, basisfunction::BasisFunction, eventfields, eventname)
    @debug "apply_basisfunction" basisfunction.name eventname
    if basisfunction.name == ""
        basisfunction.name = eventname
    elseif basisfunction.name != eventname && eventname != Any

        error(
            "since unfold 0.7 basisfunction names need to be equivalent to the event.name (or = \"\" for autofilling).",
        )
    end
    return FormulaTerm(form.lhs, TimeExpandedTerm(form.rhs, basisfunction, eventfields))
end

function apply_basisfunction(form, basisfunction::Nothing, eventfields, eventname)
    # in case of no basisfunctin, do nothing
    return form
end


function designmatrix!(uf::UnfoldModel{T}, evts; kwargs...) where {T}
    X = designmatrix(uf, evts; kwargs...)
    uf.designmatrix = X
    return uf
end


"""
    modelmatrix(uf::UnfoldLinearModel)
returns the modelmatrix of the model. Concatenates them, except in the MassUnivariate cases, where a vector of modelmatrices is return

Compare with `modelmatrices` which returns a vector of modelmatrices, one per event

"""
function StatsModels.modelmatrix(uf::UnfoldLinearModel, basisfunction)
    if basisfunction
        @warn("basisfunction not defined for this kind of model")
    else
        return modelmatrix(uf)
    end
end


"""
    $(SIGNATURES)

typically takes two modelmatrices, a vector of modelmatrices or a tuple of modelmatrices and horizontally concatenates them. Because for ContinuousTime models the matrices are typically sparse, we try to be more efficient in concatenating them

- For `UnfoldLinearModel` the matrices have to be equal in size already, and we just return the vector
- In case of the `UnfoldLinearMixedModelContinuousTime` case, we recursively call `extend_to_larger` on the FeMat and the ReMats and iteratively adjust the sizes one after another until no `ReMats` are left

"""

# catch all case
extend_to_larger(modelmatrix::AbstractMatrix) = modelmatrix

# UnfoldLinearMixedModelContinuousTime case
extend_to_larger(modelmatrix::Tuple) =
    (extend_to_larger(modelmatrix[1]), modelmatrix[2:end]...)

# UnfoldLinearModel - they have to be equal already
extend_to_larger(modelmatrix::Vector{<:AbstractMatrix}) = modelmatrix

#UnfoldLinearModelContinuousTime
extend_to_larger(modelmatrix::Vector{<:SparseMatrixCSC}) = extend_to_larger(modelmatrix...)
extend_to_larger(modelmatrix1::SparseMatrixCSC, modelmatrix2::SparseMatrixCSC, args...) =
    extend_to_larger(extend_to_larger(modelmatrix1, modelmatrix2), args...)

function extend_to_larger(modelmatrix1::SparseMatrixCSC, modelmatrix2::SparseMatrixCSC)
    sX1 = size(modelmatrix1, 1)
    sX2 = size(modelmatrix2, 1)

    # append 0 to the shorter designmat
    if sX1 < sX2
        modelmatrix1 = SparseMatrixCSC(
            sX2,
            modelmatrix1.n,
            modelmatrix1.colptr,
            modelmatrix1.rowval,
            modelmatrix1.nzval,
        )
    elseif sX2 < sX1
        modelmatrix2 = SparseMatrixCSC(
            sX1,
            modelmatrix2.n,
            modelmatrix2.colptr,
            modelmatrix2.rowval,
            modelmatrix2.nzval,
        )
    end
    return hcat(modelmatrix1, modelmatrix2)
end


"""
    StatsModels.modelmatrix(uf::UnfoldLinearModelContinuousTime, basisfunction = true)
Setting the optional second args to false, will return the modelmatrix without the timeexpansion / basisfunction applied.
"""
function StatsModels.modelmatrix(uf::UnfoldLinearModelContinuousTime, basisfunction = true)
    if basisfunction
        return modelmatrix(designmatrix(uf))
        #return hcat(modelmatrix(designmatrix(uf))...)
    else
        # replace basisfunction with non-timeexpanded one
        f = formulas(uf)

        # probably a more julian way to do this...
        #if isa(f, AbstractArray)
        return modelcols_nobasis.(f, events(uf))
        #else
        #    return modelcols_nobasis(f, events(uf))
        #end

    end
end

modelcols_nobasis(f::FormulaTerm, tbl::AbstractDataFrame) = modelcols(f.rhs.term, tbl)
StatsModels.modelmatrix(uf::UnfoldModel) = modelmatrix(designmatrix(uf))#modelmatrix(uf.design,uf.designmatrix.events)
StatsModels.modelmatrix(d::AbstractDesignMatrix) = modelmatrices(d)
StatsModels.modelmatrix(d::Vector{<:AbstractDesignMatrix}) =
    extend_to_larger(modelmatrices.(d))


"""
    formulas(design::Vector{<:Pair})
returns vector of formulas, no schema has been applied (those formulas never saw the data). Also no timeexpansion has been applied (in the case of timecontinuous models)
"""
formulas(design::Vector{<:Pair}) = formulas.(design)
formulas(design::Pair{<:Any,<:Tuple}) = last(design)[1]



"""
    formulas(uf::UnfoldModel)
    formulas(d::Vector{<:AbstractDesignMatrix})
returns vector of formulas, **after** timeexpansion / apply_schema has been used.
"""

formulas(uf::UnfoldModel) = formulas(designmatrix(uf))
formulas(d::AbstractDesignMatrix) = d.formula
formulas(d::Vector{<:AbstractDesignMatrix}) = formulas.(d)

events(uf::UnfoldModel) = events(designmatrix(uf))
events(d::AbstractDesignMatrix) = d.events
events(d::Vector{<:AbstractDesignMatrix}) = events.(d)

design(uf::UnfoldModel) = uf.design

function formulas(d::Dict) #TODO Specify Dict better
    if length(values(d)) == 1
        return [c[1] for c in collect(values(d))][1]
    else
        return [c[1] for c in collect(values(d))]
    end

end

"""
$(SIGNATURES)
calculates in which rows the individual event-basisfunctions should go in Xdc

timeexpand_rows timeexpand_vals
"""
function timeexpand_rows(onsets, bases, shift, ncolsX)
    # generate rowindices
    rows = copy(rowvals.(bases))

    # this shift is necessary as some basisfunction time-points can be negative. But a matrix is always from 1:τ. Thus we have to shift it backwards in time.
    # The onsets are onsets-1 XXX not sure why.
    for r in eachindex(rows)
        rows[r] .+= floor(onsets[r] - 1) .+ shift
    end


    rows_red = reduce(vcat, rows)
    rows_red = repeat(rows_red, ncolsX)
    return rows_red
end

"""
$(SIGNATURES)
calculates the actual designmatrix for a timeexpandedterm. Multiple dispatch on StatsModels.modelcols
"""
function StatsModels.modelcols(term::TimeExpandedTerm, tbl)
    X = modelcols(term.term, tbl)

    time_expand(X, term, tbl)
end


# helper function to get the ranges from where to where the basisfunction is added
function get_timeexpanded_time_range(onset, basisfunction)
    npos = sum(times(basisfunction) .>= 0)
    nneg = sum(times(basisfunction) .< 0)

    #basis = kernel(basisfunction)(onset)

    fromRowIx = floor(onset) - nneg
    toRowIx = floor(onset) + npos

    range(fromRowIx, stop = toRowIx)
end


function timeexpand_cols_allsamecols(bases, ncolsBasis::Int, ncolsX)
    repeatEach = length(nzrange(bases[1], 1))
    cols_r = UnitRange{Int64}[
        ((1:ncolsBasis) .+ ix * ncolsBasis) for ix in (0:(ncolsX-1)) for b = 1:length(bases)
    ]

    cols = reduce(vcat, cols_r)

    cols = repeat(cols, inner = repeatEach)
    return cols
end

"""
$(SIGNATURES)


calculates in which rows the individual event-basisfunctions should go in Xdc

see also timeexpand_rows timeexpand_vals
"""
function timeexpand_cols(basisfunction, bases, ncolsBasis, ncolsX)
    # we can generate the columns much faster, if all bases output the same number of columns
    fastpath = time_expand_allBasesSameCols(basisfunction, bases, ncolsBasis)
    #fastpath = false

    @debug "fastpath" fastpath
    fastpath = false
    #    @debug isa(basisfunction.scale_duration, Bool)
    if fastpath
        return timeexpand_cols_allsamecols(bases, ncolsBasis, ncolsX)
    else
        return timeexpand_cols_generic(bases, ncolsBasis, ncolsX)
    end
end

function timeexpand_cols_generic(bases, ncolsBasis, ncolsX)
    # it could happen, e.g. for bases that are duration modulated, that each event has different amount of columns
    # in that case, we have to go the slow route
    cols = Vector{Int64}[]
    for Xcol = 1:ncolsX # TODO: potential speedup is to calculate the loop only once and duplicate, only the coloffset is different per 1:ncolsX
        for b = 1:length(bases)
            coloffset = (Xcol - 1) * ncolsBasis
            for c = 1:ncolsBasis
                push!(cols, repeat([c + coloffset], length(nzrange(bases[b], c))))
            end
        end
    end
    return reduce(vcat, cols)


end

function timeexpand_vals(bases, X, nTotal, ncolsX)
    # generate values
    #vals = []
    vals = Array{Union{Missing,Float64}}(undef, nTotal)
    ix = 1

    for Xcol = 1:ncolsX
        for (i, b) in enumerate(bases)
            b_nz = nonzeros(b)
            l = length(b_nz)

            vals[ix:(ix+l-1)] .= b_nz .* @view X[i, Xcol]
            ix = ix + l
            #push!(vals, )
        end
    end
    return vals

end
"""
$(SIGNATURES)
performs the actual time-expansion in a sparse way.

 - Get the non-timeexpanded designmatrix X from StatsModels.
 - evaluate the basisfunction kernel at each event
 - calculate the necessary rows, cols and values for the sparse matrix
 Returns SparseMatrixCSC
"""

function time_expand(Xorg::AbstractArray, term::TimeExpandedTerm, tbl)

    tbl = DataFrame(tbl)

    # this extracts the predefined eventfield, usually "latency"
    onsets = tbl[!, term.eventfields]
    time_expand(Xorg, term.basisfunction, onsets)
end
function time_expand(Xorg::AbstractArray, basisfunction::FIRBasis, onsets)
    #    @debug basisfunction
    if basisfunction.interpolate || size(onsets, 2) > 1
        # code doubling, but I dont know how to do the multiple dispatch otherwise
        # this is the "old" way, pre 0.7.1
        #bases = kernel.(Ref(basisfunction), eachrow(onsets))
        #if size(onsets, 2) != 1
        #    error("not implemented")
        #end
        #bases = kernel.(Ref(basisfunction), onsets[!, 1])
        bases = kernel.(Ref(basisfunction), eachrow(onsets))
        return time_expand(Xorg, basisfunction, onsets[!, 1], bases)

    else
        # fast way
        return time_expand_firdiag(Xorg, basisfunction, onsets)
    end

end

function time_expand_firdiag(Xorg::AbstractMatrix{T}, basisfunction, onsets) where {T}

    #    @debug "Xorg eltype" T
    @assert width(basisfunction) == height(basisfunction)
    w = width(basisfunction)
    adjusted_onsets = Int.(floor.(onsets[!, 1] .+ shift_onset(basisfunction)))

    @assert (minimum(onsets[!, 1]) + shift_onset(basisfunction) .- 1 + w) > 0

    neg_fix = sum(adjusted_onsets[adjusted_onsets .< 1] .- 1)
    ncoeffs = w * size(Xorg, 1) * size(Xorg, 2) .+ neg_fix * size(Xorg, 2)


    I = Vector{Int}(undef, ncoeffs)
    colptr = Vector{Int}(undef, w * size(Xorg, 2) + 1)

    #V = Vector{T}(undef, ncoeffs)
    V = similar(Xorg, ncoeffs)


    for ix_X = 1:size(Xorg, 2)
        #I, J, V, m, n = spdiagm_diag(w, (.-onsets[!, 1] .=> Xorg[:, Xcol])...)

        #ix = (1+size(Xorg, 1)*(ix_X-1)):size(Xorg, 1)*ix_X
        ol = w * size(Xorg, 1) + neg_fix
        #ix = 1:ol
        ix = (1+(ix_X-1)*ol):(ix_X*ol)
        ix_col = (1+(ix_X-1)*w):(ix_X*w+1)

        #@debug ix, ix_col size(I) size(colptr) size(V) w
        @views spdiagm_diag_colwise!(
            I[ix],
            colptr[ix_col],
            V[ix],
            w,
            adjusted_onsets,
            Xorg[:, ix_X];
            colptr_offset = (ix_X - 1) * (w * size(Xorg, 1) + neg_fix),
        )

    end
    colptr[end] = length(V) + 1



    m = adjusted_onsets[end, 1] + w - 1
    n = w * size(Xorg, 2)

    SparseArrays._goodbuffers(Int(m), Int(n), colptr, I, V) || throw(
        ArgumentError(
            "Invalid buffers for SparseMatrixCSC construction n=$n, colptr=$(summary(colptr)), rowval=$(summary(rowval)), nzval=$(summary(nzval))",
        ),
    )

    Ti = promote_type(eltype(colptr), eltype(I))
    SparseArrays.sparse_check_Ti(m, n, Ti)
    SparseArrays.sparse_check(n, colptr, I, V)
    # silently shorten rowval and nzval to usable index positions.
    maxlen = abs(widemul(m, n))
    isbitstype(Ti) && (maxlen = min(maxlen, typemax(Ti) - 1))

    return SparseMatrixCSC(m, n, colptr, I, V)
    #return reduce(hcat, Xdc_list)


end



"""
Even more speed improved version of spdiagm, takes a single float value instead of a vector, like a version of spdiagm that takes in a UniformScaling.
This directly constructs the CSC required fields, thus it is possible to do SparseMatrixCSC(spdiagm_diag_colwise(...)).

Use at your own discretion!!


e.g.
> sz = 5
> ix = [1,3,10]
> vals = [4,1,2]
> spdiagm_diag(sz,ix,vals)
"""

function spdiagm_diag_colwise!(
    I,
    colptr,
    V,
    diagsz,
    onsets,
    values::AbstractVector{T};
    colptr_offset = 0,
) where {T}
    #ncoeffs = diagsz * length(onsets)
    #    @debug size(I) size(colptr) size(V) size(onsets) size(values)
    #colptr .= colptr_offset .+ (1:sum(>(0),onsets):ncoeffs+1)
    r = 1
    _c = 1
    for c = 1:diagsz
        colptr[c] = _c + colptr_offset
        c_add = 0
        for (i, o) in enumerate(onsets)
            #     @debug i, r
            row_val = o + c - 1
            if row_val < 1
                continue
            end
            I[r] = row_val
            V[r] = values[i]
            r = r + 1
            c_add = c_add + 1
        end
        _c = _c + c_add
    end
    return (onsets[end] + diagsz, diagsz, colptr, I, V)
end



function spdiagm_diag_colwise(diagsz, onsets, values::AbstractVector{T}) where {T}

    ncoeffs = diagsz * length(onsets)
    I = Vector{Int}(undef, ncoeffs)
    #J = Vector{Int}(undef, ncoeffs)
    colptr = 1:length(onsets):(ncoeffs+1) |> collect
    V = Vector{T}(undef, ncoeffs)
    r = 1
    for c = 1:diagsz
        for (i, o) in enumerate(onsets)
            v = values[i]


            I[r] = o + c
            #J[r] = c

            V[r] = v
            r = r + 1
        end
    end
    return (onsets[end] + diagsz, diagsz, colptr, I, V)
end

"""
Speed improved version of spdiagm, takes a single float value instead of a vector, like a version of spdiagm that takes in a UniformScaling


e.g.
> sz = 5
> ix = [1,3,10]
> spdiagm_diag(sz,(.-ix.=>1)...)
"""
function spdiagm_diag(diagsz, kv::Pair{<:Integer,T}...) where {T}

    ncoeffs = diagsz * length(kv)
    I = Vector{Int}(undef, ncoeffs)
    J = Vector{Int}(undef, ncoeffs)
    V = Vector{T}(undef, ncoeffs)
    i = 0
    m = 0
    n = 0
    for p in kv
        k = p.first # e.g. -1
        v = p.second # e.g. [1,2,3,4,5]
        if k < 0
            row = -k
            col = 0
        elseif k > 0
            row = 0
            col = k
        else
            row = 0
            col = 0
        end

        r = (1+i):(diagsz+i)
        I[r], J[r] = ((row+1):(row+diagsz), (col+1):(col+diagsz))
        #copyto!(view(V, r), v)
        view(V, r) .= v

        m = max(m, row + diagsz)
        n = max(n, col + diagsz)
        i += diagsz
    end
    return I, J, V, m, n
end


function time_expand(Xorg::AbstractArray, basisfunction::BasisFunction, onsets)
    bases = kernel.(Ref(basisfunction), eachrow(onsets))
    return time_expand(Xorg, basisfunction, onsets, bases)
end

function time_expand(Xorg, basisfunction::BasisFunction, onsets, bases)
    ncolsBasis = width(basisfunction)#size(kernel(basisfunction, 0), 2)::Int64

    X = reshape(Xorg, size(Xorg, 1), :) # why is this necessary?
    ncolsX = size(X)[2]::Int64

    rows = timeexpand_rows(onsets, bases, shift_onset(basisfunction), ncolsX)
    @debug "cols" size(bases) ncolsBasis ncolsX
    cols = timeexpand_cols(basisfunction, bases, ncolsBasis, ncolsX)

    @debug "sizes" size(X) size(rows) size(cols) ncolsX size(bases) ncolsBasis
    vals = timeexpand_vals(bases, X, size(rows), ncolsX)


    #vals = vcat(vals...)
    ix = rows .> 0 #.&& vals .!= 0.
    A = @views sparse(
        rows[ix],
        cols[ix],
        vals[ix],
        maximum(rows),
        ncolsX * width(basisfunction), # we need this, else we loose potential last empty columns
    )
    dropzeros!(A)

    return A
end

"""
Helper function to decide whether all bases have the same number of columns per event
"""
time_expand_allBasesSameCols(b::FIRBasis, bases, ncolBasis) =
    isa(b.scale_duration, Bool) ? true : false  # FIRBasis is always fast if scale_duration = false
function time_expand_allBasesSameCols(basisfunction, bases, ncolsBasis)
    @debug "time_expand_allBasesSameCols generic reached"
    fastpath = true
    for b in eachindex(bases)
        if length(unique(length.(nzrange.(Ref(bases[b]), 1:ncolsBasis)))) != 1
            return false
        end
    end
    return true
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
    return string(name) .* " : " .* kron(terms .* " : ", string.(colnames))
end

function StatsModels.termnames(term::TimeExpandedTerm)
    terms = coefnames(term.term)
    colnames = Unfold.colnames(term.basisfunction)
    if typeof(terms) == String
        terms = [terms]
    end
    return vcat(repeat.([[t] for t in terms], length(colnames))...)
end


function colname_basis(term::TimeExpandedTerm)
    terms = coefnames(term.term)
    colnames = Unfold.colnames(term.basisfunction)
    if typeof(terms) == String
        terms = [terms]
    end
    return repeat(colnames, length(terms))
end


function StatsModels.coefnames(terms::AbstractArray{<:FormulaTerm})
    return coefnames.(Base.getproperty.(terms, :rhs))
end
