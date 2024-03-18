
function StatsModels.coefnames(term::RandomEffectsTerm)
    coefnames(term.lhs)
end

check_groupsorting(r::MatrixTerm) = check_groupsorting(r.terms)
function check_groupsorting(r::Tuple)
    @debug "checking group sorting"
    ix = findall(isa.(r, MixedModels.AbstractReTerm))

    rhs(x::RandomEffectsTerm) = x.rhs
    rhs(x::MixedModels.ZeroCorr) = rhs(x.term)
    groupvars = [map(x -> rhs(x).sym, r[ix])...]


    @assert groupvars == sort(groupvars) "random effects have to be alphabetically ordered. e.g. (1+a|X) + (1+a|A) is not allowed. Please reorder"
end
function Unfold.unfold_apply_schema(
    type::Type{<:Union{<:UnfoldLinearMixedModel,<:UnfoldLinearMixedModelContinuousTime}},
    f,
    schema,
)
    @debug "LMM apply schema"
    f_new = apply_schema(f, schema, MixedModels.LinearMixedModel)
    check_groupsorting(f_new.rhs)
    return f_new

end



function StatsModels.coefnames(term::MixedModels.ZeroCorr)
    coefnames(term.term)
end

function lmm_combine_modelmatrix!(Xcomb, X1, X2)
    # we have random effects                
    # combine REMats in single-eventtpe formulas ala y ~ (1|x) + (a|x)
    modelmatrix1 = MixedModels._amalgamate([X1.modelmatrix[2:end]...], Float64)
    modelmatrix2 = MixedModels._amalgamate([X2.modelmatrix[2:end]...], Float64)

    Xcomb = (Xcomb, modelmatrix1..., modelmatrix2...)

    # Next we make the ranefs all equal size
    equalize_ReMat_lengths!(Xcomb[2:end])

    # check if ranefs can be amalgamated. If this fails, then MixedModels tried to amalgamate over different eventtypes and we should throw the warning
    # if it success, we have to check if the size before and after is identical. If it is not, it tried to amalgamize over different eventtypes which were of the same length

    try
        reterms = MixedModels._amalgamate([Xcomb[2:end]...], Float64)
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

    return Xcomb
end


function change_ReMat_size!(remat::MixedModels.AbstractReMat, m::Integer)

    n = m - length(remat.refs) # missing elements
    if n < 0
        deleteat!(remat.refs, range(m + 1, stop = length(remat.refs)))
        remat.adjA = remat.adjA[:, 1:(m+1)]
        remat.wtz = remat.wtz[:, 1:(m+1)]
        remat.z = remat.z[:, 1:(m+1)]
    elseif n > 0
        append!(remat.refs, repeat([remat.refs[end]], n))
        remat.adjA = [remat.adjA sparse(repeat([0], size(remat.adjA)[1], n))]
        remat.wtz = [remat.wtz zeros((size(remat.wtz)[1], n))]
        remat.z = [remat.z zeros((size(remat.z)[1], n))]
    end

    #ReMat{T,S}(remat.trm, refs, levels, cnames, z, wtz, λ, inds, adjA, scratch)
end


"""
$(SIGNATURES)
Get the timeranges where the random grouping variable was applied
"""
function get_timeexpanded_random_grouping(tblGroup, tblLatencies, basisfunction)
    ranges = Unfold.get_timeexpanded_time_range.(tblLatencies, Ref(basisfunction))
end



function equalize_ReMat_lengths!(remats::NTuple{A,MixedModels.AbstractReMat}) where {A}
    # find max length
    m = maximum([x[1] for x in size.(remats)])
    @debug print("combining lengths: $m")
    # for each reMat
    for k in range(1, length = length(remats))
        remat = remats[k]
        if size(remat)[1] == m
            continue
        end
        # prolong if necessary
        change_ReMat_size!(remat, m)

    end
end


mutable struct SparseReMat{T,S} <: MixedModels.AbstractReMat{T}
    trm::Any
    refs::Vector{Int32}
    levels::Any
    cnames::Vector{String}
    z::SparseMatrixCSC{T}
    wtz::SparseMatrixCSC{T}
    λ::LowerTriangular{T,Matrix{T}}
    inds::Vector{Int}
    adjA::SparseMatrixCSC{T,Int32}
    scratch::Matrix{T}
end

"""
$(SIGNATURES)
This function timeexpands the random effects and generates a ReMat object
"""
function StatsModels.modelcols(
    term::Unfold.TimeExpandedTerm{
        <:Union{<:RandomEffectsTerm,<:MixedModels.AbstractReTerm},
    },
    tbl,
)
    # exchange this to get ZeroCorr to work


    tbl = DataFrame(tbl)
    # get the non-timeexpanded reMat
    reMat = modelcols(term.term, tbl)

    # Timeexpand the designmatrix
    z = transpose(Unfold.time_expand(transpose(reMat.z), term, tbl))

    z = disallowmissing(z) # can't have missing in here


    # First we check if there is overlap in the timeexpanded term. If so, we cannot continue. Later implementations will remedy that
    #println(dump(term,))
    if hasfield(typeof(term.term), :rhs)
        rhs = term.term.rhs


    elseif hasfield(typeof(term.term.term), :rhs)
        # we probably have something like zerocorr, which does not need to show a .rhs necessarily
        rhs = term.term.term.rhs
    else

        println("term.term: $(dump(term.term))")
        error("unknown RE structure, has no field .rhs:$(typeof(term.term))")

    end

    group = tbl[!, rhs.sym]
    time = tbl[!, term.eventfields[1]]

    # get the from-to onsets of the grouping varibales
    onsets = get_timeexpanded_random_grouping(group, time, term.basisfunction)
    #print(size(reMat.z))
    refs = zeros(size(z)[2]) .+ 1
    for (i, o) in enumerate(onsets[2:end])
        # check for overlap
        if (minimum(o) <= maximum(onsets[i+1])) & (maximum(o) <= minimum(onsets[i+1]))
            error("overlap in random effects structure detected, not currently supported")
        end
    end

    # From now we can assume no overlap
    # We want to fnd which subject is active when
    refs = zeros(size(z)[2]) .+ 1
    uGroup = unique(group)

    for (i, g) in enumerate(uGroup[1:end])

        ix_start = findfirst(g .== group)
        ix_end = findlast(g .== group)
        if i == 1
            time_start = 1
        else
            time_start = time[ix_start]
            # XXX Replace this functionality by shiftOnset?
            time_start = time_start - sum(Unfold.times(term.basisfunction) .<= 0)
        end
        if i == length(uGroup)
            time_stop = size(refs, 1)
        else
            time_stop = time[ix_end]
            time_stop = time_stop + sum(Unfold.times(term.basisfunction) .> 0)
        end
        if time_start < 0
            time_start = 1
        end

        if time_stop > size(refs, 1)
            time_stop = size(refs, 1)
        end


        #println("$g,$time_start,$time_stop")
        refs[Int64(time_start):Int64(time_stop)] .= i
    end

    # Other variables with implementaions taken from the LinerMixedModel function
    wtz = z
    trm = term

    S = size(z, 1)
    T = eltype(z)
    λ = LowerTriangular(Matrix{T}(I, S, S))

    inds = MixedModels.sizehint!(Int[], (S * (S + 1)) >> 1) # double trouble? can this line go?
    m = reshape(1:abs2(S), (S, S))
    inds = sizehint!(Int[], (S * (S + 1)) >> 1)
    for j = 1:S
        for i = j:S
            # We currently restrict to diagonal entries
            # Once mixedmodels#293 is pushed, we can relax this and use zerocorr()
            if !(typeof(term.term) <: MixedModels.ZeroCorr) || (i == j) # for diagonal
                push!(inds, m[i, j])
            end
        end
    end

    levels = reMat.levels
    refs = refs


    # reMat.levels doesnt change
    cnames = coefnames(term)
    #print(refs)
    adjA = MixedModels.adjA(refs, z)
    scratch = Matrix{T}(undef, (S, length(uGroup)))

    # Once MixedModels.jl supports it, can be replaced with: SparseReMat
    ReMat{T,S}(rhs, refs, levels, cnames, z, wtz, λ, inds, adjA, scratch)
end
function change_modelmatrix_size!(m, fe::AbstractSparseMatrix, remats)
    change_ReMat_size!.(remats, Ref(m))
    @debug "changemodelmatrix" typeof(fe)
    fe = SparseMatrixCSC(m, fe.n, fe.colptr, fe.rowval, fe.nzval)
    return (fe, remats...)
end
function change_modelmatrix_size!(m, fe::AbstractMatrix, remats)
    change_ReMat_size!.(remats, Ref(m))
    fe = fe[1:m, :]
    return (fe, remats...)
end



#modelcols(rhs::MatrixTerm, tbl) = modelcols.(rhs, Ref(tbl))