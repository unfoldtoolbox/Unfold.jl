"""
See FIRBasis for an examples

a BasisFunction should implement:
- kernel() # kernel(b::BasisFunction,sample) => returns the designmatrix for that event
- height() # number of samples in continuous time, NaN if not defined
- width()  # number of coefficient columns (e.g. HRF 1 to 3, FIR=height(),except if interpolate=true )

- colnames() # unique names of expanded columns
- times() # vector of times along expanded columns, length = height()

- name() # name of basisfunction
- collabel() [default "colname_basis"] # name for coeftable
- shift_onset() [default 0]
"""
abstract type BasisFunction end


"""
Defines a FIRBasisfunction which can be called for each event, defining the time-expanded basis kernel


$(TYPEDEF)
$(TYPEDSIGNATURES)
$(FIELDS)

(tipp: most users would you want to call firbasis, not generate it manually)
# Examples
```julia-repl
julia>  b = FIRBasis(range(0,1,length=10),"basisA",-1)
```
"""
mutable struct FIRBasis <: BasisFunction
    "vector of times along rows of kernel-output (in seconds)"
    times::Vector
    "name of the event, should be the actual eventName in `eventcolumn` of the dataframes later"
    name::Any
    "by how many samples do we need to shift the event onsets? This number is determined by how many 'negative' timepoints the basisfunction defines"
    shift_onset::Int64
    "should we linearly interpolate events not on full samples?"
    interpolate::Bool
    "should we scale kernel to the duration? If yes, with which method"
    scale_duration::Union{Bool,Interpolations.Degree,Interpolations.InterpolationType}
    #FIRBasis(times, name, shift_onset, interpolate, scale_duration::Bool) = new(
    #    times,
    #    name,
    #    shift_onset,
    #    interpolate,
    #    scale_duration ? Interpolations.Linear() : false,
    #)
end
# default dont interpolate
FIRBasis(times, name, shift_onset) = FIRBasis(times, name, shift_onset, false)
FIRBasis(times, name, shift_onset, interpolate) =
    FIRBasis(times, name, shift_onset, interpolate, false)
@deprecate FIRBasis(kernel::Function, times, name, shift_onset) FIRBasis(
    times,
    name,
    shift_onset,
)
collabel(basis::FIRBasis) = :time
colnames(basis::FIRBasis) = basis.interpolate ? basis.times[1:(end-1)] : basis.times[1:end]


mutable struct SplineBasis <: BasisFunction
    kernel::Function

    "vector of names along columns of kernel-output"
    colnames::AbstractVector
    " vector of times along rows of kernel-output (in seconds)"
    times::Vector
    "name of the event, random 1:1000 if unspecified"
    name::String
    "by how many samples do we need to shift the event onsets? This number is determined by how many 'negative' timepoints the basisfunction defines"
    shift_onset::Int64
end


mutable struct HRFBasis <: BasisFunction
    kernel::Function
    "vector of names along columns of kernel-output"
    colnames::AbstractVector
    " vector of times along rows of kernel-output (in seconds)"
    times::Vector
    "name of the event, random 1:1000 if unspecified"
    name::String
end





"""
$(SIGNATURES)
Generate a sparse FIR basis around the *τ* timevector at sampling rate *sfreq*. This is useful if you cannot make any assumptions on the shape of the event responses. If unrounded events are supplied, they are split between samples. E.g. event-latency = 1.2 will result in a "0.8" and a "0.2" entry.


Advanced: second input can be duration in samples - careful: `times(firbasis)` always assumes duration = 1. Therefore,
issues with LMM and predict will appear!

# keyword arguments
`interpolate` (Bool, default false): if true, interpolates events between samples linearly. This results in `predict` functions to return a trailling 0`
`scale_duration` (Union{Bool,Interpolations-Interpolator}, default false):
    if true, scales the response by the fit-kwargs `eventfields` second entry. That is, the FIR becomes a stepfunction instead of a impulse response.
    if Interpolations.interpolator, e.g. `Interpolations.Linear()` - uses the fit-kwargs `eventfields` second entry to stretch the FIR kernel based on `imresize`. This implements Hassall

# Examples
Generate a FIR basis function from -0.1s to 0.3s at 100Hz
```julia-repl
julia>  f = firbasis([-0.1,0.3],100)
```
Evaluate at an event occuring at sample 103.3
```julia-repl
julia>  f(103.3)
```

"""
function firbasis(
    τ,
    sfreq,
    name = "";
    interpolate = false,
    #max_duration = nothing,
    scale_duration = false,
)
    τ = round_times(τ, sfreq)
    if interpolate
        # stop + 1 step, because we support fractional event-timings
        τ = (τ[1], τ[2] + 1 ./ sfreq)
    end
    #if !isnothing(max_duration)
    #    τ = (τ[1], max_duration)
    #τ = (τ[1], τ[2] + max_height./sfreq)

    #end
    times = range(τ[1], stop = τ[2], step = 1 ./ sfreq)

    shift_onset = Int64(floor(τ[1] * sfreq))

    return FIRBasis(times, name, shift_onset, interpolate, scale_duration)

end
# cant multiple dispatch on optional arguments
#firbasis(;τ,sfreq)           = firbasis(τ,sfreq)
firbasis(; τ, sfreq, name = "", kwargs...) = firbasis(τ, sfreq, name; kwargs...)



"""
$(SIGNATURES)
Calculate a sparse firbasis

second input can be duration in samples - careful: `times(firbasis)` always assumes duration = 1. Therefore,
issues with LMM and predict will appear!

# Examples

```julia-repl
julia>  f = firkernel(103.3,range(-0.1,step=0.01,stop=0.31))
julia>  f_dur = firkernel([103.3 4],range(-0.1,step=0.01,stop=0.31))
```
"""
function firkernel(ev, times; interpolate = false, scale_duration = false)
    @assert ndims(ev) <= 1 #either single onset or a row vector where we will take the first one :)

    # duration is 1 for FIR and duration is duration (in samples!) if e is vector
    e = interpolate ? ev[1] : 1
    dur = (scale_duration == true && size(ev, 1) > 1) ? Int.(ceil(ev[2])) : 1

    #scale_duration == false => 1
    #scale_duration => true
    #dur = (scale_duration == true || (scale_duration != false || size(ev, 1) > 1) ? Int.(ceil(ev[2])) : 1


    ksize = interpolate ? length(times) - 1 : length(times) #isnothing(maxsize) ? length(times) : max_height # kernelsize

    # if interpolatethe first and last entry is split, depending on whether we have "half" events
    eboth = interpolate ? [1 .- (e .% 1) e .% 1] : [1]
    #@show eboth dur
    values = [eboth[1]; repeat([e], dur - 1); eboth[2:end]]
    values[isapprox.(values, 0, atol = 1e-15)] .= 0 # keep sparsity pattern

    # without interpolation, remove last entry again, which should be 0 anyway
    # values = interpolate ? values : values[1:end-1]  # commented out if the eboth[2:end] trick works

    # build the matrix, we define it by diagonals which go from 0, -1, -2 ...
    pairs = [x => repeat([y], ksize) for (x, y) in zip(.-range(0, dur), values)]
    #return pairs, ksize
    #@debug pairs
    x_single = spdiagm(ksize + length(pairs) - 1, ksize, pairs...)
    if scale_duration == false || scale_duration == true
        return x_single
    else
        #@show "imresize"
        return imresize(x_single, ratio = (ev[2] / ksize, 1), method = scale_duration)
    end

end





"""
$(SIGNATURES)
Generate a Hemodynamic-Response-Functio (HRF) basis with inverse-samplingrate "TR" (=1/FS)


Optional Parameters p:
                                                           defaults
                                                          {seconds}
        p(1) - delay of response (relative to onset)          6
        p(2) - delay of undershoot (relative to onset)       16
        p(3) - dispersion of response                         1
        p(4) - dispersion of undershoot                       1
        p(5) - ratio of response to undershoot                6
        p(6) - onset {seconds}                                0
        p(7) - length of kernel {seconds}                    32

# Examples
Generate a HRF basis function object with Sampling rate 1/TR. And evaluate it at an event occuring at TR 103.3 with duration of 4.1 TRs
```julia-repl
julia>  f = hrfbasis(2.3)
julia>  f(103.3,4.1)

```


"""
function hrfbasis(TR::Float64; parameters = [6.0 16.0 1.0 1.0 6.0 0.0 32.0], name = "")
    # Haemodynamic response function adapted from SPM12b "spm_hrf.m"
    # Parameters:
    #                                                           defaults
    #                                                          {seconds}
    #        p(1) - delay of response (relative to onset)          6
    #        p(2) - delay of undershoot (relative to onset)       16
    #        p(3) - dispersion of response                         1
    #        p(4) - dispersion of undershoot                       1
    #        p(5) - ratio of response to undershoot                6
    #        p(6) - onset {seconds}                                0
    #        p(7) - length of kernel {seconds}                    32
    kernel = e -> hrfkernel(e, TR, parameters)
    return HRFBasis(
        kernel,
        ["f(x)"],
        range(0, (length(kernel([0, 1])) - 1) * TR, step = TR),
        name,
    )
end

shift_onset(basis::HRFBasis) = 0

collabel(basis::HRFBasis) = :derivative
collabel(basis::SplineBasis) = :splineTerm

collabel(uf::UnfoldModel) = collabel(formulas(uf))
collabel(form::FormulaTerm) = collabel(form.rhs)

collabel(term::Array{<:AbstractTerm}) = collabel(term[1].rhs)  # in case of combined formulas
#collabel(form::MatrixTerm) = collabel(form[1].rhs)

# typical defaults

shift_onset(basis::BasisFunction) = basis.shift_onset
colnames(basis::BasisFunction) = basis.colnames
kernel(basis::BasisFunction, e) = basis.kernel(e)
@deprecate kernel(basis::BasisFunction) basis.kernel


basisname(fs::Vector{<:FormulaTerm}) = [name(f.rhs.basisfunction) for f in fs]
basisname(uf::UnfoldModel) = basisname(formulas(uf))
basisname(uf::UnfoldLinearModel) = first.(design(uf)) # for linear models we dont save it in the formula

kernel(basis::FIRBasis, e) = firkernel(
    e,
    basis.times;
    interpolate = basis.interpolate,
    scale_duration = basis.scale_duration,
)


times(basis::BasisFunction) = basis.times
name(basis::BasisFunction) = basis.name

StatsModels.width(basis::BasisFunction) = height(basis)
function StatsModels.width(basis::FIRBasis)
    if basis.scale_duration == false#isa(basis.scale_duration, Bool)
        if basis.interpolate
            return height(basis) - 1
        else
            return height(basis)
        end
    else
        return length(times(basis))
    end
end
height(basis::BasisFunction) = length(times(basis))

height(basis::FIRBasis) = isa(basis.scale_duration, Bool) ? length(times(basis)) : 0

StatsModels.width(basis::HRFBasis) = 1
times(basis::HRFBasis) = NaN # I guess this could also return 1:32 or something?

"""
$(SIGNATURES)
Calculate a HRF kernel. Input e can be [onset duration]
# Examples

```julia-repl
julia>  f = hrfkernel(103.3,2.3,[6. 16. 1. 1. 6. 0. 32.])
```
"""
function hrfkernel(e, TR, p)
    @assert ndims(e) <= 1 #either single onset or a row vector where we will take the first one :)

    # code adapted from SPM12b
    mt = 32
    dt = TR / mt
    firstElem = Int(round((1 - (e[1] % 1)) .* mt))

    if size(e, 1) == 2 && e[2] > 0.0
        duration = e[2]
    else
        duration = 1.0 / mt
    end
    duration_mt = Int(ceil(duration .* mt))
    box_mt = [zeros(firstElem)' ones(duration_mt)']'
    #box_mt = box_mt ./ (sum(box_mt))
    # u   = [0:ceil(p(7)/dt)] - p(6)/dt;
    u = range(0, stop = ceil((p[7] + duration) / dt)) .- p[6] / dt

    #hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3))
    # Note the inverted scale parameter compared to SPM12.
    g1 = Gamma(p[1] ./ p[3], p[3] ./ dt)
    #spm_Gpdf(u,p(2)/p(4),dt/p(4))
    g2 = Gamma(p[2] ./ p[4], p[4] ./ dt)
    # g1 - g2/p(5);
    hrf = pdf.(g1, u) .- pdf.(g2, u) ./ p[5]

    # hrf = hrf([0:floor(p(7)/RT)]*fMRI_T + 1);
    hrf = hrf ./ sum(hrf)
    hrf = conv1D(box_mt, hrf)'
    #println(box_mt)
    hrf = hrf[((range(0, stop = Int(floor(p[7] ./ TR+duration)))*mt)) .+ 1]



    return (SparseMatrixCSC(sparse(hrf)))

end



# taken from https://codereview.stackexchange.com/questions/284537/implementing-a-1d-convolution-simd-friendly-in-julia
# to replace the DSP.conv function
function conv1D!(vO::Array{T,1}, vA::Array{T,1}, vB::Array{T,1})::Array{T,1} where {T<:Real}

    lenA = length(vA)
    lenB = length(vB)

    fill!(vO, zero(T))
    for idxB = 1:lenB
        for idxA = 1:lenA
            @inbounds vO[idxA+idxB-1] += vA[idxA] * vB[idxB]
        end
    end

    return vO

end

function conv1D(vA::Array{T,1}, vB::Array{T,1})::Array{T,1} where {T<:Real}

    lenA = length(vA)
    lenB = length(vB)

    vO = Array{T,1}(undef, lenA + lenB - 1)

    return conv1D!(vO, vA, vB)

end
