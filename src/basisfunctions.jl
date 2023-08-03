"""
See FIRBasis for an examples

    a BasisFunction should implement:
    kernel() 
    collabel() [default "colname_basis"] # name for 
    colnames() # unique names of expanded columns
    times() # vector of times along expanded columns
    name() # name of basis
    width() # expansion to how many columns

    shiftOnset() [default 0]
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
julia>  b = FIRBasis(kernelfunction,"derivative",["f(x)"],range(0,(length(kernelfunction([0, 1]))-1)*TR,step=TR),"hrf_kernel","basis_A",0)
```
"""
struct FIRBasis <: BasisFunction
    "a design-matrix kernel function used to timeexpand, given a timepoint in 'sample' timeunits"
    kernel::Function

    " vector of times along rows of kernel-output (in seconds)"
    times::AbstractVector
    "name of the event, random 1:1000 if unspecified"
    name::String
    "by how many samples do we need to shift the event onsets? This number is determined by how many 'negative' timepoints the basisfunction defines"
    shiftOnset::Integer
end

collabel(basis::FIRBasis) = :time
colnames(basis::FIRBasis) = basis.times[1:end-1]


struct SplineBasis <: BasisFunction
    "a design-matrix kernel function used to timeexpand, given a timepoint in 'sample' timeunits"
    kernel::Function
    "name of column dimension (e.g 'time' for FIR, 'derivative', for HRF etc.)"

    "vector of names along columns of kernel-output"
    colnames::AbstractVector
    " vector of times along rows of kernel-output (in seconds)"
    times::AbstractVector
    "name of the event, random 1:1000 if unspecified"
    name::String
    "by how many samples do we need to shift the event onsets? This number is determined by how many 'negative' timepoints the basisfunction defines"
    shiftOnset::Integer
end


struct HRFBasis <: BasisFunction
    "a design-matrix kernel function used to timeexpand, given a timepoint in 'sample' timeunits"
    kernel::Function
    "vector of names along columns of kernel-output"
    colnames::AbstractVector
    " vector of times along rows of kernel-output (in seconds)"
    times::AbstractVector
    "name of the event, random 1:1000 if unspecified"
    name::String
end




function Base.show(io::IO, obj::BasisFunction)
    println(io, "name: $(name(obj))")
    println(io, "collabel: $(collabel(obj))")
    println(io, "colnames: $(colnames(obj))")
    println(io, "kerneltype: $(typeof(obj))")
    println(io, "times: $(times(obj))")
    println(io, "shiftOnset: $(shiftOnset(obj))")
end




"""
$(SIGNATURES)
Generate a sparse FIR basis around the *τ* timevector at sampling rate *sfreq*. This is useful if you cannot make any assumptions on the shape of the event responses. If unrounded events are supplied, they are split between samples. E.g. event-latency = 1.2 will result in a "0.8" and a "0.2" entry.

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
function firbasis(τ, sfreq, name::String)
    τ = round_times(τ, sfreq)
    times = range(τ[1], stop = τ[2]+1 ./ sfreq, step = 1 ./ sfreq) # stop + 1 step, because we support fractional event-timings

    kernel = e -> firkernel(e, times[1:end-1])


    shiftOnset = Int64(floor(τ[1] * sfreq))

    return FIRBasis(kernel, times, name, shiftOnset)
end
# cant multiple dispatch on optional arguments
#firbasis(;τ,sfreq)           = firbasis(τ,sfreq)
firbasis(; τ, sfreq, name) = firbasis(τ, sfreq, name)
firbasis(τ, sfreq) = firbasis(τ, sfreq, "basis_" * string(rand(1:10000)))

"""
$(SIGNATURES)
Calculate a sparse firbasis
# Examples

```julia-repl
julia>  f = firkernel(103.3,range(-0.1,step=0.01,stop=0.31))
```
"""
function firkernel(e, times)
    @assert ndims(e) <= 1 #either single onset or a row vector where we will take the first one :)
    if size(e, 1) > 1
        # XXX we will soon assume that the second entry would be the duration
        e = Float64(e[1])
    end
    e = [1 .- (e .% 1) e .% 1]
    e[isapprox.(e, 0, atol = 1e-15)] .= 0
    ksize = length(times) # kernelsize

    kernel =
        spdiagm(ksize + 1, ksize, 0 => repeat([e[1]], ksize), -1 => repeat([e[2]], ksize))
    
    #    dropzeros!(kernel) # we often get implicit 0, especially if the latencies are rounded
    return (kernel)

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
function hrfbasis(
    TR::Float64;
    parameters = [6.0 16.0 1.0 1.0 6.0 0.0 32.0],
    name::String = "basis_" * string(rand(1:10000)),
)
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
    return HRFBasis(kernel,["f(x)"],range(0, (length(kernel([0, 1])) - 1) * TR, step = TR),name)
end

shiftOnset(basis::HRFBasis) = 0

collabel(basis::HRFBasis) = :derivative
collabel(basis::SplineBasis) = :splineTerm

collabel(uf::UnfoldModel) = collabel(formula(uf))
collabel(form::FormulaTerm) = collabel(form.rhs)
collabel(t::Tuple) = collabel(t[1]) # MixedModels has Fixef+ReEf
collabel(term::Array{<:AbstractTerm}) = collabel(term[1].rhs)  # in case of combined formulas
#collabel(form::MatrixTerm) = collabel(form[1].rhs)

# typical defaults

shiftOnset(basis::BasisFunction) = basis.shiftOnset
colnames(basis::BasisFunction) = basis.colnames
kernel(basis::BasisFunction) = basis.kernel

times(basis::BasisFunction) = basis.times
name(basis::BasisFunction) = basis.name

StatsModels.width(basis::BasisFunction) = length(times(basis))

StatsModels.width(basis::FIRBasis) = length(times(basis))


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
    hrf = hrf[((range(0, stop = Int(floor(p[7] ./ TR + duration)))*mt)).+1]



    return (SparseMatrixCSC(sparse(hrf)))

end



# taken from https://codereview.stackexchange.com/questions/284537/implementing-a-1d-convolution-simd-friendly-in-julia
# to replace the DSP.conv function
function conv1D!( vO :: Array{T, 1}, vA :: Array{T, 1}, vB :: Array{T, 1} ) :: Array{T, 1} where {T <: Real}

    lenA = length(vA);
    lenB = length(vB);

    fill!(vO, zero(T));
    for idxB in 1:lenB
        for idxA in 1:lenA
            @inbounds vO[idxA + idxB - 1] += vA[idxA] * vB[idxB];
        end
    end

    return vO;

end

function conv1D( vA :: Array{T, 1}, vB :: Array{T, 1} ) :: Array{T, 1} where {T <: Real}

    lenA = length(vA);
    lenB = length(vB);

    vO = Array{T, 1}(undef, lenA + lenB - 1);

    return conv1D!(vO, vA, vB);

end