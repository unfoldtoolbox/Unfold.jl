# implement different basisfunctions e.g.


struct BasisFunction
    kernel::Function # function that given an event onset, results in a kernel used to timeexpand
    colnames::AbstractVector # vector of names along columns of kernel-output
    times::AbstractVector # vector of times along rows of kernel-output (in seconds)
    type::String # type of basisfunction (bookkeeping)
    name::String # name of the event, XXX used for naming coefficients
    shiftOnset::Integer # by how many samples do we need to shift the event onsets?
end



function Base.show(io::IO, obj::BasisFunction)
    println(io, "name: $(obj.name)")
    println(io, "colnames: $(obj.colnames)")
    println(io, "kernel: $(obj.type)")
end

firbasis(;τ,sfreq)           = firbasis(τ,sfreq,"")
firbasis(;τ,sfreq,name="basis_"*string(rand(1:10000))) = firbasis(τ,sfreq,name)
firbasis(τ,sfreq)            = firbasis(τ,sfreq,"")

function firbasis(τ,sfreq,name::String)

    times =range(τ[1],stop=τ[2]+ 1 ./sfreq,step=1 ./sfreq)
    kernel=e->firkernel(e,times[1:end-1])
    type = "firkernel"

    shiftOnset = Int64(round(τ[1] * sfreq))

    return BasisFunction(kernel,times[1:end-1],times,type,name,shiftOnset)
end


function firkernel(e,times)
    @assert ndims(e) <= 1 #either single onset or a row vector where we will take the first one :)
    if size(e,1) > 1
        # XXX we will soon assume that the second entry would be the duration
        e = Float64(e[1])
    end
    e = [1 .- (e .% 1)  e .% 1]
    e[isapprox.(e,0,atol = 1e-15)] .= 0
    ksize=length(times) # kernelsize

    kernel = spdiagm(ksize+1,ksize,0 => repeat([e[1]],ksize), -1 => repeat([e[2]],ksize))

    return(kernel)

end

function hrfbasis(TR::Float64;parameters= [6. 16. 1. 1. 6. 0. 32.],name::String="basis_"*string(rand(1:10000)))
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
    kernel=e->hrfkernel(e,TR,parameters)
    type = "hrfkernel"
    return BasisFunction(kernel,["hrf"],range(0,(length(kernel([0, 1]))-1)*TR,step=TR),type,name,0)
end

function hrfkernel(e,TR,p)
    @assert ndims(e) <= 1 #either single onset or a row vector where we will take the first one :)

    # code adapted from SPM12b
    mt = 32
    dt  = TR/mt;
    firstElem = Int(round((1-(e[1]%1)).*mt))

    if size(e,1) == 2 && e[2]>0.
        duration = e[2]
    else
        duration = 1. /mt
    end
    duration_mt = Int(ceil(duration .* mt))
    box_mt = [zeros(firstElem)' ones(duration_mt)']'
    #box_mt = box_mt ./ (sum(box_mt))
    # u   = [0:ceil(p(7)/dt)] - p(6)/dt;
    u   = range(0,stop = ceil((p[7]+duration)/dt)) .- p[6]/dt;

    #hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3))
    # Note the inverted scale parameter compared to SPM12.
    g1  = Gamma(p[1] ./ p[3],p[3] ./ dt)
     #spm_Gpdf(u,p(2)/p(4),dt/p(4))
    g2  = Gamma(p[2] ./ p[4],p[4] ./ dt);
    # g1 - g2/p(5);
    hrf = pdf.(g1,u) .- pdf.(g2,u) ./ p[5]

    # hrf = hrf([0:floor(p(7)/RT)]*fMRI_T + 1);
    hrf = hrf ./ sum(hrf);
    hrf = conv(box_mt,hrf)'
    #println(box_mt)
    hrf = hrf[((range(0,stop=Int(floor(p[7] ./ TR + duration)))*mt)) .+ 1];



    return(SparseMatrixCSC(sparse(hrf)))

end
