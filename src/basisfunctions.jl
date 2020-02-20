# implement different basisfunctions e.g.


struct BasisFunction{T<:AbstractVector}
    kernel::Function
    times::T
    type::String
end

function Base.show(io::IO, obj::BasisFunction)
    println(io, "times: $(obj.times)")
    println(io, "kernel: $(obj.type)")
end


function firbasis(;τ,sfreq)
    # Helper function to call firbasis with named arguments
    return firbasis(τ,sfreq)
end

function firbasis(τ,sfreq)
    times =range(τ[1],stop=τ[2],step=1 ./sfreq)
    kernel=e->firkernel([1-(e%1),e % 1],times)
    type = "firkernel"
    return BasisFunction(kernel,times,type)
end


function firkernel(e,times)
    @assert(length(e)==2,"")
    e[isapprox.(e,0,atol = 1e-15)] .= 0
    #shift=sum(times.<0) # shift the whole thing by how many negative
    ksize=length(times) # kernelsize
    kernel = zeros(ksize+1,ksize)
    i = 1
    for lagix = 1:ksize
            kernel[lagix:lagix+1,lagix]=e
    end

    return(kernel)

end

function hrfbasis(TR::Float64;parameters= [6. 16. 1. 1. 6. 0. 32.])
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
    times = 0
    kernel=e->hrfkernel(e,TR,parameters)
    type = "hrfkernel"
    return BasisFunction(kernel,[times],type)
end

function hrfkernel(e,TR,p)
    # code adapted from SPM12b
    u   = 1-(e%TR) .+ range(0,stop = ceil(p[7]),step=TR) .- p[6];
    g1  = Gamma(p[1] ./ p[3],1 ./ p[3])
    g2  = Gamma(p[2] ./ p[4],1 ./ p[4]);
    hrf = pdf.(g1,u) .- pdf.(g2,u) ./ p[5]
    hrf = hrf[range(1,stop=Int(1+floor(p[7] ./ TR )))];
    hrf = hrf ./ sum(hrf);

    return(hrf)

end
