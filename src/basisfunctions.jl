# implement different basisfunctions e.g.


struct BasisFunction{T<:AbstractVector}
    kernel::Function
    times::T
    type::String
    name::String
end



function Base.show(io::IO, obj::BasisFunction)
    println(io, "name: $(obj.name)")
    println(io, "times: $(obj.times)")
    println(io, "kernel: $(obj.type)")
end

firbasis(;τ,sfreq)           = firbasis(τ,sfreq,"")
firbasis(;τ,sfreq,name="") = firbasis(τ,sfreq,name)
firbasis(τ,sfreq)            = firbasis(τ,sfreq,"")

function firbasis(τ,sfreq,name::String)
    times =range(τ[1],stop=τ[2],step=1 ./sfreq)
    kernel=e->firkernel([1-(e%1),e % 1],times)
    type = "firkernel"
    return BasisFunction(kernel,times,type,name)
end


function firkernel(e,times)
    @assert(length(e)==2,"")
    e[isapprox.(e,0,atol = 1e-15)] .= 0
    ksize=length(times) # kernelsize

    kernel = spdiagm(ksize+1,ksize,0 => repeat([e[1]],ksize), -1 => repeat([e[2]],ksize))

    return(kernel)

end

function hrfbasis(TR::Float64;parameters= [6. 16. 1. 1. 6. 0. 32.],name::String="")
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
    return BasisFunction(kernel,[times],type,name)
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
