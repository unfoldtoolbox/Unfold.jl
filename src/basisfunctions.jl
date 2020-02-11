# implement different basisfunctions e.g.

struct BasisFunction{T<:AbstractVector}
    kernel::Function
    times::T
end


function Base.show(io::IO, obj::BasisFunction)
    print(io, "times [")
    print(io,obj.times)
    print(io,"]\n")
    print(io,"kernel: ")
    println(io,obj.kernel)
end

function firbasis(;τ,sfreq)
    firbasis(τ,sfreq)
end

function firbasis(τ,sfreq)
    times =range(τ[1],stop=τ[2],step=1 ./sfreq)
    kernel=e->firkernel([1-(e%1),e % 1],times)
    return(BasisFunction(kernel,times))
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
function boldbasis(lambda=6.::Float64,TR=2,τ=(-.5,15))
    error("Not implemented")
end
