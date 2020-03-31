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
firbasis(;τ,sfreq,name="") = firbasis(τ,sfreq,name)
firbasis(τ,sfreq)            = firbasis(τ,sfreq,"")

function firbasis(τ,sfreq,name::String)

    times =range(τ[1],stop=τ[2]+ 1 ./sfreq,step=1 ./sfreq)
    kernel=e->firkernel(e,times)
    type = "firkernel"

    shiftOnset = Int64(τ[1] * sfreq)

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

    kernel = spdiagm(ksize+1,ksize,0 => repeat([e[:,1]],ksize), -1 => repeat([e[:,2]],ksize))

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
    kernel=e->hrfkernel(e,TR,parameters)
    type = "hrfkernel"
    return BasisFunction(kernel,["hrf"],range(0,(length(kernel(0))-1)*TR,step=TR),type,name,0)
end

function hrfkernel(e,TR,p)
    # code adapted from SPM12b
    u   = 1-(e%TR) .+ range(0,stop = ceil(p[7]),step=TR) .- p[6];
    g1  = Gamma(p[1] ./ p[3],1 ./ p[3])
    g2  = Gamma(p[2] ./ p[4],1 ./ p[4]);
    hrf = pdf.(g1,u) .- pdf.(g2,u) ./ p[5]
    hrf = hrf[range(1,stop=Int(1+floor(p[7] ./ TR )))];
    hrf = hrf ./ sum(hrf);

    return(SparseMatrixCSC(sparse(hrf)))

end
