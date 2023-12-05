import StatsBase.predict

using DocStringExtensions: format


function StatsBase.predict(model::UnfoldModel, events)
    # make a copy of it so we don't change it outside the function
    newevents = copy(events)
  
    formulas = formula(model)
    if typeof(formulas) <: FormulaTerm
        formulas = [formulas]
    end
    if typeof(model) == UnfoldLinearModel
        eff = yhat(model,formulas[1],newevents)
        timesVec = gen_timeev(times(model),size(newevents, 1)) 
    else
        fromTo,timesVec,eff = yhat(model,formulas,newevents) 
    end

   
    # init the meta dataframe
    metaData = DataFrame([:time => vcat(timesVec...), :basisname => ""])
    
    for c in names(newevents)
        metaData[:, c] .= newevents[1, c] # assign first element in order to have same column type
    end
    if typeof(model) == UnfoldLinearModel
        # for mass univariate we can make use of the knowledge that all events have the same length :)
        ntimes = size(coef(model), 2)
        for c in names(newevents)
            for row = 1:size(newevents, 1)
                rowIx = (1.0 .+ (row - 1) .* ntimes) .+ range(1.0, length = ntimes) .- 1
                metaData[Int64.(rowIx), c] .= newevents[row, c]
            end
        end

    else


        # shift variable to keep track of multiple basisfunctions
        shift = 0
        # for each basis function
        
        for (bIx, basisfun) in enumerate(fromTo)
        
            # go through all predictors
            for (i, fstart) in enumerate(basisfun[1:end])
                
                fend = basisfun[i] + basisfun.step-1
        
                # couldn't figure out how to broadcast everything directly (i.e out[fstart:fend,names(newevents)] .= newevents[i,:])
                # copy the correct metadata
                
                for j = fstart:fend
                    metaData[j+shift, names(newevents)] = newevents[i, :]
                end
                # add basisfunction name
                metaData[shift.+(fstart:fend), :basisname] .=
                    formulas[bIx].rhs.basisfunction.name
                
            end
            # the next meta data has to be at the end
            shift += basisfun[end]-1+basisfun.step

        end
    end
    

    out = DataFrame([:yhat => vec(reshape(eff, :, 1))])
    nchannel = size(eff, 2)

    out.channel = repeat(1:nchannel, inner = size(eff, 1))
    out = hcat(out, repeat(metaData, nchannel))
    return out
end

# in case the formula is not an array
#yhat(model::UnfoldModel,formulas::AbstractTerm) = yhat(model,[formulas])

# special case if one formula is defined but in an array => multiple events
function yhat(model::Union{<:UnfoldLinearMixedModel,<:UnfoldLinearModel},formulas::AbstractArray,newevents::DataFrame)
    X = modelcols.(formulas,Ref(newevents))
    
    co = coef(model)
    Xsizes = size.(X,Ref(2))
    Xsizes_cumsum = vcat(0,cumsum(Xsizes))
    
    indexes = [(Xsizes_cumsum[ix]+1):Xsizes_cumsum[ix+1] for ix = 1:length(formulas)]
    
    coArray = [co[:,:,ix] for ix in indexes]
    

    return yhat.(coArray,X)
    

end


yhat(model::UnfoldLinearModel,formulas::FormulaTerm,newevents) = yhat(model,formulas.rhs,newevents)
yhat(model::UnfoldLinearModelContinuousTime,formulas::FormulaTerm,newevents) = yhat(model,formulas.rhs,newevents)

function yhat(model::Union{<:UnfoldLinearMixedModel,<:UnfoldLinearModel},formulas::AbstractTerm,newevents)#::AbstractArray{AbstractTerm})
    X = modelcols(formulas,newevents)
    return yhat(model,X)
end


yhat(model::UnfoldLinearModelContinuousTime,formulas::MatrixTerm,events) = yhat(model,formulas.terms,events)

function yhat(model::UnfoldLinearModelContinuousTime,formulas,events)#::AbstractArray{AbstractTerm})
    X = []
    fromTo = []
    timesVec = []
    for f in formulas

        # due to many reasons, there can be several different options of how formula arrives here - not easy to catch via multiple dispatch
        if !(isa(f,TimeExpandedTerm))
            if !(isa(f,MatrixTerm))
            f = f.rhs
            #elseif !(isa(f,Effects.TypicalTerm))
                
            else
                f = f.terms[1]
            end
            
        end
        # find out how long each designmatrix is
        n_range = length(times(f.basisfunction))
        # find out how much to shift so that X[1,:] is the the first "sample"
        n_negative = f.basisfunction.shiftOnset
        # generate the correct eventfields (default: latencies)
        events[:, f.eventfields[1]] =
            range(-n_negative + 1, step = n_range, length = size(events, 1))


            
        # get the model
        Xsingle = modelcols(f, events)
        
        timesSingle = times(f.basisfunction)

        # this pertains only to FIR-models

        # remove the last time-point because it is attached due to non-integer latency/eventonsets.
        # e.g. x denotes a sample
        # x- - -x- - -x- - -x
        # e- - - -f- - - -g-
        # 
        # e is aligned (integer) with the sample
        # f&g are between two samples, thus the design matrix would interpolate between them. Thus has as a result, that the designmatrix is +1 longer than what would naively be expected
        #
        # because in "predict" we define where samples onset, we can remove the last sample, it s always 0 anyway, but to be sure we test it
        
        if typeof(f.basisfunction) <: FIRBasis
            
            keep = ones(size(Xsingle,1))
            keep[range(length(timesSingle), size(Xsingle,1),step=length(timesSingle))] .= 0
            Xsingle = Xsingle[keep.==1, :]
            timesSingle = timesSingle[1:end-1]
            n_range = n_range-1
        end
        # combine designmats
        append!(X, [Xsingle])
        
        # keep track of how long each event is
        append!(fromTo, [range(1, step = n_range, length = size(events, 1))])
        # keep track of the times
        append!(timesVec, repeat(timesSingle, size(events, 1)))

    end

    # Concat them, but without introducing any overlap, we want the isolated responses
    Xconcat = blockdiag(X...)

    # calculate yhat
    eff = yhat(model,Xconcat)

    # output is a bit ugly, but we need the other two vectors as well. Should be refactored at some point - but not right now ;-) #XXX
    return fromTo,timesVec,eff
    
end
function yhat(model::UnfoldLinearModelContinuousTime,X::AbstractArray{T,2};times=nothing) where {T<:Union{Missing, <:Number}}
    
    yhat = X * coef(model)'
    return yhat

end

# kept for backwards compatability
yhat(model::Union{<:UnfoldLinearMixedModel,<:UnfoldLinearModel},X::AbstractArray{T,2};kwargs...) where {T<:Union{Missing, <:Number}} = yhat(coef(model),X;kwargs...)


function yhat(coef::AbstractArray,X::AbstractArray{T,2};kwargs...) where {T<:Union{Missing, <:Number}}
    # function that calculates coef*designmat, but in the ch x times x coef vector
    # setup the output matrix, has to be a matrix
    # then transforms it back to 2D matrix times/coef x ch to be compatible with the timecontinuous format
    yhat = Array{Union{Missing,Float64}}(
        missing,
        size(coef, 1),
        size(X, 1),
        size(coef, 2),
    )
    for ch = 1:size(coef, 1)
        yhat[ch, :, :] = X * permutedims(coef[ch, :, :], (2, 1))
    end
    
    # bring the yhat into a ch x yhat format
    yhat = reshape(permutedims(yhat, (1, 3, 2)), size(yhat, 1), :)
    yhat = permutedims(yhat, (2, 1))
    return yhat
end

# there is only 1 times for Mass Univariate Models possible
times(model::Union{<:UnfoldLinearMixedModel,<:UnfoldLinearModel}) = times(design(model))
times(d::Dict{<:Union{Any,<:AbstractString,Symbol},<:Tuple{<:AbstractTerm,<:AbstractVector}}) = first(values(d))[2]#[k[2] for k in values(d)] # probably going for steprange would be better


function gen_timeev(timesVec,nRows)
    timesVec = repeat(timesVec, nRows)
    return timesVec
end