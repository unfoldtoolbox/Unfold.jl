import StatsBase.predict

function predict(model::UnfoldLinearModel, newdata)
    # make a copy of it so we don't change it outside the function
     data = copy(newdata)
    #data = newdata
    # if only one formulas
    if typeof(model.X.formulas) <: FormulaTerm
        formulas = [model.X.formulas]
    else
        formulas = model.X.formulas
    end    
    # add those as latencies [xxx change to basisfunction field] to newdata
    ##

    X = []
    fromTo = []
    times = []
    
    for f = formulas
        # find out how long each designmatrix is
        n_range = length(f.rhs.basisfunction.times)
        # find out how much to shift so that X[1,:] is the the first "sample"
        n_negative = f.rhs.basisfunction.shiftOnset
        # generate the correct eventfields (default: latencies)
        data[:,f.rhs.eventfields[1]] = range(-n_negative + 1, step=n_range, length=size(data, 1))
        # get the model
        Xsingle = modelcols(f.rhs, data)
        # combine designmats
        append!(X, [Xsingle])

        # keep track of how long each event is
        append!(fromTo, [range(1, step=n_range, length=1+size(data, 1))])

        # keep track of the times
        append!(times, repeat(f.rhs.basisfunction.times[1:end], size(data, 1)))

    end
    
    # Concat them, but without introducing any overlap, we want the isolated responses
    Xconcat = blockdiag(X...)
    
    # calculate yhat
    yhat =  Xconcat * model.beta'

    
    # init the meta dataframe
    metaData = DataFrame([:times => vcat(times...),:basisname=>""])
    for c = names(newdata)
        metaData[:,c] .= newdata[1,c] # assign first element in order to have same column type
    end
   
    # shift variable to keep track of multiple basisfunctions
    shift = 0
    # for each basis function
    for (bIx,basisfun) = enumerate(fromTo)
        # go through all predictors
        for (i, fstart) in enumerate(basisfun[1:end-1])

            # the last one is skipped in the enumerate above
            fend = basisfun[i + 1]-1
           
            # couldn't figure out how to broadcast everything directly (i.e out[fstart:fend,names(newdata)] .= newdata[i,:])
            # copy the correct metadata
            for j = fstart:fend
                metaData[j+shift,names(newdata)] = newdata[i,:]
            end
            # add basisfunction name
            metaData[shift.+(fstart:fend),:basisname] .= formulas[bIx].rhs.basisfunction.name
        end
        # the next meta data has to be at the end
        shift += basisfun[end]-1

    end


    out = DataFrame([:yhat => vec(reshape(yhat, :, 1))])
    nchannel = size(yhat, 2)
    out.channel = repeat(1:nchannel, inner=size(yhat, 1))
   
    out = hcat(out, repeat(metaData, nchannel))
    return out
end
