import StatsBase.predict



function StatsBase.predict(model::UnfoldModel, newdata, times = nothing)
    # make a copy of it so we don't change it outside the function
    data = copy(newdata)
    #data = newdata
    # if only one formulas
    if typeof(formula(model)) <: FormulaTerm
        formulas = [formula(model)]
    else
        formulas = formula(model)
    end
    # add those as latencies [xxx change to basisfunction field] to newdata
    ##

    if !(typeof(formulas[1].rhs) <: Unfold.TimeExpandedTerm)
        # mass univariate model. Not sure how I can multiple dispatch correctly :| Maybe I need 4 types after all

        # just a single designmat, same for all timepoints
        #X = designmatrix(UnfoldLinearModel, formulas[1], data)
        X = modelcols(Unfold.formula(model).rhs,data)
        
        # setup the output matrix, has to be a matrix
        yhat = Array{Union{Missing,Float64}}(
            missing,
            size(coef(model), 1),
            size(X, 1),
            size(coef(model), 2),
        )


        for ch = 1:size(coef(model), 1)
            yhat[ch, :, :] = X * permutedims(coef(model)[ch, :, :], (2, 1))
        end

        # fake times because we never save it in the object
        timesVec = repeat(range(1, length = size(yhat, 3)), size(data, 1))
        if !isnothing(times)
            timesVec = times[timesVec]
        end
        # but times could be optionally given

        # reshape yhat

        yhat = reshape(permutedims(yhat, (1, 3, 2)), size(yhat, 1), :)
        yhat = permutedims(yhat, (2, 1))
    else

        X = []
        fromTo = []
        timesVec = []
        for f in formulas

            # find out how long each designmatrix is
            n_range = length(f.rhs.basisfunction.times)
            # find out how much to shift so that X[1,:] is the the first "sample"
            n_negative = f.rhs.basisfunction.shiftOnset
            # generate the correct eventfields (default: latencies)
            data[:, f.rhs.eventfields[1]] =
                range(-n_negative + 1, step = n_range, length = size(data, 1))



            # get the model
            Xsingle = modelcols(f.rhs, data)

            # remove the last time-point because it is attached due to non-integer latency/eventonsets.
            # e.g. x denotes a sample
            # x- - -x- - -x- - -x
            # e- - - -f- - - -g-
            # 
            # e is aligned (integer) with the sample
            # f&g are between two samples, thus the design matrix would interpolate between them. Thus has as a result, that the designmatrix is +1 longer than what would naively be expected
            #
            # because in "predict" we define where samples onset, we can remove the last sample, it s always 0 anyway, but to be sure we test it
            @assert Xsingle[end, end] == 0.0
            Xsingle = Xsingle[1:end-1, :]

            # combine designmats
            append!(X, [Xsingle])

            # keep track of how long each event is
            append!(fromTo, [range(1, step = n_range, length = size(data, 1))])

            # keep track of the times
            append!(timesVec, repeat(f.rhs.basisfunction.times[1:end], size(data, 1)))

        end

        # Concat them, but without introducing any overlap, we want the isolated responses
        Xconcat = blockdiag(X...)

        # calculate yhat
        # yhat x channels
        yhat = Xconcat * coef(model)'
    end

    # init the meta dataframe
    metaData = DataFrame([:times => vcat(timesVec...), :basisname => ""])
    for c in names(newdata)
        metaData[:, c] .= newdata[1, c] # assign first element in order to have same column type
    end

    if !(typeof(formulas[1].rhs) <: Unfold.TimeExpandedTerm)
        # for mass univariate we can make use of the knowledge that all events have the same length :)
        ntimes = size(coef(model), 2)
        for c in names(newdata)
            for row = 1:size(data, 1)
                rowIx = (1.0 .+ (row - 1) .* ntimes) .+ range(1.0, length = ntimes) .- 1
                metaData[Int64.(rowIx), c] .= data[row, c]
            end
        end

    else


        # shift variable to keep track of multiple basisfunctions
        shift = 0
        # for each basis function
        for (bIx, basisfun) in enumerate(fromTo)
            # go through all predictors
            for (i, fstart) in enumerate(basisfun[1:end-1])

                # the last one is skipped in the enumerate above
                fend = basisfun[i+1] - 1

                # couldn't figure out how to broadcast everything directly (i.e out[fstart:fend,names(newdata)] .= newdata[i,:])
                # copy the correct metadata
                for j = fstart:fend
                    metaData[j+shift, names(newdata)] = newdata[i, :]
                end
                # add basisfunction name
                metaData[shift.+(fstart:fend), :basisname] .=
                    formulas[bIx].rhs.basisfunction.name
            end
            # the next meta data has to be at the end
            shift += basisfun[end] - 1

        end
    end


    out = DataFrame([:yhat => vec(reshape(yhat, :, 1))])
    nchannel = size(yhat, 2)

    out.channel = repeat(1:nchannel, inner = size(yhat, 1))

    out = hcat(out, repeat(metaData, nchannel))
    return out
end
