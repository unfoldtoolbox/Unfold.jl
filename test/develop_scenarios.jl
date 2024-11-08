##
using Test, StatsModels
using DataFrames
import Logging
using Unfold
include("test_utilities.jl")

Logging.global_logger(Logging.SimpleLogger(stdout, Logging.Debug))


data, evts = loadtestdata("testCase3") #
append!(data, zeros(1000))
data = reshape(data, 1, :)
data = vcat(data, data)
data = data .+ 1 * randn(size(data)) # we have to add minimal noise, else mixed models crashes.
data_missing = Array{Union{Missing,Number}}(undef, size(data))
data_missing .= data
data_missing[4500:4600] .= missing

categorical!(evts, :subject)
f = @formula 0 ~ 1 + condA + condB + (1 + condA + condB | subject)
#f  = @formula 0~1 + (1|subject)



# cut the data into epochs
# TODO This ignores subject bounds
data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
evts_e, data_e = Unfold.drop_missing_epochs(evts, data_e)

##

basisfunction = firbasis(τ = (-0.1, 0.3), sfreq = 10)

f = @formula 0 ~ 1 + condA + condB # 1
Xs = designmatrix(UnfoldLinearModel, f, evts_e)
ufModel_A = unfoldfit(UnfoldLinearModel, Xs, data_e)

basisfunction = firbasis(τ = (-0.1, 0.3), sfreq = 10)
Xs = designmatrix(UnfoldLinearModel, f, evts, basisfunction)
basisfunction = firbasis(τ = (-0.1, 0.5), sfreq = 10)
Xs2 = designmatrix(UnfoldLinearModel, f, evts, basisfunction)
Xs = Xs + Xs2

ufModel_B = unfoldfit(UnfoldLinearModel, Xs, data)


f = @formula 0 ~ 1 + condA + condB + (1 + condA | subject)
Xs = designmatrix(UnfoldLinearMixedModel, f, evts_e)

ufModel_C = unfoldfit(UnfoldLinearMixedModel, Xs, data_e)

basisfunction = firbasis(τ = (-0.1, 0.3), sfreq = 10)
Xs = designmatrix(UnfoldLinearMixedModel, f, evts, basisfunction)

ufModel_D = unfoldfit(UnfoldLinearMixedModel, Xs, data)


ufA = condense_long(ufModel_A, times)
ufB = condense_long(ufModel_B)
ufC = condense_long(ufModel_C, times)
ufD = condense_long(ufModel_D)

plot(ufA.colname_basis, ufA.estimate)
plot(ufC.colname_basis, ufC.estimate)
plot(ufB.colname_basis, ufB.estimate)
plot(ufD.colname_basis, ufD.estimate)
