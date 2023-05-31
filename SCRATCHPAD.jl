using BenchmarkTools, Profile

testarray = rand(Float64, 1000000)
testpeakidxs, testpeakvals = findmaxima(testarray, 1)

@btime @fastmath getDif($testpeakidxs, $testarray)


#< Benchmarking automatic vs stiff solvers
@benchmark autosol = solve(fullprob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
@benchmark @fastmath solve(fullprob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
plot(autosol)
@benchmark stiffsol = solve(fullprob, Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
plot!(stiffsol)

difference = sum(autosol - stiffsol)



@btime CostFunction($autosol)
@btime @fastmath CostFunction($autosol)


#< PerAmp testing 
per, amp = getPerAmp(autosol)


eval_ic_fitness(fullprob.u0, fullprob)


#< Comparing getDif_bidirectional to getDif
getDif(testpeakidxs, testarray)
getDif(testpeakvals)
getDif_bidirectional(testpeakvals)


#< Testing Supertypes 
abstract type AbstractTestType <: ConstraintType end

struct TestType1 <: AbstractTestType 
    a::Int
    b::Int
    c::Int
end

testype1 = TestType1(1,2,3)

for n in testype1
    println(n)
end


#<Costfunction testing
function old_getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) #get difference between fft peak indexes
    #=
    :param indexes: the indexes of the peaks in the Fourier transform of a solution
    :param arrayData: the normalized absolute values of the rfft of a solution
    =#
    arrLen = lengt h(indexes)
    if arrLen < 2
        return 0.0 #? If there is only one peak, the score is set to 0. May not be necessary
    end
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1))
    sum_diff += arrayData[indexes[end]]
    return sum_diff
end

"""
param peakindxs: the indexes of the peaks in the Fourier transform of a solution
:param arrayData: the normalized absolute values of the rfft of a solution
:window_ratio: percent of frequencies in arrayData around peak that the stdev is calculated for
return sum of stdev of data points around the peak divided by the number of peaks (First term of equation 4 in Pušnik et al)
"""
function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
    window = max(1, round(Int, window_ratio * length(arrayData)))
    sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(length(arrayData), ind + window)]) for ind in peakindxs)
    return sum_std / length(peakindxs)
end

function getFrequencies(y::Vector{Float64})
    res = abs.(rfft(y))
    return res ./ cld(length(y), 2) #normalize amplitudes
    #? That length is defined to be half the length +1 of the input of the data set.
    #? Instead,
end


a=[1.0,3.0,1.0,10.0,1.0,3.0,1.0,3.0,1.0]
ffta=getFrequencies(a)
b=2*a
fftb=getFrequencies(b)
indicesA = findmaxima(ffta,1)[1]
indicesB = findmaxima(fftb,1)[1]
getSTD(indicesA,ffta,0.001)
getSTD(indicesB,fftb,0.001)
old_getDif(indicesA, ffta)
old_getDif(indicesB, fftb)
