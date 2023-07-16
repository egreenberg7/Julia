using DifferentialEquations
using DiffEqCallbacks
using Statistics
using Peaks
using FFTW

function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData::Vector{Float64}; window::Int = 1)#, window_ratio::Float64) #get average standard deviation of fft peak indexes
    arrLen = length(fft_arrayData)
    sum_std = @inbounds sum(std(fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs) #* sum rolling window of standard deviations
    return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std
end

"""Return normalized FFT of solution vector"""
function getFrequencies(timeseries::Vector{Float64}) #todo fix normalization or something 
    res = abs.(rfft(timeseries))
    return res ./ cld(length(timeseries), 2) #* normalize amplitudes
end

"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunction(sol::ODESolution)::Vector{Float64}
    #*get| the fft of the solution
    if isSteady(sol)
        return [1.0, 0.0, 0.0]
    end
    fftData = getFrequencies(sol.u)
    fft_peakindexes, fft_peakvals = findmaxima(fftData,10) #* get the indexes of the peaks in the fft
    time_peakindexes, time_peakvals = findmaxima(sol.u,5) #* get the times of the peaks in the fft
    if length(fft_peakindexes) < 2 || length(time_peakindexes) < 2 #* if there are no peaks in either domain, return 0
        return [1.0, 0.0, 0.0]
    end
    std = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain

    #* Compute the period and amplitude
    period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)

    #* Return cost, period, and amplitude as a vector
    return [-std, period, amplitude]
end

function eval_fitness_catcherrors(sol::ODESolution)
    Y = nothing
    try 
        if sol.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(sol) for x in array) || any(x==1 for array in isless.(sol, 0.0) for x in array)
            return [1.0, 0.0, 0.0]
        end
    catch e 
        if isa(e, DomainError) #catch domain errors
            return [1.0,0.0,0.0]
        else
            rethrow(e) #rethrow other errors
        end
    end
    if !(sol.retcode in (ReturnCode.Success, ReturnCode.Default)) #I added this
        println("FAILURE")
        return [1.0,0.0,0.0]
    end
    return CostFunction(sol)
end


"""
Returns true if an ODE solution has more than 3 local maxima, and the height of the middle peak varies from
the height of the last peak by less than 5%.
- `Y::ODESolution` ODESolution that you are testing
"""
function hasPeaks(Y::ODESolution)
    peaks = findmaxima(Y.u, 10)[2]
    numPeaks = length(peaks)
    if numPeaks > 3 && abs((peaks[div(numPeaks,2)]-peaks[end])) / peaks[div(numPeaks,2)] < 0.05 
        return true
    else
        return false
    end
end

function peaksClassifier(sol::ODESolution)
    if(isSteady(sol) || !hasPeaks(sol))
        return [1.0,0.0,0.0]
    end
    #* Compute the period and amplitude
    time_peakindexes, time_peakvals = findmaxima(sol.u,5) #* get the times of the peaks in the fft
    period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)
    return [-1.0, period, amplitude]
end

function peaksClassifierNoErrors(sol::ODESolution)
    Y = nothing
    try 
        if sol.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(sol) for x in array) || any(x==1 for array in isless.(sol, 0.0) for x in array)
            return [1.0, 0.0, 0.0]
        end
    catch e 
        if isa(e, DomainError) #catch domain errors
            return [1.0,0.0,0.0]
        else
            rethrow(e) #rethrow other errors
        end
    end
    if !(sol.retcode in (ReturnCode.Success, ReturnCode.Default)) #I added this
        #println("FAILURE")
        return [1.0,0.0,0.0]
    end
    return peaksClassifier(sol)
end