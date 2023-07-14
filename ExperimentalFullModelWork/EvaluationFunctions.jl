using DifferentialEquations
using DiffEqCallbacks
using Statistics
using Peaks
using FFTW

"""    
Checks if from time:end of the solution the standard deviation is greater than 0.01 * the mean for the data points.
Default time is 80 seconds. Makes it so fft unneeded on clearly nonoscillatory solutions.
"""
function isSteady(Y::ODESolution; time::Float64 = 80.0)
    timeIndex = findfirst(x -> x > time, Y.t) #Find index of time greater than 80
    solutionInterval = Y.u[timeIndex:end]
    meanVal = mean(solutionInterval)
    stdVal = std(solutionInterval)
    return stdVal > 0.01 * meanVal
end

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

function getPerAmp(sol::ODESolution, indx_max::Vector{Int}, vals_max::Vector{Float64})
    #* Find peaks of the minima too 
    indx_min, vals_min = findminima(sol[1,:], 1)

    if length(indx_max) < 1 || length(indx_min) < 1 #todo need to fix this, either keep or move check into cost function
        return 0.0, 0.0
    end

    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))
    # @inbounds amps = 0.5 .* (vals_max .- vals_min)[1:min(length(indx_max), length(indx_min))]

    return mean(pers), mean(amps)
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

#// TODO Determine tspanInc and maximum number of iterations
#// TODO Determine minimum prominence for peaks
#// TODO Determine what needs to be returned
#TODO Implement isDamped
#TODO Finalize isSteady conditions
#// TODO Refactor code to minimize calls to findmaxima
#// TODO Find different between findmaxima and argmaxima
#TODO Have a fun time :)
"""
Given an ODEProblem and parameters, it will first try to classify the solution as oscillatory or not using `finalClassifier`
by integrating from `0` to `shortSpan`. If unsuccessful, it will integrate again from `0` to `longSpan` and then returns 
the values from `finalClassifier` given on this second run.
- `prob` ODEProblem to solve
- `u0` Initial concentrations
- `shortSpan` time on first integration run
- `longSpan` time on second integration run
- `p` Parameter values
"""
function adaptiveSolve(prob::ODEProblem, u0, shortSpan, longSpan, p)
    num_iters = 1
    sol = solve(remake(prob, u0=u0, tspan=(0.0, shortSpan), p=p), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false, callback=TerminateSteadyState())
    ret1 = finalClassifier(sol, shortSpan)#, u0[1])
    if ret1[1] != 4.0
        return ret1
    else
        sol = solve(remake(prob, u0=u0, tspan=(0.0, longSpan), p=p), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false, callback=TerminateSteadyState())
        return finalClassifier(sol, longSpan)#, u0[1])
    end 
end

"""
Evaluates the solution to a differential equation with the following return values:
- `[-tspan, period, amplitude]` The solution was classified as oscillatory when integrated from `0` to `tspan,` the period and amplitude are also returned
- `[1.0, 0.0, 0.0]` Retcode is a failure 
- `[1.1, 0.0, 0.0]` Retcode is `:Terminated` (presumable becase `TerminateSteadyState` was used from `DiffEqCallbacks`)
- `[1.5, 0.0, 0.0]` A lipid concentration higher than the initial concentration `L0` was found in the solution (mass conservation broken), indicating numerical errors
- `[2.0, 0.0, 0.0]` The solution was marked as steady by the `isSteady` function
- `[3.0, 0.0, 0.0]` The solution was marked as damped oscillations by the `isDamped` function
- `[4.0, 0.0, 0.0]` None of the above applied
"""
function finalClassifier(sol::ODESolution, tspan)#, L0)
    if !SciMLBase.successful_retcode(sol) #sol.retcode in (ReturnCode.Unstable, ReturnCode.InitialFailure, ReturnCode.ConvergenceFailure, ReturnCode.Failure)
        return [1.0,0.0,0.0]
    elseif sol.retcode == ReturnCode.Terminated
        return [1.1,0.0,0.0]
    elseif findfirst(x -> x > sol.u[1], sol.u) !== nothing #If lipids not conserved due to numerical issues
        return [1.5, 0.0, 0.0]
    elseif isSteady(sol; time = tspan * 0.8) #Check if last 20% of solution is steady
        return [2.0, 0.0, 0.0]
    else
        maxindices, maxima = findmaxima(sol.u, 10)
        if length(maxima) > 4 #Must have at least 5 local maxima
            peaks, proms = peakproms(maxindices, sol.u; minprom = sol.u[1] * 0.03)
            if(length(peaks) > 4) #Check threshold for prominence to be valid maxima
                if isDamped(proms)
                    return [3.0,0.0,0.0]
                end
            else
                per, amp = getPerAmp(sol, maxindices, maxima)
                return [-tspan, per, amp]
            end
        end
    end
    return [4.0,0.0,0.0]
end


"""
Determines if a given set of peak prominences represents a damped oscillation. Specifically,
it compares the median prominence with the last prominence times 1.02 and the maximum prominence of the three peaks after the median
- `proms` vector of peak prominences from `peakproms` function call
"""
function isDamped(proms)
    midpoint = div(length(proms),2) 
    return proms[midpoint] > 1.02 * proms[end] && (proms[midpoint] > maximum(proms[(midpoint + 1):(midpoint+3)]))
end