module Evaluator
    using DifferentialEquations
    using FFTW
    using Statistics
    include("CustomPeakFinder.jl")

    export fitness_function


    """Takes in an ODEProblem and returns solution excluding first 10% of tspan"""
    function solve_odeprob(prob::OP, idx=[6, 9, 10, 11, 12, 15, 16]) where OP <: ODEProblem
        #* calculate first 10% of the tspan
        tstart = prob.tspan[2] / 10

        #* solve the ODE and only save the last 90% of the solution
        savepoints = tstart+0.1:0.1:prob.tspan[2]
        solve(prob, Rosenbrock23(), saveat=savepoints, save_idxs=idx, verbose=false, maxiters=1e6)
    end


    function fitness_function(params, inits, odeprob, idxs = [])

        newprob = remake(odeprob, u0=inits, p=params)
        sol = solve_odeprob(newprob, idxs)

        if sol.retcode != ReturnCode.Success
            return [0, 0.0, 0.0]
        end

        #* Calculate Amem by summing all the AP2-bound species on the membrane
        Amem_sol = map(sum, sol.u)

        #* Normalize signal to be relative to total AP2 concentration
        Amem_sol ./= sum(inits[idxs]) #replace 17 with get_index(L)

        #* Find the extrema of the time series
        max_idxs, max_vals, min_idxs, min_vals = findextrema(Amem_sol, min_height=0.1)

        period = 0.0
        amplitude = 0.0
        fitness = 0.0

        if is_oscillatory(Amem_sol, sol.t, max_idxs, min_idxs)
            # period, amplitude = getPerAmp(sol.t, max_idxs, max_vals, min_idxs, min_vals)
            fitness += log10(period)
        end

        fitness += get_fitness!(Amem_sol)

        return [fitness, period, amplitude]
    end



    function get_fitness!(solu::Vector{Float64})

        #* Reuse the same time array to preallocate the fft array
        fftData = @view solu[1:(length(solu) ÷ 4) + 1] 

        #* Get the rfft of the solution and normalize it
        fftData = getFrequencies!(fftData, solu) #|> normalize_time_series!

        #* get the indexes of the peaks in the fft
        fft_peakindexes, fft_peakvals = findmaxpeaks(fftData) 

        if isempty(fft_peakvals)
            return 0.0
        end

        #* get the summed standard deviation of the peaks in frequency domain
        standard_deviation = getSTD(fft_peakindexes, fftData) 

        #* get the summed difference between the first and last peaks in frequency domain
        sum_diff = getDif(fft_peakvals) 
        # sum_diff = getWeightedAvgPeakDiff(fft_peakvals, fft_peakindexes)

        #* add the log of the period to the standard deviation and summed difference to calculate fitness and privelage longer periods
        return standard_deviation + sum_diff
    end

    function is_oscillatory(solu::Vector{Float64}, solt::Vector{Float64}, max_idxs::Vector{Int}, min_idxs::Vector{Int})
        if !is_steadystate(solu, solt) && length(max_idxs) > 1 && length(min_idxs) > 1
            return true
        else
            return false
        end
    end

    function is_steadystate(solu::Vector{Float64}, solt::Vector{Float64})
        tstart = cld(length(solt),10) 

        #* Check if last tenth of the solution array is steady state
        testwindow = solu[end-tstart:end]
        if std(testwindow; mean=mean(testwindow)) < 0.01  
            return true
        else
            return false
        end 
    end

    """
        getFrequencies(timeseries)
    Return the real-valued FFT of a timeseries, will be half the length of the timeseries
    """
    function getFrequencies(timeseries::Vector{Float64}; jump::Int = 2) 
        # rfft_result = rfft(@view timeseries[1:2:end])
        sampled_timeseries = @view timeseries[1:jump:end]
        rfft_result = rfft(sampled_timeseries)
        norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
        abs.(rfft_result) ./ norm_val
    end

    """
        getFrequencies!(fft_array, timeseries)
    Computes the real-valued FFT and returns it in-place to the preallocated fft_array, which is half the length of `timeseries`.
    """
    function getFrequencies!(fft_array, timeseries::Vector{Float64}; jump::Int = 2) 
        rfft_result = rfft(@view timeseries[1:jump:end])
        norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
        fft_array .= abs.(rfft_result) ./ norm_val
    end

    """Get summed average standard deviation of peaks values from the FFT of the solution"""
    function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData; window::Int =1) #get average standard deviation of fft peak indexes
        arrLen = length(fft_arrayData)

        #window = max(1,cld(arrLen,window_ratio)) #* window size is 1% of array length, or 1 if array length is less than 100
        sum_std = sum(std(@view fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs; init=0.0) #* sum rolling window of standard deviations

        return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std, add 1 to avoid divide by zero
    end 

    """Get average difference of the first and last peak values from the FFT of the solution"""
    function getDif(peakvals::Vector{Float64})
        (peakvals[begin] - peakvals[end])/length(peakvals)
    end
end