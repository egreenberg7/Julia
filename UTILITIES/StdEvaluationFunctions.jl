#! Helper functions for cost function ## 

#< START 
"""Get summed average standard deviation of peaks in the frequency domain"""
function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData::Vector{Float64}; window::Int = 1)#, window_ratio::Float64) #get average standard deviation of fft peak indexes
    arrLen = length(fft_arrayData)
    sum_std = @inbounds sum(std(fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs) #* sum rolling window of standard deviations
    return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std
end

"""    
Checks if from time:end of the solution the standard deviation is greater than 0.001 * the mean for the data points.
Default time is 80 seconds.
"""
function isSteady(Y::ODESolution; time::Float64 = 80.0)
    timeIndex = findfirst(x -> x > time, Y.t) #Find index of time greater than 80
    solutionInterval = Y.u[timeIndex:end]
    meanVal = mean(solutionInterval)
    stdVal = std(solutionInterval)
    return stdVal > 0.001 * meanVal
end


"""Return normalized FFT of solution vector"""
function getFrequencies(timeseries::Vector{Float64}) #todo fix normalization or something 
    res = abs.(rfft(timeseries))
    return res ./ cld(length(timeseries), 2) #* normalize amplitudes
end

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution, indx_max::Vector{Int}, vals_max::Vector{Float64})
    #* Find peaks of the minima too 
    indx_min, vals_min = findminima(sol.u, 5)

    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))
    # @inbounds amps = 0.5 .* (vals_max .- vals_min)[1:min(length(indx_max), length(indx_min))]

    return mean(pers), mean(amps)

end

"""
Cost function based on stdev term and not having a steady state solution. Based on idea that just looking
at the height of the first peak seems somewhat arbitrary.
Changed minimum length of fft_peakindexes from 2 to 1
"""
function stdCostFunction(Y::ODESolution)
    #*check for steady state
    if isSteady(Y)
        return [0.0, 0.0, 0.0]
    end   
    #*get the fft of the solution
    fftData = getFrequencies(sol.u)
    fft_peakindexes, fft_peakvals = findmaxima(fftData,10) #* get the indexes of the peaks in the fft
    time_peakindexes, time_peakvals = findmaxima(sol.u,5) #* get the times of the peaks in the fft
    std = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain
    #* Compute the period and amplitude
    period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)
    #* Return cost, period, and amplitude as a vector
    return [-std, period, amplitude]
end


#! EVALUATION FUNCTIONS ## 
"""Evaluate the fitness of an individual with new parameters"""
function std_eval_param_fitness(params::Vector{Float64},  prob::ODEProblem)
    # remake with new parameters
    new_prob = remake(prob, p=params)
    return solve_for_fitness_peramp(new_prob)
end

function std_eval_param_fitness(params::Vector{Float64},  prob::ODEProblem)
    # remake with new parameters
    new_prob = remake(prob, p=params)
    return std_solve_for_fitness_peramp(new_prob)
end

"""Evaluate the fitness of an individual with new initial conditions"""
function eval_ic_fitness(initial_conditions::Vector{Float64}, prob::ODEProblem)
    # remake with new initial conditions
    new_prob = remake(prob, u0=[initial_conditions; zeros(length(prob.u0)-length(initial_conditions))])
    return std_solve_for_fitness_peramp(new_prob)
end

"""Utility function to solve the ODE and return the fitness and period/amplitude"""
function std_solve_for_fitness_peramp(prob)

    sol = solve(prob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)

    return stdCostFunction(sol)
end