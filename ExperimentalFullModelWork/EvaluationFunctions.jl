using DifferentialEquations
#using DiffEqCallbacks
using Statistics
using Peaks

"""    
Checks if for last 20% of the solution the standard deviation is less than 0.01 * the mean for the data points.
This eliminates damped oscillations, oscillations with unobservably small amplitudes, and other solutions approaching
a steady state.
"""
function isSteady(Y::ODESolution)
    solutionInterval = Y.u[4 * length(Y.u) ÷ 5:end]
    meanVal = mean(solutionInterval)
    stdVal = std(solutionInterval)
    return stdVal < 0.01 * meanVal
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


#// TODO Determine tspanInc and maximum number of iterations
#// TODO Determine minimum prominence for peaks
#// TODO Determine what needs to be returned
#TODO Implement isDamped
#// TODO Finalize isSteady conditions
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

Return codes from `finalClassifier` are as follows:
- `[-tspan, period, amplitude]` The solution was classified as oscillatory when successfully integrated from `0` to `tspan,` the period and amplitude are also returned
- `[-1.0, period, amplitude]` The solution terminated due to reaching `maxiters` and was classified as oscillatory
- `[1.0, 0.0, 0.0]` Retcode is a failure 
- `[1.1, 0.0, 0.0]` Retcode is `:Terminated` (presumable becase `TerminateSteadyState` was used from `DiffEqCallbacks`)
- `[1.5, 0.0, 0.0]` A lipid concentration higher than the initial concentration `L0` was found in the solution (mass conservation broken), indicating numerical errors
- `[2.0, 0.0, 0.0]` The solution was marked as steady by the `isSteady` function
- `[3.0, 0.0, 0.0]` The solution was marked as damped oscillations by the `isDamped` function
- `[4.0, 0.0, 0.0]` None of the above applied
"""
function adaptiveSolve(prob::ODEProblem, u0, shortSpan, longSpan, p; abstol = 1e-6, reltol=1e-8)
    sol = solve(remake(prob, u0=u0, tspan=(0.0, shortSpan), p=p), Rodas4(), abstol=abstol, reltol=reltol, saveat=0.1, save_idxs=1, maxiters=150.0 * shortSpan, verbose=false)#, callback=TerminateSteadyState(1e-8, 1e-12)
    ret1 = finalClassifier(sol, shortSpan)
    if ret1[1] != 4.0
        return ret1
    else
        sol = solve(remake(prob, u0=u0, tspan=(0.0, longSpan), p=p), Rodas4(), abstol=abstol, reltol=reltol, saveat=0.1, save_idxs=1, maxiters=150.0 * longSpan, verbose=false)#, callback=TerminateSteadyState(1e-8,1e-12)
        return finalClassifier(sol, longSpan)
    end 
end

"""
Evaluates the solution to a differential equation with the following return values:
- `[-tspan, period, amplitude]` The solution was classified as oscillatory when successfully integrated from `0` to `tspan,` the period and amplitude are also returned
- `[-1.0, period, amplitude]` The solution terminated due to reaching `maxiters` and was classified as oscillatory
- `[1.0, 0.0, 0.0]` Retcode is a failure 
- `[1.1, 0.0, 0.0]` Retcode is `:Terminated` (presumable becase `TerminateSteadyState` was used from `DiffEqCallbacks`)
- `[1.5, 0.0, 0.0]` A lipid concentration higher than the initial concentration `L0` was found in the solution (mass conservation broken), indicating numerical errors
- `[2.0, 0.0, 0.0]` The solution was marked as steady by the `isSteady` function
- `[3.0, 0.0, 0.0]` The solution was marked as damped oscillations by the `isDamped` function
- `[4.0, 0.0, 0.0]` None of the above applied
"""
function finalClassifier(sol::ODESolution, tspan)
    if !(SciMLBase.successful_retcode(sol) || sol.retcode == ReturnCode.MaxIters) #sol.retcode in (ReturnCode.Unstable, ReturnCode.InitialFailure, ReturnCode.ConvergenceFailure, ReturnCode.Failure)
        return [1.0,0.0,0.0]
    #elseif sol.retcode == ReturnCode.Terminated
        #return [1.1,0.0,0.0]
    elseif findfirst(x -> x > sol.u[1], sol.u) !== nothing #If lipids not conserved due to numerical issues
        return [1.5, 0.0, 0.0]
    elseif isSteady(sol) #Check if last 20% of solution is steady
        return [2.0, 0.0, 0.0]
    else
        maxindices, maxima = findmaxima(sol.u, 10)
        if length(maxima) > 4 #Must have at least 5 local maxima
            peaks, proms = peakproms(maxindices, sol.u; minprom = sol.u[1] * 0.03)
            if(length(peaks) > 4) #Check threshold for prominence to be valid maxima
                if isDamped(proms)
                    return [3.0,0.0,0.0]
                else
                    per, amp = getPerAmp(sol, maxindices, maxima)
                    if sol.retcode == ReturnCode.MaxIters
                        return [-1.0, per, amp]
                    else
                        return [-tspan, per, amp]
                    end
                end
            end
        end
    end
    return [4.0,0.0,0.0]
end


"""
Determines if a given set of peak prominences represents a damped oscillation. Specifically,
it compares the median prominence with the last prominence times 1.25
- `proms` vector of peak prominences from `peakproms` function call
"""
function isDamped(proms)
    midpoint = length(proms) ÷ 2 
    return proms[midpoint] > 1.25 * proms[end]# && (proms[midpoint] > maximum(proms[(midpoint + 1):(midpoint+3)]))
end