using DifferentialEquations
using Plots
using Catalyst
using DataFrames
using CSV
using Statistics
using Peaks
using FFTW

#Known parameters
const Km1exp = 15.0
const ka2exp = 0.7*10^-3
const kb2exp = 2.0 * 10^-3
const ka3exp = 0.118
const kb3exp = 0.0609
const Kd4exp =29.0
const Km7exp = 39.0
const kcat7exp = 85.3

#Estimates for unknown parameters
ka1est = 0.9 * 10^-3 #Slightly less than ka2
kb1est = 0.2 * 10^-3
ka4est = 0.12 #Estimating LpA bindking to K (ka3) similar to LpA binding to L
ka7est = 0.9 #Pretty much guess, around ka1
dfest = 750

#parameter list Changed this around
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
        :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
        :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
p = [x[2] for x in psym]

    
#initial condition list
usym = [:L => 5, :K => 0.1, :P => 0.1, :A => 0.1, :Lp => 0.1, :LpA => 0.0, :LK => 0, 
        :LpP => 0, :LpAK => 0., :LpAP => 0., :LpAKL => 0., :LpAPLp => 0, :AK => 0, :AP => 0, 
        :AKL => 0, :APLp => 0]
u0 = [x[2] for x in usym]

const fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) K(t) P(t) A(t) Lp(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
    # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
    # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    (ka1,kb1), L + K <--> LK # L binding to kinase
    kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
    kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

    #previously excluded reactions, all possible combinations possible in vitro
    (ka2,kb2), Lp + AK <--> LpAK
    (ka2*y,kb2), Lp + AKL <--> LpAKL
    (ka2,kb2), Lp + AP <--> LpAP
    (ka2*y,kb2), Lp + APLp <--> LpAPLp
    (ka3,kb3), A + K <--> AK
    (ka4,kb4), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3*y,kb3), LpA + LK <--> LpAKL
    (ka4*y,kb4), LpA + LpP <--> LpAPLp
    (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end

const tspan = (0., 100.)

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

dfRange = 10.0 .^ (-1:3) #exponential range from 0.1 to 1000
kaRange = 10.0 .^ (-3:1)
kbRange = 10.0 .^ (-3:3)

Lrange = 10.0 .^ (-2:2)
Krange = 10.0 .^ (-3:2)
Prange = 10.0 .^ (-3:2)
Arange = 10.0 .^ (-2:2)

u0combos = Array{Float64}(undef, length(Lrange)*length(Krange)*length(Arange)*length(Prange),4)

#Set up array of initial conditions that will be tested in advance for ease of later code
count = 1
for L in Lrange
    for K in Krange
        for P in Prange
            for A in Arange
                u0combos[count, 1] = L
                u0combos[count, 2] = K
                u0combos[count, 3] = P
                u0combos[count, 4] = A
                count+=1
            end
        end
    end
end

#Set up ODE problem
const prob = ODEProblem(fullrn, u0, tspan, p)

oscData = DataFrame(DF=Float64[], ka7=Float64[], ka4=Float64[], ka1=Float64[], kb1=Float64[],
                    fit=Float64[],per=Float64[], amp=Float64[],u0=Vector{Float64}[])
nonoscData = DataFrame(DF=Float64[], ka7=Float64[], ka4=Float64[], ka1=Float64[], kb1=Float64[])

numU0combos = length(u0combos[:,1])

for df in dfRange
    dfest = df
    for ka7 in kaRange
        ka7est = ka7
        for ka4 in kaRange
            ka4est = ka4
            for ka1 in kaRange
                ka1est = ka1
                for kb1 in kbRange
                    kb1est = kb1
                    #Set up parameter values
                    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
                    p = [x[2] for x in psym]
                    #Count how many oscillatory solutions
                    oscFound = 0
                    for i in 1:numU0combos
                        curU0 = u0combos[i,:]
                        u0[1:4] = curU0[1:4]
                        currentSol = solve(remake(prob, u0=u0, p=p), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
                        costPerAmp = eval_fitness_catcherrors(currentSol)
                        if(costPerAmp[1] < 0)
                            oscFound+=1
                            if oscFound >= 3
                                push!(oscData, Dict(:DF => df, :ka7 => ka7, :ka4 => ka4, :ka1 => ka1, :kb1 => kb1,
                                                    :fit => costPerAmp[1], :per => costPerAmp[2], :amp => costPerAmp[3],:u0 => curU0))
                                break
                            end
                        end
                    end
                    if oscFound < 3
                        push!(nonoscData, Dict(:DF => df, :ka7 => ka7, :ka4 => ka4, :ka1 => ka1, :kb1 => kb1))
                    end
                    oscFound = 0
                end
            end
        end
        println("You are at ka7=$(ka7) and df=$(df)")
    end
    oscFileName = "OscCsvs/Osc_df_$(df).csv"
    nonoscFileName = "OscCsvs/No_osc_df$(df).csv"
    CSV.write(oscFileName,oscData)
    CSV.write(nonoscFileName,nonoscData)
end



                            
