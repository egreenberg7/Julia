#=
-----------
My code to generate 10000 graphs with Fourier transforms of oscillation plots with random initial concentrations
of Lp, K, P, and A varied on a log-uniform scale within experimentally realizable values. All other initial concentrations
were set to 0, and the rate constants were kept constant at values from one of your optimization runs.

I used the cost function (score of less than -0.01) followed by a check of if the last 20 seconds of simulation time the standard
deviation fo the data points varied by more than 0.1% of the mean of the data points. The revised cost function
where getdiff is divided by the number of peaks was used. Additionally, all Fourier data is divided by the height
of the highest peak in order to make the size of the amplitude of oscillations irrelevent; rather, only
the shape of the plot will matter. I also cleaned up some of the top part.

As it is now, it only takes the getSTD function into account and does not do getdif. This is because
the getdif often just returned -1 since the highest peak was often the first peak, which
is what the getdif function was evaluating. 54/60 solutions were correctly classified as oscillatory. Only 
7 plots passed the cost function but were then found by the simple check to be damped oscillations,
making this method the best at filtering out ddamped oscillations.

I also ran it with the getdif funciton (it is not set that way now). That way 56/60 solutions were
correctly classified as oscillatory. However, a very large number (I did not document the exact number, oops)
of damped solutions had to be filtered out by my simple check. This is due to the issue mentioned above
where if the first peak in the fft was the tallest, a high value of getdif was returned. This could 
not be solved simply by increasing the cost function threshold since several oscillatory solutions
did not have their first peak as the highest and therefore had cost function values with magnitude
lower than 1.
-----------
=#

begin
    using Plots 
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
    using Unitful
    using BenchmarkTools, Profile, ProgressMeter
    using MultivariateStats, UMAP, TSne, StatsPlots
    using GlobalSensitivity, QuasiMonteCarlo
    using BifurcationKit, Setfield, LinearAlgebra, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes
    using Printf #And this
    default(lw = 2, size = (1000, 600))
end
"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
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

#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
        :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :y => 3631.050539219606]
p = [x[2] for x in psym]
    
#initial condition list
usym = [:L => 0, :Lp => 10^0.488788, :K => 10^-0.2895987, :P => 0.820348, :A => 10^0.42483, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]

#timespan for integration
const tspan = (0., 100.)
#solve the reduced ODEs
const prob = ODEProblem(fullrn, u0, tspan, p)
sol = solve(prob, saveat=0.1, save_idxs=1) #solve adaptively

"""
Gets the Fourier transform data and divides it by the height of the tallest peak. This is different 
from other programs where this is used, where it only divides by the length of the data set. As a result,
it returns the rfft data, but the height of the tallest peak is set to 1 for each data set.
- `y::Vector{Float64}` The set of data points in the solution to an ODE (sol.u)
"""
function getFrequencies(y::Vector{Float64})
    res = abs.(rfft(y))
    peaks = findmaxima(res,10)[2]
    if length(peaks) > 0
        maxPeak = maximum(peaks)
        return res ./ maxPeak #standardize amplitudes to have height of 1 for all 
    else
        return res
    end
end

"""
For each peak found in the Fourier transform, calculates the standard deviation of the height of the peak 
and the data points one index above and below. The standard deviations are then summed and divided by the number of peaks
in the Fourier transform. This is equivalent to the first term of the cost function in equation 4 of Pušnik et al. Very 
narrow and tall peaks in the Fourier transform are hence rewarded in the final cost function.
- `peakindxs` the indexes of the peaks in the Fourier transform of a solution
- `arrayData` the normalized absolute values of the rfft of a solution
"""
function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64})
    #Note: I got rid of the calculating of the window and just set the window size to three.
    sum_std = @inbounds sum(std(arrayData[max(1, ind - 1):min(length(arrayData), ind + 1)]) for ind in peakindxs)
    return sum_std / length(peakindxs)
end

#I got rid of the old get_dif function since it literally just evaluated the height of the first peak
#as long as there were two or more peaks.

"""Cost function based on Pušnik et al equation 4 with slight differences"""
function CostFunction(Y::ODESolution)
    fftData = getFrequencies(Y.u) #get the standardized rfft of the solution
    fftindexes, fftpeaks = findmaxima(fftData,10) #get the indexes and values of the peaks in the fft
    timeindexes = findmaxima(Y.u,10)[1] #get the times of the peaks in the fft
    if isempty(fftindexes) || length(timeindexes) < 2 #if there are no peaks in fft or less than 2 peaks in the untransformed solution return 0
        return 0.0
    end
    std = newgetSTD(fftindexes, fftData) #get the standard deviation of the peaks
    #diff = fftpeaks[1] #What Jonathan's cost function really was doing...
    return std #+ diff
end

"""Cost function that catches errors and returns 1.0 if there is an error"""
function eval_fitness_catcherrors(sol::ODESolution)
    Y = nothing
    try 
        if sol.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(sol) for x in array) || any(x==1 for array in isless.(sol, 0.0) for x in array)
            return 1.0
        end
    catch e 
        if isa(e, DomainError) #catch domain errors
            return 1.0
        else
            rethrow(e) #rethrow other errors
        end
    end
    return -CostFunction(sol)
end
"""
Determines if from time to end of the solution the standard deviation of the data points
is greater than 0.1% of the mean for the data points. This eliminates damped solutions
that make it through the cost function.
- `Y::ODESolution` the solution you want to check
- `time::Float64` the time point before the end that you want the algorithm to check from
"""
function simpleOscillationCheck(Y::ODESolution, time::Float64)
    timeIndex = findfirst(x -> x > time, Y.t) #Find index of time greater than 80
    solutionInterval = Y.u[timeIndex:end]
    meanVal = mean(solutionInterval)
    stdVal = std(solutionInterval)
    return stdVal > 0.001 * meanVal
end

#Currently set to Lp
numIterations = 10000
successCount = 0


for i in 1:numIterations
    u0[2] = rand(Random.seed!(i),Distributions.LogUniform(0.01, 100)) #Lp #1 for L, 2 for Lp
    u0[3] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.001, 100)) #K
    u0[4] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.01, 100)) #P
    u0[5] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.001, 100)) #A
    currentSol = solve(remake(prob, u0=u0), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
    name = "Score=$(round(eval_fitness_catcherrors(currentSol),digits=5)) Lp=$(u0[2]) K=$(u0[3]) P=$(u0[4])) A=$(u0[5]).png"
    solPlot = plot(currentSol, title="L vs t")
    fftSolPlot = Plots.plot((1:div(length(currentSol.t),2)),(abs.(rfft(currentSol.u))[2:end]), title="Fourier Transform")
    plot(solPlot,fftSolPlot, layout=(2,1))
    write("/Users/ezragreenberg/StandardizedCostFilter/SeedConcentrations.txt","$([i for i in u0[2:5]])\n","a")
    if eval_fitness_catcherrors(currentSol) < -0.01
        if simpleOscillationCheck(currentSol, 80.)
           savefig("/Users/ezragreenberg/StandardizedCostFilter/Pass/$name")
           successCount += 1
        else
            savefig("/Users/ezragreenberg/StandardizedCostFilter/Fail80/$name")
        end
        
    else
        savefig("/Users/ezragreenberg/StandardizedCostFilter/FailCost/$name")
    end
    if i % 500 == 0 
        println("$i iterations completed")
    end
end

println("Number of solutions classified as oscillatory: $successCount")