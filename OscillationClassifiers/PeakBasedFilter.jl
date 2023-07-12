#=
-----------
My code to generate 10000 graphs with Fourier transforms of oscillation plots with random initial concentrations
of Lp, K, P, and A varied on a log-uniform scale within experimentally realizable values. All other initial concentrations
were set to 0, and the rate constants were kept constant at values from one of your optimization runs.

I dispensed of the cost function and only looked at the real space solution. If there were greater than 3 peaks 
(which is a potential issue if we used different rate constants that could promote smaller frequencies) and the last
maxima varies by more than 0.1 from a peak in the middle of the timespan (which could be an issue if there is amplitude
mdoulation), and the last 20 seconds of solution are not steady state, then the solution is classified as oscillatory.

All solutions were correctly classified of the 10,000. Out of 100,000, there was 1 or 2 false positives (damped solutions), but false 
negatives were not checked.
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
Returns true if an ODE solution has more than 3 local maxima, and the height of the middle peak varies from
the height of the last peak by less than 0.1.
- `Y::ODESolution` ODESolution that you are testing
"""
function peakCheck(Y::ODESolution)
    peaks = findmaxima(Y.u, 10)[2]
    numPeaks = length(peaks)
    if numPeaks > 3 && abs((peaks[div(numPeaks,2)]-peaks[end])) < 0.1
        return true
    else
        return false
    end
end

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
    name = "Lp=$(u0[2]) K=$(u0[3]) P=$(u0[4])) A=$(u0[5]).png"
    solPlot = plot(currentSol, title="L vs t")
    fftSolPlot = Plots.plot((1:div(length(currentSol.t),2)),(abs.(rfft(currentSol.u))[2:end]), title="Fourier Transform")
    plot(solPlot,fftSolPlot, layout=(2,1))
    write("/Users/ezragreenberg/PeakFilter/SeedConcentrations.txt","$([i for i in u0[2:5]])\n","a")
    if peakCheck(currentSol)
        if simpleOscillationCheck(currentSol, 80.)
           savefig("/Users/ezragreenberg/PeakFilter/Pass/$name")
           successCount += 1
        else
            savefig("/Users/ezragreenberg/PeakFilter/Fail80/$name")
        end
        
    else
        savefig("/Users/ezragreenberg/PeakFilter/FailPeak/$name")
    end
    if i % 500 == 0 
        println("$i iterations completed")
    end
end

println("Number of solutions classified as oscillatory: $successCount")


#second check for only false positives over 100,000 iterations
for i in 1:100000
    u0[2] = rand(Random.seed!(i),Distributions.LogUniform(0.01, 100)) #Lp #1 for L, 2 for Lp
    u0[3] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.001, 100)) #K
    u0[4] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.01, 100)) #P
    u0[5] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.001, 100)) #A
    currentSol = solve(remake(prob, u0=u0), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)   
    if simpleOscillationCheck(currentSol, 80.) && peakCheck(currentSol)
        name = "Lp=$(u0[2]) K=$(u0[3]) P=$(u0[4])) A=$(u0[5]).png"
        solPlot = plot(currentSol, title="L vs t")
        fftSolPlot = Plots.plot((1:div(length(currentSol.t),2)),(abs.(rfft(currentSol.u))[2:end]), title="Fourier Transform")
        plot(solPlot,fftSolPlot, layout=(2,1))
        savefig("/Users/ezragreenberg/PeakFilter/LargeTrialPass/$name")
        successCount += 1
    end
    if i % 1000 == 0 
        println("$i iterations completed")
    end
end