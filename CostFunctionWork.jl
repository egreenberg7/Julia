#My work on testing the cost function.

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

#Block of code setting up reactions.

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
    Plots.plot((2:length(sol.t)),(abs.(fft(sol.u[2:end]))))


## Genetic algorithm to find the best parameters for the reduced model ## 
function old_getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) #get difference between fft peak indexes
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0
    end
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1))
    sum_diff += arrayData[indexes[end]]
    return sum_diff/length(indexes)
end

function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
    #=

    :param peakindxs: the indexes of the peaks in the Fourier transform of a solution
    :param arrayData: the normalized absolute values of the rfft of a solution
    =#
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

"""Old cost function to be plugged into eval_fitness wrapper"""
function CostFunction(Y::ODESolution)
    #get the fft of the solution
    fftData = getFrequencies(Y.u)
    fftindexes = findmaxima(fftData,10)[1] #get the indexes of the peaks in the fft
    timeindexes = findmaxima(Y.u,10)[1] #get the times of the peaks in the original solution
    if isempty(fftindexes) || length(timeindexes) < 2 #if there are no peaks, return 0
        return 0.0
    end
    std = getSTD(fftindexes, fftData, 0.001) #get the standard deviation of the peaks
    println("Std: $std")
    diff = old_getDif(fftindexes, fftData) #get the difference between the peaks
    println("Diff $diff")
    return std + diff
end
|