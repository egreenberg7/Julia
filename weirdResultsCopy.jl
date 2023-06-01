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
#Currently set to an antidamped solution
usym = [:L => 0, :Lp => 10^1.486928914407112 , :K => 10^-0.20450631277896825, :P => 10.848484361685005, :A => 10^0.6597142624102299, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]


"""
Plot a the solution to a differential equation along with its Fourier transform, without the frequency of 0
- `sol` Solution being plotted
"""
function fftPlot(sol::ODESolution)
    solutionPlot = Plots.plot(sol)
    fourierPlot = Plots.plot(abs.(rfft(sol.u))[2:end])
    Plots.plot(solutionPlot,fourierPlot)
end

#timespan for integration
const tspan = (0., 100.)
#solve the reduced ODEs
const prob = ODEProblem(fullrn, u0, tspan, p)
sol = solve(prob, Rosenbrock23(), saveat=0.1, save_idxs=1) #solve adaptively

"""
Takes array of all initial concentrations (u0) and updates the initial concentrations that can be changed
- `curConcentrations` Array of all initial concetrations as currently set (u0)
- `newConcentrations` Array of new desired initial concentrations [L/Lp K P A]
- `isLp` Boolean, true if Lp in newConcentrations, false if L in newConcentrations
"""
function changeInitialConcentration(curConcentrations::Vector{Float64}, newConcentrations::Vector{Float64},isLp::Bool)
    if isLp
        curConcentrations[2] = newConcentrations[1] #Change Lp
    else
        curConcentrations[1] = newConcentrations[1] #Change L
    end
    curConcentrations[3] = newConcentrations[2]
    curConcentrations[4] = newConcentrations[3]
    curConcentrations[5] = newConcentrations[4]
end

function old_getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) #get difference between fft peak indexes
    #=
    :param indexes: the indexes of the peaks in the Fourier transform of a solution
    :param arrayData: the normalized absolute values of the rfft of a solution
    =#
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0 #? If there is only one peak, the score is set to 0. May not be necessary
    end
    println("The first peak is $(arrayData[indexes[1]])")
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1)) #? This is the negative of what is in the original paper
    sum_diff += arrayData[indexes[end]] #? The original paper also did not have a term just adding the one thing that was never added
    return sum_diff #/ length(indexes)
end

function getFrequencies(y::Vector{Float64})
    res = abs.(rfft(y))
    return res ./ cld(length(y), 2) #normalize amplitudes
    #? That length is defined to be half the length +1 of the input of the data set.
    #? Instead,
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

function CostFunction(Y::ODESolution)
    #get the fft of the solution
    fftData = getFrequencies(Y.u)
    fftindexes = findmaxima(fftData,10)[1] #get the indexes of the peaks in the fft
    timeindexes = findmaxima(Y.u,10)[1] #get the times of the peaks in the original solution
    if isempty(fftindexes) || length(timeindexes) < 2 #if there are no peaks, return 0
        return 0.0
    end
    std = getSTD(fftindexes, fftData) #get the standard deviation of the peaks
    println("Std: $std")
    diff = old_getDif(fftindexes, fftData) #get the difference between the peaks
    println("Diff $diff")
    return std + diff
end

#timespan for integration
const tspan = (0., 100.)
#solve the reduced ODEs
const prob = ODEProblem(fullrn, u0, tspan, p)
sol = solve(prob, Rosenbrock23(), saveat=0.1, save_idxs=1) #solve adaptively
fftPlot(sol)
CostFunction(sol)

dampedConcentrations = Vector{Float64}([10^0.5508305830906121, 10^0.4880537201023163, 0.6083519943817587, 10^1.4822439628752706])
changeInitialConcentration(u0,dampedConcentrations,true)
dampedSol = solve(remake(prob, u0=u0),Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
fftPlot(dampedSol)

falseNegative = [1.1503063547416548, 2.170245919643596, 0.5300251756760845, 4.861619177259957]
changeInitialConcentration(u0,falseNegative,true)
falseSol = solve(remake(prob, u0=u0),Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
fftPlot(falseSol)
falsefftData = abs.(rfft(falseSol.u))
peakIndices,peakVals = findmaxima(falsefftData,20)
std(falsefftData[40:42])
