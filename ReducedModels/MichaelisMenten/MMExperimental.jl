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
    #using BenchmarkTools, Profile, ProgressMeter
    #using MultivariateStats, UMAP, TSne, StatsPlots
    #using GlobalSensitivity, QuasiMonteCarlo
    using BifurcationKit, Setfield, LinearAlgebra, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes
    #using Printf #And this
    default(lw = 2, size = (1000, 600))
end

# """Full oscillator model, exluded reactions"""
# fullrn = @reaction_network fullrn begin
#     @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
#     @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
#     # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
#     # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
#     (ka1,kb1), L + K <--> LK # L binding to kinase
#     kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
#     (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
#     (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
#     (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
#     kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
#     (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
#     kcat7, LpP --> L + P # L dephosphorylation by phosphatase
#     (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
#     (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
#     kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality
# end

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
kcat1est = 15. * ka1est - kb1est
ka4est = 0.12 #Estimating LpA bindking to K (ka3) similar to LpA binding to L
ka7est = 0.9 #Pretty much guess, around ka1
yest = 750

#parameter list Changed this around
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
        :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
        :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => yest]
p = [x[2] for x in psym]

    
#initial condition list
usym = [:L => 5, :Lp => 0.1, :K => 0.1, :P => 0.1, :A => 0.1, :LpA => 0.0, :LK => 0, 
        :LpP => 0, :LpAK => 0., :LpAP => 0., :LpAKL => 0., :LpAPLp => 0, :AK => 0, :AP => 0, 
        :AKL => 0, :APLp => 0]
u0 = [x[2] for x in usym]


"""Only setting dLK/dt to 0 and dLpAKL/dt = 0"""
mm_rn = @reaction_network mm_rn begin
    @parameters kcat1 Km1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LpP(t) LpAK(t) LpAP(t) LpAPLp(t)

    kcat1 / Km1, L + K --> Lp + K # L phosphorylation by kinase into Lp using Michaelis-Menten
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    y*kcat1 / Km1, LpAK + L --> Lp + LpAK # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A) using Michaelis-Menten
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    (kcat7), LpP --> L + P # L dephosphorylation by phosphatase
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    (kcat7), LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionalityend
end
"""
Finds the values of the Michaelis-Menten parameters given
a list of parameteres in the full model.
"""
function findMMParams(fullParams::Vector{Float64})
    mmParams = zeros(12)
    mmParams[1] = fullParams[3] #kcat1
    mmParams[2] = sum(fullParams[2:3]) / fullParams[1] #Km1 = (kb1 + kcat) / ka1
    mmParams[3:12] = fullParams[4:13] #ka2,kb2,ka3,kb3,ka4,kb4, ka7, kb7, kcat7, y
    return mmParams
end

"""
Finds the initial concentrations of the Michaelis-Menten species given
the initial concentrations in the full model.
"""
function findMMConc(fullConc::Vector{Float64})
    mmConc = zeros(10)
    mmConc[1:6] = fullConc[1:6] #Sets initial conditions equal
    mmConc[7:9] = fullConc[8:10]
    mmConc[10] = fullConc[12] 
    return mmConc
end

mmp = findMMParams(p)
mmu0 = findMMConc(u0)

#timespan for integration
const tspan = (0., 100.)
#solve the reduced ODEs
fullProb1 = ODEProblem(fullrn, u0, tspan, p)
mmProb = ODEProblem(mm_rn, mmu0, tspan, mmp)
fullsol = solve(fullProb1, saveat=0.1, save_idxs=1)
mmsol = solve(mmProb, saveat=0.1, save_idxs=1)
a = plot(fullsol, label = "full")
plot!(mmsol, label = "MM")

#I changed the ranges of concentrations to try to meet the
#substrate concentration assumption.
#TODO Fix 2D Michaelis Menten approximation since the assumptions are not met when the enzyme concentration is variable even if it is low.
numIterations = 100
for i in 1:numIterations
    u0[1] = rand(Random.seed!(4 * numIterations + i),Distributions.LogUniform(1, 100)) #L
    u0[2] = 0.0 #rand(Random.seed!(i),Distributions.LogUniform(0, 0)) #Lp
    u0[3] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.001, 1)) #K
    u0[4] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.01, 1)) #P
    u0[5] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.001, 1)) #A
    mmu0 = findMMConc(u0)
    currentFullSol = solve(remake(fullProb1, u0=u0), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
    currentMMSol = solve(remake(mmProb, u0=mmu0), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
    name = "L=$(u0[1]) Lp=$(u0[2]) K=$(u0[3]) P=$(u0[4]) A=$(u0[5]).png"
    SolPlot = plot(currentFullSol, label="Full Reaction", ylims=(0,1.25 * u0[1]))
    SolPlot = plot!(currentMMSol, label="MM Reaction", ylims=(0,1.25*u0[1]))
    savefig("~/mmTest/$name")
end