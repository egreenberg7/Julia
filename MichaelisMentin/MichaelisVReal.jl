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

#parameter list Changed this around
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => 2.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
        :kb7 => 50.04411989797898752, :kcat7 => 42.288085868394326, :y => 3631.050539219606]
p = [x[2] for x in psym]
    
#initial condition list
usym = [:L => 5, :Lp => 5, :K => 0.1, :P => 0.1, :A => 0.1, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]

mm_rn = @reaction_network mm_rn begin
    @parameters kcat1 Km1 ka2 kb2 ka3 kb3 ka4 kb4 kcat7 Km7 DF
    @species L(t) K(t) P(t) A(t) Lp(t) LpA(t) LpAK(t) LpAP(t) 

    kcat1*(K + LpAK)/(K*(Km1+L)), L + K --> Lp + K # L phosphorylation by kinase into Lp using Michaelis-Menten
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (kcat1*(K + LpAK) * DF)/(LpAK*(Km1+L)), LpAK + L --> Lp + LpAK # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A) using Michaelis-Menten
    (kcat7*(P+LpAP))/((Km7+Lp) * P), Lp + P --> L + P # Lp dephosphorylation by phosphatase using Michaelis-Menten
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    DF * (kcat7*(P+LpAP))/((Km7+Lp) * LpAP), Lp + LpAP --> L + LpAP # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A) using Michaelis-Menten
end

"""
Finds the values of the Michaelis-Menten parameters given
a list of parameteres in the full model.
"""
function findMMParams(fullParams::Vector{Float64})
    mmParams = zeros(11)
    mmParams[1] = fullParams[3] #kcat1
    mmParams[2] = sum(fullParams[2:3]) / fullParams[1] #Km1 = (kb1 + kcat) / ka1
    mmParams[3:8] = fullParams[4:9] #ka2,kb2,ka3,kb3,ka4,kb4
    mmParams[9] = fullParams[12] #kcat7
    mmParams[10] = sum(fullParams[11:12]) / fullParams[10] #Km7 = (kb7+kcat7) / ka7
    mmParams[11] = fullParams[13] #df
    return mmParams
end

"""
Finds the initial concentrations of the Michaelis-Menten species given
the initial concentrations in the full model.
"""
function findMMConc(fullConc::Vector{Float64})
    mmConc = zeros(8)
    mmConc[1:5] = fullConc[1:5] #Sets initial conditions equal
    return mmConc
end

mmp = findMMParams(p)
mmu0 = findMMConc(u0)

#timespan for integration
const tspan = (0., 100.)
#solve the reduced ODEs
const fullProb = ODEProblem(fullrn, u0, tspan, p)
const mmProb = ODEProblem(mm_rn, mmu0, tspan, mmp)
fullsol = solve(fullProb, saveat=0.1, save_idxs=1)
mmsol = solve(mmProb, saveat=0.1, save_idxs=1)
a,b = plot(fullsol), plot(mmsol)
#plot(a,b)

#I changed the ranges of concentrations to try to meet the
#substrate concentration assumption.
#TODO Fix 2D Michaelis Menten approximation since the assumptions are not met when the enzyme concentration is variable even if it is low.
numIterations = 100
for i in 1:numIterations
    u0[1] = rand(Random.seed!(4 * numIterations + i),Distributions.LogUniform(1, 100)) #L
    u0[2] = rand(Random.seed!(i),Distributions.LogUniform(1, 100)) #Lp
    u0[3] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.001, 0.1)) #K
    u0[4] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.01, 0.1)) #P
    u0[5] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.001, 0.1)) #A
    mmu0 = findMMConc(u0)
    currentFullSol = solve(remake(fullProb, u0=u0), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
    currentMMSol = solve(remake(mmProb, u0=mmu0), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
    name = "L=$(u0[1]) Lp=$(u0[2]) K=$(u0[3]) P=$(u0[4]) A=$(u0[5]).png"
    fullSolPlot = plot(currentFullSol, title="Full Reaction")
    mmSolPlot = plot(currentMMSol, title="MM Reaction")
    plot(fullSolPlot,mmSolPlot, layout=(2,1))
    savefig("~/mmTest/$name")
end

mm_eqs = convert(ODESystem, mm_rn)
txt = latexify(mm_eqs)
print(txt)
render(txt)