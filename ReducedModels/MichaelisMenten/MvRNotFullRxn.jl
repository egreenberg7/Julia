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
    using Latexify
    #using BenchmarkTools, Profile, ProgressMeter
    #using MultivariateStats, UMAP, TSne, StatsPlots
    #using GlobalSensitivity, QuasiMonteCarlo
    using BifurcationKit, Setfield, LinearAlgebra, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes
    #using Printf #And this
    default(lw = 2, size = (1000, 600))
end

"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t)
    # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
    # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    (ka1,kb1), L + K <--> LK # L binding to kinase
    (kcat1), LK --> Lp + K # L phosphorylation by kinase into Lp
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
    kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    (kcat7), LpP --> L + P # L dephosphorylation by phosphatase
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    (kcat7), LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality
end

#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => 0.055, :kb1 => 19.8, :kcat1 => 241, :ka2 => 1, :kb2 => 0.95,
        :ka3 => 41, :kb3 => 193, :ka4 => 0.19, :kb4 => 0.13, :ka7 => 0.62, 
        :kb7 => 3.39, :kcat7 => 4.6, :y => 750]
p = [x[2] for x in psym]
    
#initial condition list
usym = [:L => 0.2, :Lp => 3.0, :K => 0.2, :P => 0.3, :A => 0.6, :LpA => 0.01, :LK => 0., 
        :LpP => 0., :LpAK => 0.1, :LpAP => 0.01, :LpAKL => 0.0, :LpAPLp => 0.0 ]
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


#={ dLK/dt = 0 and dLpP/dt = 0
mm_rn = @reaction_network mm_rn begin
    @parameters kcat1 Km1 ka2 kb2 ka3 kb3 ka4 kb4 kcat7 Km7 DF
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LpAK(t) LpAP(t) 

    kcat1 / Km1, L + K --> Lp + K # L phosphorylation by kinase into Lp using Michaelis-Menten
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    DF*kcat1 / Km1, LpAK + L --> Lp + LpAK # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A) using Michaelis-Menten
    kcat7 / Km7, Lp + P --> L + P # Lp dephosphorylation by phosphatase using Michaelis-Menten
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    DF * kcat7 / Km7, Lp + LpAP --> L + LpAP # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A) using Michaelis-Menten
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
    mmConc[1:6] = fullConc[1:6] #Sets initial conditions equal
    mmConc[7] = fullConc[9]
    mmConc[8] = fullConc[10]
    return mmConc
end
}=#



mmp = findMMParams(p)
mmu0 = findMMConc(u0)

#timespan for integration
const tspan = (0., 100.)
#solve the reduced ODEs
fullProb1 = ODEProblem(fullrn, u0, tspan, p)
mmProb = ODEProblem(mm_rn, mmu0, tspan, mmp)
fullsol = solve(fullProb1, saveat=0.1,save_idxs=[1,8])
mmsol = solve(mmProb, saveat=0.1,save_idxs=1)
a = plot(fullsol)
plot!(mmsol, label = "MM")

#I changed the ranges of concentrations to try to meet the
#substrate concentration assumption.
#TODO Fix 2D Michaelis Menten approximation since the assumptions are not met when the enzyme concentration is variable even if it is low.
numIterations = 100
for i in 1:numIterations
    u0[1] = rand(Random.seed!(4 * numIterations + i),Distributions.LogUniform(1, 100)) #L
    u0[2] = 0 #rand(Random.seed!(i),Distributions.LogUniform(0, 0)) #Lp
    u0[3] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.1, 10.0)) #K
    u0[4] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.1, 10.0)) #P
    u0[5] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.1, 10.0)) #A
    mmu0 = findMMConc(u0)
    currentFullSol = solve(remake(fullProb1, u0=u0), Rosenbrock23(), saveat=0.1, save_idxs=2, maxiters=10000, verbose=false)
    currentMMSol = solve(remake(mmProb, u0=mmu0), Rosenbrock23(), saveat=0.1, save_idxs=2, maxiters=10000, verbose=false)
    name = "L=$(u0[1]) Lp=$(u0[2]) K=$(u0[3]) P=$(u0[4]) A=$(u0[5]).png"
    SolPlot = plot(currentFullSol, label="Full Reaction", ylims=(0,1.25 * u0[1]))
    SolPlot = plot!(currentMMSol, label="MM Reaction", ylims=(0,1.25*u0[1]))
    savefig("~/mmTest/$name")
end


using Latexify
mm_eqs = convert(ODESystem, mm_rn)
txt = latexify(mm_eqs)
print(txt)
render(txt)

full_eqs = convert(ODESystem,fullrn)
txt2 = latexify(full_eqs)
render(txt2)

#=
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
=#