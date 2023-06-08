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
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    # using BifurcationKit, Setfield, LinearAlgebra, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes
    using Printf #And this
    default(lw = 2, size = (1000, 600))
    include("../UTILITIES/ODEProbMaker.jl")
    include("../UTILITIES/GA_functions.jl")
    include("../UTILITIES/EvaluationFunctions.jl")
    include("../UTILITIES/EvolutionaryOverloads.jl")
    include("../UTILITIES/ReactionNetwork.jl")
end

# """Full oscillator model"""
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

#     #previously excluded reactions, all possible combinations possible in vitro
#     (ka2,kb2), Lp + AK <--> LpAK
#     (ka2*y,kb2), Lp + AKL <--> LpAKL
#     (ka2,kb2), Lp + AP <--> LpAP
#     (ka2*y,kb2), Lp + APLp <--> LpAPLp
#     (ka3,kb3), A + K <--> AK
#     (ka4,kb4), A + P <--> AP
#     (ka3,kb3), A + LK <--> AKL
#     (ka4,kb4), A + LpP <--> APLp
#     (ka3*y,kb3), LpA + LK <--> LpAKL
#     (ka4*y,kb4), LpA + LpP <--> LpAPLp
#     (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
#     kcat1, AKL --> Lp + AK #phosphorylation of lipid
#     (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
#     kcat7, APLp --> L + AP #dephosphorylation of lipid
# end
fullrn = make_fullrn()

#parameter list with experimental values
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
# psym = [:ka1 => 0., :kb1 => 0., :kcat1 => 15. * ka1 - kb1, :ka2 => 0.7*10^-3, :kb2 => 2.0*10^-3,
#         :ka3 => 0.118, :kb3 => 0.0609, :ka4 => 0., :kb4 => 29. * ka4, :ka7 => 0, 
#         :kb7 => 39. * ka7 - kcat7, :kcat7 => 85.3, :y => 0.]
# p = [x[2] for x in psym]

psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
        :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :y => 3631.050539219606]
p = [x[2] for x in psym]

function constructParams(ka1, kb1, ka4, ka7, y)
end


#initial condition list
usym = [:L => 5, :K => 0.1, :P => 0.1, :A => 0.1, :Lp => 5, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.01, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]

exp_prob = make_nerdssODEProb(fullrn; p=p, u0=u0)

#Play with make_fitness_function if you want to make your own and pass it into eval_fitness_function_factor_thing

IC_Constraints = define_initialcondition_constraints()
myGA_problem = GAProblem(IC_Constraints,exp_prob)
record = run_GA(myGA_problem)

recordDF=DataFrame(record)
argmin(recordDF[!,:fit])
mostFit = recordDF[20000,:]
sol = solve(remake(exp_prob, u0=vcat(mostFit[:ind],[5,0,0,0,0,0.01,0,0,0,0,0,0])), Rosenbrock23(),saveat=0.1)
plot(sol)
fftsol = abs.(rfft(sol[1,:]))
plot(fftsol[2:end])