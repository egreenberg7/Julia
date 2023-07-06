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

rxn = @reaction_network rxn begin
    @parameters ka kb kcat k3 k4
    @species A(t) B(t) C(t) D(t) ABD(t) AB(t)
    (ka, kb), A + B <--> AB
    (kcat), AB --> C
    (k3, k4), AB + D <--> ABD
end

#This one works better
mmrxn1 = @reaction_network mmrxn begin
    @parameters Km kcat k3 k4
    @species A(t) B(t) C(t) D(t) ABD(t)
    (kcat / Km), A + B --> C
    (k3 / Km, k4), A + B + D <--> ABD
end

mmrxn2 = @reaction_network mmrxn2 begin
    @parameters Km kcat k3 k4
    @species A(t) B(t) C(t) D(t) ABD(t)
    (kcat / Km), A + B --> C
    (k3 * A * B / Km, k4), D <--> ABD
end

u0 = [:A => 4, :B => 2, :C => 0, :D => 2, :ABD => 0, :AB => 0]
mmu0 = u0[1:5]


params = [5, 2 ,10,2 ,2]
mmparams = [(params[2] + params[3]) / params[1], params[3], params[4],params[5]]

tspan2 = (0., 10.)
smallProb = ODEProblem(rxn, u0, tspan2, params)
mmProb1 = ODEProblem(mmrxn1, mmu0, tspan2, mmparams)
mmProb2 = ODEProblem(mmrxn2, mmu0, tspan2, mmparams)
smallsol = solve(smallProb, Rosenbrock23())
mmSmallsol1 = solve(mmProb1, Rosenbrock23())
mmSmallsol2 = solve(mmProb2, Rosenbrock23())

a,b,c = plot(smallsol, title = "Full"), plot(mmSmallsol1, title="mm1"), plot(mmSmallsol2,title="mm2")
plot(a,b,c)

diffVals = smallsol[1:2, :] - mmSmallsol[1:2, :]
absDiffs = abs.(diffVals)
sum1 = sum(absDiffs[1,:])

function latexRxn(rxnNet)
    rxnEqs = convert(ODESystem, rxnNet)
    txt = latexify(rxnEqs)
    render(txt)
end

latexRxn(mmrxn1)
latexRxn(mmrxn2)