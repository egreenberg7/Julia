using DifferentialEquations
using DataFrames
using CSV
using Catalyst
using Evolutionary
using FFTW
using ProgressMeter
using LinearAlgebra
using Distributions
using Peaks
using Plots
include("../UTILITIES/ODEProbMaker.jl")
include("../UTILITIES/ReactionNetwork.jl")
include("../UTILITIES/EvaluationFunctions.jl")
include("../UTILITIES/EvolutionaryOverloads.jl")
include("../UTILITIES/GA_functions.jl")

#Experimental values, see ipynb file
p=[0.1, 0.31622776601683794, 1.183772233983162, 0.0007, 0.002, 0.118, 0.0609, 0.0031622776601683794, 0.091706052144883, 6.309573444801933, 160.7733643472754, 85.3, 10000.0]
u0=[2.51189,0.251189,0.01,0.251189,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.0]
fullrn = make_fullrn()
exp_prob = make_nerdssODEProb(fullrn; p=p,u0=u0, tspan=600.0)

#Play with make_fitness_function if you want to make your own and pass it into eval_fitness_function_factor_thing
IC_Constraints = define_initialcondition_constraints()
myGA_problem = GAProblem(IC_Constraints,exp_prob)
record = run_GA(myGA_problem)

for i in eachrow(record)
    curFit = i 
    sol = solve(remake(exp_prob, u0=vcat(curFit[:ind],[0,0,0,0,0,0.0,0,0,0,0,0,0])), Rodas4(),saveat=0.1, save_idxs = 1)
    solutionInterval = sol.u[4800:end]
    meanVal = mean(solutionInterval)
    stdVal = std(solutionInterval)
    p1 = plot(sol[1,:])
    display(p1)
    savefig(p1, "graphStorage/Graph$(round(stdVal / meanVal, sigdigits = 3, base = 10)),$(i.fit).png")
end