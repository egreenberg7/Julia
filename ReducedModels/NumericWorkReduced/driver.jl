include("ReactionUtilities.jl")
using .SimpleReactions
using .FullReaction
using DifferentialEquations
include("CustomPeakFinder.jl")
include("Evaluator.jl")
using .Evaluator
using Plots
include("/Users/ezragreenberg/Documents/GitHub/Julia/ExperimentalFullModelWork/EvaluationFunctions.jl")

function makeProblem(reactionSystem, m)
    prob = ODEProblem(reactionSystem, m.randomU(), (0,600), m.randomP())
end

const Pprob = makeProblem(simplePRxn, SimpleReactions)
const Kprob = makeProblem(simpleKRxn, SimpleReactions)
const fullProb = makeProblem(fullRxn, FullReaction)


"""
    testReaction(prob, m, repetitions=100,tspan=600)
    Solves the reaction system over repetition times
    #Arguments
    - `prob` ODEProblem to solve
    - `m` module from which to draw parameter generators
"""
function testReaction(problem, m; repetitions=100000, tspan=600)
    for i in 1:repetitions
        p = m.randomP()
        p[end] = 10000  
        u = m.randomU()
        prob = remake(problem, u0=u, p=p, save_idxs=1)
        sol = solve(prob)
        if finalClassifier(sol, tspan)[1] < 0#fitness_function(p, u, Pprob, [1])[1] ≠ 0
            display(plot(sol, xlims=(0,tspan)))
            println("P: " * string(p))
            println("u0: " * string(u))
        end
        if i % (repetitions / 10) == 0
            println("$i repetitions completed")
        end 
    end
end
