begin 
    using Plots; #theme(:juno)
    # using Compose
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using Statistics
    # using Peaks
    using FindPeaks1D
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    # using Unitful
    # using Unitful: µM, M, nm, µm, s, μs, Na, L, 𝐍
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra

    using Setfield
    
    using ColorSchemes, Plots.PlotMeasures
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # plotlyjs()
    # import CairoMakie as cm 
    # gr()
    # push!(LOAD_PATH, "../../UTILITIES")

    #* import the overloads for Evolutionary.jl
    include("../../UTILITIES/EvolutionaryOverloads.jl")

    #* import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    #* import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    #* import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    #* import the plotting functions
    include("../../UTILITIES/TestBenchPlotUtils.jl")

    # include("../../UTILITIES/UnitTools.jl")

    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

function make_ODE_problem()
    tspan = (0., 2000.0)

    fullrn = make_fullrn()

    ogprob = ODEProblem(fullrn, [], tspan, [])

    de = modelingtoolkitize(ogprob)

    ogprobjac::ODEProblem = ODEProblem(de, [], tspan, jac=true)
    return ogprobjac, ogprob
end

ogprobjac, ogprob = make_ODE_problem();


ogjacsol = solve(ogprobjac, Rosenbrock23(), saveat=0.1, save_idxs= [6, 9, 10, 11, 12, 15, 16])

@code_warntype eval_all_fitness(rand(17), ogprobjac)

@code_warntype CostFunction(ogjacsol)

@btime CostFunction($ogjacsol)


solu = map(sum, ogjacsol.u)



@btime CostFunction($solu, $ogjacsol.t)

# tstart = cld(length(ogjacsol[1,:]), 1000)

# std(ogjacsol[1,end-tstart:end])
# # plot(ogjacsol)

# maxpeaks = findextrema(ogjacsol[1,:]; height=1e-2, distance=2)
# # minpeaks = findextrema(ogjacsol[1,:]; height=-1e-2, distance=2, find_maxima=false)
# # minpeaks = findextrema(flip_about_mean(ogjacsol[1,:]); height=-1e-2, distance=2, find_maxima=false)


# @btime CostFunction($ogjacsol)
# @code_warntype CostFunction(ogjacsol)


ogprobjac, ogprob = make_ODE_problem();



param_constraints = ParameterConstraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
ic_constraints = InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

allconstraints = AllConstraints(param_constraints, ic_constraints)



gaproblem = GAProblem(allconstraints, ogprobjac)

@code_warntype GAProblem(allconstraints, ogprobjac)






@btime generate_population($allconstraints, 5000)
pop = generate_population(allconstraints, 5000)
@btime generate_population!($pop, $allconstraints)

@code_llvm generate_population(allconstraints, 5000)



# Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int})
    function fitness_function(input::Vector{Float64})
        # Create a new input vector that includes the fixed inputs.
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_triplet))

        # Keep track of the number of fixed inputs that have been inserted.
        fixed_inputs_inserted = 0

        for i in eachindex(new_input)
            if i in triplet_idxs
                # If the current index matches the index of a fixed input, insert the fixed input.
                new_input[i] = fixed_input_triplet[fixed_inputs_inserted + 1]
                fixed_inputs_inserted += 1
            else
                # Otherwise, insert the next value from the input vector.
                new_input[i] = input[i - fixed_inputs_inserted]
            end
        end

        return evalfunc(new_input, prob)
    end
    return fitness_function
end


# Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs_bothparamsIC(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int}, fixedDF=1000.)
    function fitness_function(input::Vector{Float64})
        # Create a new input vector that includes the fixed inputs.
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_triplet))

        # Keep track of the number of fixed inputs that have been inserted.
        fixed_inputs_inserted = 0

        for i in eachindex(new_input)
            if i in triplet_idxs
                # If the current index matches the index of a fixed input, insert the fixed input.
                new_input[i] = fixed_input_triplet[fixed_inputs_inserted + 1]
                fixed_inputs_inserted += 1
            else
                # Otherwise, insert the next value from the input vector.
                new_input[i] = input[i - fixed_inputs_inserted]
            end
        end


        insert!(new_input, 13, fixedDF)

        return evalfunc(new_input, prob)
    end
    return fitness_function
end


fixed_inputs = (L = 100.0, K = 1.0, P = 1.0, A = 10.0)

set_fixed_constraints!(gaproblem.constraints, fixed_inputs)



"""
    Takes a named tuple of fixed inputs and a GAProblem and sets the constraints of the GAProblem to the fixed inputs.
    Returns `DataFrame` of the optimization results.
"""
function test_fixedparam(gaprob::GAProblem, fixedDF=1000.; fixed_inputs)

    constraints = gaprob.constraints

    set_fixed_constraints!(gaprob; fixed_inputs..., DF=fixedDF)


    Random.seed!(1234)

    initial_population = generate_population(constraints, 10000)

    ga_results = run_GA(gaprob, initial_population; iterations = 5)

    oscillatory_points_df = make_ga_dataframe(ga_results, constraints) 
    num_oscillatory_points = nrow(oscillatory_points_df)
    @info num_oscillatory_points

    return oscillatory_points_df
end


fixed_inputs = (L = 100.0, K = 1.0, P = 1.0, A = 10.0)

testfixed_df = test_fixedparam(gaproblem; fixed_inputs)



set_fixed_constraints!(allconstraints, fixed_inputs)

fitfunc = gaproblem.fitness_function

input = ogprob.p

fitfunc(input)
#########################

#* OBSERVATIONS FROM THE TESTBENCH
#* No false positives or negatives currently with no fixed param GA optimization
#* Bias towards high frequency I think
#* Last half STD check seems necessary, not sure about trim
#* Not much difference between Rosenbrock23 and Tsit5, or dt = 0.1 and dt = 0.01
#* Peak height threshold seems to be important





#!NOTES 
#* Need to play around with PM mutation scheme. Not working as well as BGA, but has benefits
#* Adapt some sort of PM scheme for the initial population generation. Need to make sure it's diverse enough
#* Need to make this GA more rigorous as far as what it tells about a regime. Having more oscillatory points but where each point isn't that different from the next isn't a good comparitive metric.
#* DE instead of GA? Want global optimization and search, not trapped in local minima
function testbench(constraints::ConstraintType, prob::ODEProblem)
    test_gaproblem = GAProblem(constraints, prob)
    Random.seed!(1234)
    test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 5, show_trace=true)

    # avg_fitness = mean(test_results.fitvals)
    # @info "Average fitness: $avg_fitness"
    # avg_period = mean(test_results.periods)
    # @info "Average period: $avg_period"
    # avg_amplitude = mean(test_results.amplitudes)
    # @info "Average amplitude: $avg_amplitude"
    if length(test_results.fitvals) == 0
        @info "No oscillatory points found."
    else 
        test_results_df = make_ga_dataframe(test_results, prob)
        return test_results_df#, avg_fitness, avg_period, avg_amplitude
    end
end

#* now testing Rosenbrock23 and new peak finder in getPerAmp for the ringing solutions
@code_warntype testbench(param_constraints, ogprobjac)

ogprobjac = remake(ogprobjac, u0 = [[100., 0.2, 0.2, 4.64]; zeros(12)])

test_results_df = testbench(param_constraints, ogprobjac)

test_results_df = testbench(allconstraints, ogprobjac)


function testBGA(valrange::Vector, m::Int = 2)
    prob = 1.0 / m
    function mutation(recombinant::T;
                      rng::AbstractRNG=Random.default_rng()
                     ) where {T <: AbstractVector}
        d = length(recombinant)
        @assert length(valrange) == d "Range matrix must have $(d) columns"
        δ = zeros(m)
        for i in 1:length(recombinant)
            for j in 1:m
                δ[j] = (rand(rng) < prob) ? δ[j] = 2.0^(-j) : 0.0
            end
            if rand(rng, Bool)
                recombinant[i] += sum(δ)*valrange[i]
            else
                recombinant[i] -= sum(δ)*valrange[i]
            end
        end
        return recombinant
    end
    return mutation
end

valrange = fill(2.0, 13)

mutationfunc = testBGA(valrange, 1)

params = copy(ogprob.p)

mutationfunc(params)


plotboth(test_results_df[9,:], ogprob)

@btime testbench($param_constraints, $ogprobjac)

plot_everything(test_results_df, ogprob; setnum=16, label="Testing", jump = 10)


stats_df = describe(test_results_df)
show(stats_df, allrows=true)

using StatsPlots

@df stats_df plot(:min, :max)

test_results_df.amp_percentage = test_results_df.amp./test_results_df.A


# split_dataframe!(test_results_df, ogprobjac)

CSV.write("OscillatorPaper/FigureGenerationScripts/high_amp_Amem.csv", test_results_df)

@btime testbench($param_constraints, $ogprob)
@btime testbench($param_constraints, $ogprobjac)

plot_everything(test_results_df, ogprob; setnum=15, label="TestingSTDWindow", jump = 10)


plotboth(test_results_df[1,:], ogprob)



#* measure A in solution vs A membrane 
#* quadruplet fixed search initial conditions 


# get data 
testdf = CSV.read("OscillatorPaper/FigureGenerationScripts/test.csv", DataFrame)

#combine all parameter columns into one column of vectors



p = [param for param in test_results_df[1763, Between(:ka1, :DF)]]
u0 = [ic for ic in test_results_df[1763, Between(:L, :A)]]

reprob = remake(ogprob, p = p, u0 = [u0; zeros(length(ogprob.u0) - length(u0))])

sol = solve(reprob, Rosenbrock23(), saveat=0.1, save_idxs = [6, 9, 10, 11, 12, 15, 16])

Amem = map(sum, sol.u)

findextrema(Amem; height=1e-2, distance=2)
findextrema(Amem; height=0.0, distance=2, find_maxima=false)

CostFunction(sol)
plot(sol)




testdf = CSV.read("/Users/jonathanfischer/Desktop/PhD_ThesisWork/Julia/OscillatorPaper/FigureGenerationScripts/DF=100.0.csv", DataFrame)

plot_everything(testdf, ogprob; setnum=15, label="TestingSTDWindow", jump = 10)

#! diversity metric
#* PCA?
#* minimum distance between points in the population
#* need to figure out convergence metric to compare 