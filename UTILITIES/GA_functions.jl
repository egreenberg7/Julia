#< CONSTRAINT RANGE STRUCTS

abstract type ConstraintType end

"""
    ConstraintRange

Struct for defining parameter ranges. Each instance contains a name, and a range defined by a minimum and maximum value.

# Fields
- `name::String`: Name of the parameter.
- `symbol::Symbol`: Symbol of the parameter.
- `min::Float64`: Minimum value for the parameter.
- `max::Float64`: Maximum value for the parameter.
- `nominal::Float64`: Nominal value for the parameter.
"""
struct ConstraintRange
    name::String
    symbol::Symbol
    min::Float64
    max::Float64
    nominal::Float64
end

"""
    ParameterConstraints

Struct encapsulating parameter constraints. Each field represents a different parameter, holding a `ConstraintRange` object that defines the valid range for that parameter.
"""
mutable struct ParameterConstraints <: ConstraintType
    ranges::Vector{ConstraintRange}
end

"""
    InitialConditionConstraints

Struct encapsulating initial condition constraints. Each field represents a different initial condition, holding a `ConstraintRange` object that defines the valid range for that initial condition.
"""
mutable struct InitialConditionConstraints <: ConstraintType 
    ranges::Vector{ConstraintRange} 
end
#> END 



#< CONSTRAINT RANGE CONSTRUCTORS
"""
    define_parameter_constraints(; kwargs...)

Define parameter constraints. Each keyword argument represents a different parameter, where the value is a tuple defining the valid range for that parameter.

# Example
```julia
constraints = define_parameter_constraints(
    karange = (-3.0, 1.0), 
    kbrange = (-3.0, 3.0), 
    kcatrange = (-3.0, 3.0), 
    dfrange = (1.0, 5.0)
)
```
"""
function define_parameter_constraints(; karange = (1e-3, 1e1), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e3, 1e5),
    nominalvals = (;ka1 = 0.009433439939827041, kb1 = 2.3550169939427845, kcat1 = 832.7213093872278, ka2 = 12.993995997539924, kb2 = 6.150972501791291,
            ka3 = 1.3481451097940793, kb3 = 0.006201726090609513, ka4 = 0.006277294665474662, kb4 = 0.9250191811994848, ka7 = 57.36471615394549, 
            kb7 = 0.04411989797898752, kcat7 = 42.288085868394326, DF = 3631.050539219606)
            )
    #* Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale


    return ParameterConstraints(
        [
        ConstraintRange("ka1", :ka1, ka_min, ka_max, nominalvals[1]),
        ConstraintRange("kb1", :kb1, kb_min, kb_max, nominalvals[2]),
        ConstraintRange("kcat1", :kcat1, kcat_min, kcat_max, nominalvals[3]),
        ConstraintRange("ka2", :ka2, ka_min, ka_max, nominalvals[4]),
        ConstraintRange("kb2", :kb2, kb_min, kb_max, nominalvals[5]),
        ConstraintRange("ka3", :ka3, ka_min, ka_max, nominalvals[6]),
        ConstraintRange("kb3", :kb3, kb_min, kb_max, nominalvals[7]),
        ConstraintRange("ka4", :ka4, ka_min, ka_max, nominalvals[8]),
        ConstraintRange("kb4", :kb4, kb_min, kb_max, nominalvals[9]),
        ConstraintRange("ka7", :ka7, ka_min, ka_max, nominalvals[10]),
        ConstraintRange("kb7", :kb7, kb_min, kb_max, nominalvals[11]),
        ConstraintRange("kcat7", :kcat7, kcat_min, kcat_max, nominalvals[12]),
        ConstraintRange("DF", :DF, df_min, df_max, nominalvals[13])
        ]
    )
end

define_parameter_constraints(prob::ODEProblem; karange = (1e-3, 1e1), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e3, 1e5)) = 
                        define_parameter_constraints(; karange=karange, kbrange=kbrange, kcatrange=kcatrange, dfrange=dfrange, nominalvals = prob.p)


"""
    define_initialcondition_constraints(; kwargs...)

Define initial condition constraints. Each keyword argument represents a different initial condition, where the value is a tuple defining the valid range for that initial condition.

# Example
```julia
constraints = define_initialcondition_constraints(
    lipidrange = (0.1, 10.0), 
    kinaserange = (0.1, 10.0), 
    phosphataserange = (0.1, 10.0), 
    ap2range = (0.1, 10.0)
)
```
"""
function define_initialcondition_constraints(;lipidrange = (0.1, 10.0), kinaserange = (0.1, 5.0), phosphataserange = (0.1, 5.0), ap2range = (0.1, 10.0),
                                                nominalvals = (;L = 3.0, K = 0.5, P = 0.3, A = 2.0)
                                            )
    # Define parameter constraint ranges
    lipid_min, lipid_max = lipidrange  # uM
    kinase_min, kinase_max = kinaserange  # uM
    phosphatase_min, phosphatase_max = phosphataserange # uM
    ap2_min, ap2_max = ap2range # uM

    return InitialConditionConstraints(
        [
        ConstraintRange("PIP+PIP2", :L, lipid_min, lipid_max, nominalvals[1]),
        ConstraintRange("Kinase", :K, kinase_min, kinase_max, nominalvals[2]),
        ConstraintRange("Phosphatase", :P, phosphatase_min, phosphatase_max, nominalvals[3]),
        ConstraintRange("AP2", :A, ap2_min, ap2_max, nominalvals[4])
        ]
    )
end

define_initialcondition_constraints(prob::ODEProblem; lipidrange = (0.1, 10.0), kinaserange = (0.1, 5.0), phosphataserange = (0.1, 5.0), ap2range = (0.1, 10.0)) = 
                                    define_initialcondition_constraints(;lipidrange=lipidrange, kinaserange=kinaserange, phosphataserange= phosphataserange, ap2range=ap2range, nominalvals = prob.u0)
#> END




#< GA PROBLEM TYPE
"""
    GAProblem{T <: ConstraintType}

Struct encapsulating a Genetic Algorithm (GA) optimization problem. It holds the constraints for the problem, the ODE problem to be solved, the fitness function, and any additional keyword arguments.

# Fields
- `constraints::T`: Constraints for the problem. Either `ParameterConstraints` or `InitialConditionConstraints`.
- `ode_problem::ODEProblem`: ODE problem to be solved.
- `fitness_function::Function`: Fitness function, automatically generated with constructor
- `options::NamedTuple`: Additional keyword arguments for `Evolutionary.Options`.
"""
struct GAProblem{T <: ConstraintType}
    constraints::T
    ode_problem::ODEProblem
    eval_function::Function

    function GAProblem(constraints::ParameterConstraints, ode_problem::ODEProblem) 
        new{ParameterConstraints}(constraints, ode_problem, eval_param_fitness)
    end

    function GAProblem(constraints::InitialConditionConstraints, ode_problem::ODEProblem) 
        new{InitialConditionConstraints}(constraints, ode_problem, eval_ic_fitness)
    end
end


function Base.show(io::IO, ::MIME"text/plain", prob::GAProblem) #TODO add labels for nominal values
    printstyled(io, "GAProblem with constraints:\n"; bold = true, underline=true, color = :green)
    printstyled(io, prob.constraints, "\n")
    printstyled(io, "\nNominal parameter values:\n"; bold = true, color = :blue)
    printstyled(io, prob.ode_problem.p, "\n")
    printstyled(io, "\nNominal initial conditions:\n"; bold = true, color = :blue)
    printstyled(io, prob.ode_problem.u0, "\n")
end

#> END

#< POPULATION GENERATION METHODS
"""
    generate_population(constraints::ParameterConstraints, n::Int)

Generate a population of `n` individuals for the given parameter constraints. Each individual is sampled from a log-uniform distribution within the valid range for each parameter.

# Example
```julia
constraints = define_parameter_constraints(karange = (-3.0, 1.0), kbrange = (-3.0, 3.0))
population = generate_population(constraints, 100)
```
"""
function generate_population(constraint::ParameterConstraints, n::Int)
    population = [exp10.(rand(Uniform(log10(conrange.min), log10(conrange.max)), n)) for conrange in constraint.ranges]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end

"""
generate_population(constraints::InitialConditionConstraints, n::Int)

Generate a population of `n` individuals for the given initial condition constraints. Each individual is sampled from a uniform distribution within the valid range for each initial condition.

# Example
```julia
constraints = define_initialcondition_constraints(lipidrange = (0.1, 10.0), kinaserange = (0.1, 10.0))
population = generate_population(constraints, 100)
```
"""
function generate_population(constraint::InitialConditionConstraints, n::Int)
    population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint.ranges]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end

"""For calculating volume when optimizing for NERDSS solutions"""
function generate_population(constraint::Vector{ConstraintRange}, n::Int)
    population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end
#> END


#< DEFAULT FITNESS FUNCTION FACTORY
"""Returns the `fitness function(input)` for the cost function, referencing the ODE problem with closure"""
function make_fitness_function(evalfunc::Function, prob::ODEProblem; fitidx = 1)
    function fitness_function(input::Vector{Float64})
        #? Returns a cost function method that takes in just a vector of parameters/ICs and references the ODE problem 
        return evalfunc(input, prob; idx = fitidx)
    end
    return fitness_function
end
#> END


"""### Callback function that terminates the GA if the number of oscillations exceeds the threshold, and updates the progress bar"""
function ga_callback(trace::Evolutionary.OptimizationTrace, progressbar::Progress, threshold::Int)
    #? Callback function for the GA, updating the progress bar
    num_oscillation = trace[end].metadata["num_oscillatory"]
    if num_oscillation >= threshold 
        finish!(progressbar)
        return true
    else
        next!(progressbar, step = num_oscillation)
        return false
    end
end


#< RUN GENETIC ALGORITHM OPTIMIZATION ##
"""
Runs the genetic algorithm, returning the `result`, and the `record` named tuple
"""
function run_GA(ga_problem::GAProblem, fitnessfunction_factory::Function=make_fitness_function; 
                                            threshold=10000, population_size = 5000, abstol=1e-4, reltol=1e-2, successive_f_tol = 2, iterations=5, parallelization = :thread, fitidx = 1)
    blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    #* Generate the initial population.
    pop = generate_population(ga_problem.constraints, population_size)

    #* Create constraints using the min and max values from param_values.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints.ranges], [constraint.max for constraint in ga_problem.constraints.ranges])

    # *Create Progress bar and callback function
    # ga_progress = Progress(threshold; desc = "GA Progress")
    # callback_func = (trace) -> ga_callback(trace, ga_progress, threshold)

    #* Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=true, show_every=1, parallelization=parallelization)#, callback=callback_func)

    #* Define the range of possible values for each parameter when mutated, and the mutation scalar.
    mutation_scalar = 0.5; mutation_range = fill(mutation_scalar, length(ga_problem.constraints.ranges))

    #* Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
    crossover = TPX, crossoverRate = 1.0, # Two-point crossover event
    mutation  = BGA(mutation_range, 2), mutationRate = 1.0)

    #* Make fitness function. Makes closure of evaluation function and ODE problem
    fitness_function = fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem; fitidx = fitidx)

    #* Run the optimization.
    result = Evolutionary.optimize(fitness_function, [0.0,0.0,0.0], boxconstraints, mthd, pop, opts)
    # return result

    #* Get the individual, fitness, period/amplitude, for each oscillatory individual evaluated
    # record::Vector{NamedTuple{(:ind,:fit,:per,:amp),Tuple{Vector{Float64},Float64, Float64, Float64}}} = reduce(vcat,[gen.metadata["staterecord"] for gen in result.trace[2:end]])
    # num_oscillatory = sum([gen.metadata["num_oscillatory"] for gen in result.trace[2:end]])

    BLAS.set_num_threads(blas_threads)
    return trace_to_df(result)
end
#> END


function trace_to_df(results)
    #* make a dataframe from the trace
    df = DataFrame(ind = [], fit = [], per = [], amp = [])
    for gen in results.trace
        push!(df.ind, gen.metadata["population"]...)
        push!(df.fit, gen.metadata["fitvals"]...)
        push!(df.per, gen.metadata["periods"]...)
        push!(df.amp, gen.metadata["amplitudes"]...)
    end
    return df
end


#< MISCELLANEOUS FUNCTIONS ##
"""Defines logspace function for sampling parameters"""
logrange(start, stop, length) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))

"""Extract solution of a row from the dataframe"""
function extract_solution(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)), tspan = (0.0, 100.0))
    reprob = length(df.ind[row]) > 4 ? remake(prob, p = df.ind[row], tspan = tspan) : remake(prob, u0 = [df.ind[row]; zeros(length(prob.u0) - length(df.ind[row]))], tspan = tspan)
    solve(reprob, Rosenbrock23(), save_idxs = vars)
end


"""Splits ind column into separate columns for each parameter, adds initial conditions for writing DataFrame to CSV"""
function split_dataframe!(df, prob)

    paramsymbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka4,:kb4,:ka7,:kb7,:kcat7,:DF]

    #* Split ind column into separate columns for each parameter
    for (i,param) in enumerate(paramsymbols)
        df[!, param] .= [x[i] for x in df.ind]
    end

    select!(df, Not(:ind))

    #* Add initial conditions
    df.L .= prob.u0[1]
    df.K .= prob.u0[2]
    df.P .= prob.u0[3]
    df.A .= prob.u0[4]
end


"""Find the indices of the inputs in a `NAME` array"""
function find_indices(combination::Vector{String}, NAMES::Vector{String})
    p1idx = findfirst(isequal(combination[1]), NAMES)
    p2idx = findfirst(isequal(combination[2]), NAMES)
    return p1idx, p2idx
end

function find_indices(combination::Vector{Symbol}, NAMES::Vector{String})
    str1 = string(combination[1])
    str2 = string(combination[2])
    p1idx = findfirst(isequal(str1), NAMES)
    p2idx = findfirst(isequal(str2), NAMES)
    return p1idx, p2idx
end

function find_indices(combination::Vector{ConstraintRange}, constraints::Vector{ConstraintRange})::Tuple{Int,Int}
    p1idx = findfirst(x -> x.name == combination[1].name, constraints)
    p2idx = findfirst(x -> x.name == combination[2].name, constraints)
    return p1idx, p2idx
end

"""Triplet version for 3FixedParamCSVMaker"""
function find_indices(param1::String, param2::String, param3::String, constraints::Vector{ConstraintRange})::Tuple{Int,Int,Int}
    p1idx = findfirst(x -> x.name == param1, constraints)
    p2idx = findfirst(x -> x.name == param2, constraints)
    p3idx = findfirst(x -> x.name == param3, constraints)

    return p1idx, p2idx, p3idx
end
#> END
