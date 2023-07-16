using BenchmarkTools
using CSV
using DataFrames
include("EvaluationFunctions.jl")
include("Constants.jl")
function testClassifier(numIterations; saveCSV = true, outputdirectory = "/Users/ezragreenberg/Julia/Someplots/",filename="mytest1.csv", abstol=1e-8, reltol=1e-12)
    #df = DataFrame(u0=Vector{Float64}[], p=Vector{Float64}[], retcode=Float64[], per=Float64[], amp=Float64[])
    #CSV of array gets read in as string, will have to generalize this later
    df = DataFrame(L=Float64[],K=Float64[],P=Float64[],A=Float64[], retcode=Float64[], per=Float64[], amp=Float64[])
    for i in 1:numIterations
        u0[1] = rand(Random.seed!(i),Distributions.LogUniform(0.01, 100)) #Lp #1 for L, 2 for Lp
        u0[2] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.001, 100)) #K
        u0[3] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.01, 100)) #P
        u0[4] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.001, 100)) #A
        retcode = adaptiveSolve(prob, u0, shortspan, longspan, p; abstol = abstol, reltol = reltol)
        push!(df, Dict(:L=>u0[1],:K=>u0[2],:P=>u0[3],:A=>u0[4], :retcode => retcode[1], :per => retcode[2], :amp =>retcode[3]))
        if i % 500 == 0 
            println("$i iterations completed")
        end
    end
    if saveCSV
        CSV.write(outputdirectory*filename, df)
    end
    return df
end
@btime testClassifier(500)