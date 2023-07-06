using DifferentialEquations
using Plots
using Catalyst
using DataFrames
using CSV
using Statistics
using Peaks
using FFTW
include("Constants.jl")
include("EvaluationFunctions.jl")

osc_data = DataFrame(CSV.File("ExperimentalFullModelWork/OscCsvs/Osc_df_1.0.csv"))

function entryToSol(df, row)
    currow = df[row,:]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    ka1est = currow[:ka1]
    kb1est = currow[:kb1]
    ka4est = currow[:ka4]
    ka7est = currow[:ka7]
    dfest = currow[:DF]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
    p = [x[2] for x in psym]
    return solve(remake(prob, u0=u0, p=p), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
end

for i in 1:size(osc_data)[1]
    cursol = entryToSol(osc_data, i)
    solPlot = plot(cursol, title="L vs t")
    fftSolPlot = plot((2:length(cursol.t)/2 + 1),(abs.(rfft(cursol.u)))[2:end], title="Fourier Transform")
    plot(solPlot,fftSolPlot, layout=(2,1))
    if(isSteady(cursol))
        println("We have a problem")
    end
    savefig("Someplots/$i.png")
end