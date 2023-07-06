using DifferentialEquations
using Plots
using Catalyst
using DataFrames
using CSV
using Statistics
using Peaks
using FFTW
include("EvaluationFunctions.jl")
include("Constants.jl")

oscData = DataFrame(DF=Float64[], ka7=Float64[], ka4=Float64[], ka1=Float64[], kb1=Float64[],
                    fit=Float64[],per=Float64[], amp=Float64[],L=Float64[], K=Float64[], P=Float64[], A=Float64[])
nonoscData = DataFrame(DF=Float64[], ka7=Float64[], ka4=Float64[], ka1=Float64[], kb1=Float64[])

for df in dfRange
    dfest = df
    for ka7 in kaRange
        ka7est = ka7
        for ka4 in kaRange
            ka4est = ka4
            for ka1 in kaRange
                ka1est = ka1
                for kb1 in kbRange
                    kb1est = kb1
                    #Set up parameter values
                    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
                    p = [x[2] for x in psym]
                    #Count how many oscillatory solutions
                    oscFound = 0
                    for i in 1:numU0combos
                        curU0 = u0combos[i,:]
                        u0[1:4] = curU0[1:4]
                        currentSol = solve(remake(prob, u0=u0, p=p), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
                        costPerAmp = eval_fitness_catcherrors(currentSol)
                        if(costPerAmp[1] < 0)
                            oscFound+=1
                            if oscFound >= 3
                                push!(oscData, Dict(:DF => df, :ka7 => ka7, :ka4 => ka4, :ka1 => ka1, :kb1 => kb1,
                                                    :fit => costPerAmp[1], :per => costPerAmp[2], :amp => costPerAmp[3],
                                                    :L => curU0[1],:K => curU0[2],:P => curU0[3],:A => curU0[4]))
                                break
                            end
                        end
                    end
                    if oscFound < 3
                        push!(nonoscData, Dict(:DF => df, :ka7 => ka7, :ka4 => ka4, :ka1 => ka1, :kb1 => kb1))
                    end
                    oscFound = 0
                end
            end
        end
        println("You are at ka7=$(ka7) and df=$(df)")
    end
    oscFileName = "OscCsvs/Osc_df_$(df).csv"
    nonoscFileName = "OscCsvs/No_osc_df_$(df).csv"
    CSV.write(oscFileName,oscData)
    CSV.write(nonoscFileName,nonoscData)
end



                            
