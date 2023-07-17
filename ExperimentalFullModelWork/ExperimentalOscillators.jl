using DataFrames
using CSV
include("EvaluationFunctions.jl")
include("Constants.jl")

dfRange = 10.0 .^ (-1:4) #exponential range from 0.1 to 1000
kaRange = 10.0 .^ (-3:1) 
kbRange = 10.0 .^ (-3:3)

function setOutputDir()
    mkdir("osc_dfs")
    cd("osc_dfs")

    ranges = [dfRange, kaRange, kbRange, Lrange, Krange, Prange, Arange]
    rangenames = ["dfRange", "kaRange", "kbRange", "Lrange", "Krange", "Prange", "Arange"]
    open("README", "a") do io
        for i in zip(rangenames, ranges)
            println(io, "$(i[1]): ($(i[2][1]),$(i[2][end])) with $(length(i[2])) elements")
        end
    end      
end

function makeSolutionCSVs()
    setOutputDir()
    for df in dfRange
        println("#################")
        println("df = $(df)")
        println("#################")
        flush(stdout)
        #Initialize new dataframe at each df value
        oscData = DataFrame(DF=Float64[], ka7=Float64[], ka4=Float64[], ka1=Float64[], kb1=Float64[],
                    fit=Float64[],per=Float64[], amp=Float64[],L=Float64[], K=Float64[], P=Float64[], A=Float64[])
        nonoscData = DataFrame(DF=Float64[], ka7=Float64[], ka4=Float64[], ka1=Float64[], kb1=Float64[])
        dfest = df
        numericalErrorCount = 0
        for ka7 in kaRange
            ka7est = ka7
            println("You are at ka7=$(ka7) and df=$(df)")
            flush(stdout)
            for ka4 in kaRange
                ka4est = ka4
                println("ka4 = $(ka4)")
                flush(stdout)
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
                            costPerAmp = adaptiveSolve(prob,u0,shortSpan,longSpan,p; abstol=1e-8,reltol=1e-6)
                            retcode = costPerAmp[1]
                            if retcode < 0.0
                                oscFound+=1
                                if oscFound >= 3
                                    push!(oscData, Dict(:DF => df, :ka7 => ka7, :ka4 => ka4, :ka1 => ka1, :kb1 => kb1,
                                                        :fit => costPerAmp[1], :per => costPerAmp[2], :amp => costPerAmp[3],
                                                        :L => curU0[1],:K => curU0[2],:P => curU0[3],:A => curU0[4]))
                                    break
                                end
                            elseif retcode in [1.0, 1.5]
                                numericalErrorCount += 1
                            end
                        end
                        if oscFound < 3
                            push!(nonoscData, Dict(:DF => df, :ka7 => ka7, :ka4 => ka4, :ka1 => ka1, :kb1 => kb1))
                        end
                        oscFound = 0
                    end
                end
            end
        end
        println("At df = $(df), $(numericalErrorCount) numerical errors occurred.")
        flush(stdout)
        oscFileName = "/Users/ezragreenberg/Julia/ExperimentalFullModelWork/OscCsvs/Osc_df_$(df).csv"
        nonoscFileName = "/Users/ezragreenberg/Julia/ExperimentalFullModelWork/OscCsvs/No_osc_df_$(df).csv"
        CSV.write(oscFileName,oscData)
        CSV.write(nonoscFileName,nonoscData)
    end
end

makeSolutionCSVs()


                            
