using DataFrames
using CSV
include("EvaluationFunctions.jl")
include("FinalConstants.jl")


#This file takes in the value of df and the minimum and maximum of the ka7 range you are exploring.
dfRange = [parse(Float64, ARGS[1])]
#ka7 inputs can go from 1 to 9, consistent with kaRange
ka7minExponent, ka7maxExponent = parse(Int64, ARGS[2]), parse(Int64, ARGS[3])

"""
Make directory where output CSVs where be stored. Also makes README 
file detialing ranges used on the program run and tests that CSVs can successfully
be created.
"""
function setOutputDir(df=dfRange[1], ka7=ka7minExponent)
    mkdir("$(df)DF_$(ka7)ka7min")
    cd("$(df)DF_$(ka7)ka7min")
    mkdir("AllData")
    ranges = [dfRange, kaRange[ka7minExponent:ka7maxExponent], kaRange, kbRange, Lrange, Krange, Prange, Arange]
    rangenames = ["dfRange", "ka7Range", "kaRange", "kbRange", "Lrange", "Krange", "Prange", "Arange"]
    open("README", "a") do io
        for i in zip(rangenames, ranges)
            println(io, "$(i[1]): ($(i[2][1]),$(i[2][end])) with $(length(i[2])) elements")
        end
    end
    testdf = DataFrame(Dict("X" => 0.0))
    CSV.write("AllData/test.csv", testdf)
end

function printProgressMessage(message)
    open("progress.txt", "a") do file
        println(file, message)
    end
end


"""
Loop over all ranges and evaluate if oscillatory.
"""
function makeSolutionCSVs(dfRange=dfRange, ka7minExponent=ka7minExponent, ka7maxExponent=ka7maxExponent, kaRange=kaRange, kbRange=kbRange, Lrange=Lrange, Prange=Prange, Krange=Krange, Arange=Arange)
    setOutputDir()
    numRateCombos = length(dfRange) * (ka7maxExponent - ka7minExponent + 1) * length(kaRange) * length(kaRange) * length(kbRange)
    for dfest in dfRange
        #Initialize dictionary to store data in, will be converted to dataframe at termination
        data = Dict(:df => zeros(Float64, numRateCombos), :ka7 => zeros(Float64, numRateCombos),
            :ka4 => zeros(Float64, numRateCombos), :ka1 => zeros(Float64, numRateCombos), :kb1 => zeros(Float64, numRateCombos),
            :fit => zeros(Float64, numRateCombos), :per => zeros(Float64, numRateCombos), :amp => zeros(Float64, numRateCombos),
            :L => zeros(Float64, numRateCombos), :K => zeros(Float64, numRateCombos), :P => zeros(Float64, numRateCombos), :A => zeros(Float64, numRateCombos),
            :oscFound => zeros(Int64, numRateCombos))            
        printProgressMessage("Initalizing: df = $(dfest)")
        #Initialize new dataframe at each df value
        count = 1
        numericalErrorCount = 0
        for ka7est in kaRange[ka7minExponent:ka7maxExponent]
            printProgressMessage("You are at ka7=$(ka7est).")
            for ka4est in kaRange
                printProgressMessage("ka4 = $(ka4est)")
                for ka1est in kaRange
                    for kb1est in kbRange
                        #Set up parameter values
                        psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est,
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
                        p = [x[2] for x in psym]
                        #Count how many oscillatory solutions
                        oscFound = 0
                        for i in 1:numU0combos
                            curU0 = u0combos[i, :]
                            u0[1:4] = curU0[1:4]
                            costPerAmp = adaptiveSolve(prob, u0, shortSpan, longSpan, p; abstol=1e-8, reltol=1e-6)
                            retcode = costPerAmp[1]
                            if retcode < 0.0
                                oscFound += 1
                                if oscFound >= 3
                                    data[:df][count] = dfest
                                    data[:ka7][count] = ka7est
                                    data[:ka1][count] = ka1est
                                    data[:ka4][count] = ka4est
                                    data[:kb1][count] = kb1est
                                    data[:fit][count] = retcode
                                    data[:per][count] = costPerAmp[2]
                                    data[:amp][count] = costPerAmp[3]
                                    data[:L][count] = curU0[1]
                                    data[:K][count] = curU0[2]
                                    data[:P][count] = curU0[3]
                                    data[:A][count] = curU0[4]
                                    data[:oscFound][count] = oscFound
                                    count += 1
                                    break
                                end
                            elseif retcode in [1.0, 1.5]
                                numericalErrorCount += 1
                            end
                        end
                        if oscFound < 3
                            data[:df][count] = dfest
                            data[:ka7][count] = ka7est
                            data[:ka1][count] = ka1est
                            data[:ka4][count] = ka4est
                            data[:kb1][count] = kb1est
                            data[:fit][count] = retcode
                            data[:oscFound][count] = oscFound
                            count += 1
                        end
                        oscFound = 0
                    end
                end
            end
            datadf = DataFrame(data)
            allDataName = "AllData/df_($dfest)_ka7_$(ka7est).csv"
            CSV.write(allDataName, datadf)
        end
        printProgressMessage("At df = $(dfest), $(numericalErrorCount) numerical errors occurred.")
    end
end

makeSolutionCSVs()