include("OutputHandling.jl")
using StatsPlots

oscdata = DataFrame(CSV.File("/Users/ezragreenberg/JLab/Julia/ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))

# """
# Double check if my oscillatory data is actually oscillatory by evaluating
# over 5000 seconds. Returns array of 1s for oscillatory solutions and 0s for
# non-oscillatory solutions
# """
# function checkLongerTimeSpan(longTspan)
#     numOscillatory = size(oscdata)[1]
#     trulyOscillatory = zeros(numOscillatory)
#     for i in 1:numOscillatory
#         longsol = entryToSol(oscdata,i; tspan = longTspan)
#         if finalClassifier(longsol, longTspan)[1] < 0
#             trulyOscillatory[i] = 1
#         end
#         if i % 50 == 0
#             println("Checked $i solutions")
#         end
#     end
#     return trulyOscillatory
# end

# trulyOscillatory = checkLongerTimeSpan(5000)

# trulyOscData = oscdata[trulyOscillatory .== 1.0, :]
# CSV.write("ExperimentalFullModelWork/TrulyOscillatoryData.csv", trulyOscData)

trulyOscData = CSV.read("ExperimentalFullModelWork/TrulyOscillatoryData.csv",DataFrame)

function manipulateDFentryToSol(row, df; dframe = trulyOscData, tspan = shortSpan, save_idxs = 1)
    currow = dframe[row,:]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    ka1est = currow[:ka1]
    kb1est = currow[:kb1]
    ka4est = currow[:ka4]
    ka7est = currow[:ka7]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => df]
    p = [x[2] for x in psym]
    return solve(remake(prob, u0=u0, tspan=(0.0, tspan), p=p), Rodas4(), abstol=1e-8, reltol=1e-12, saveat=0.1, save_idxs=save_idxs, maxiters=200 * tspan, verbose=false) 
end

"""
Calculate the minimum df at which oscillations for a given oscillatory parameter set.
"""
function getMinDf(row; dframe = trulyOscData, tspan = 5000)
    previousOscillatory = true
    increment = -1000
    curDF = dframe[row, :df]
    previousDF = curDF

    #Continueally jump up or down increments until overshooting, then narrow jumps
    while abs(increment) > 1
        curSol = manipulateDFentryToSol(row, curDF; tspan = tspan)
        curOscillatory = finalClassifier(curSol, tspan)[1] < 0
        if (curOscillatory && !previousOscillatory) || (!curOscillatory && previousOscillatory)
            println(curDF)
            previousDF = curDF
            increment /= -10
            previousOscillatory = curOscillatory
        end
        curDF += increment
        println(curDF)
    end
    println("Made it!")

    #Only incrementing by 1 until we reach oscillatory solution
    while !previousOscillatory
        curSol = manipulateDFentryToSol(row, curDF; tspan = tspan)
        if finalClassifier(curSol, tspan)[1] < 0
            previousOscillatory = true
        else
            curDF += increment
        end
    end
    return curDF
end

function getAllMinDfs!(dFrame=trulyOscData; savefile = true, filename = "ExperimentalFullModelWork/paramsWithMinDF.csv")
    numRows = size(dFrame)[1]
    minDFs = zeros(numRows)
    for i in 1:numRows
        minDFs[i] = getMinDf(i)
    end
    dFrame.minimumDF = minDFs
    if savefile
        CSV.write(filename, dFrame)
    end
end

#getAllMinDfs!() 

moreData = CSV.read("ExperimentalFullModelWork/paramsWithMinDF.csv", DataFrame)
@df moreData scatter(log10.(:L), log10.(:minimumDF))

function scatterParams(dFrame, param1, param2)
    scatter(log10.(dFrame[:,param1]), log10.(dFrame[:,param2]), xlabel=param1, ylabel=param2, title="Truly Log10")
end

function pairwiseScatters(dFrame)
    numParams = size(names(dFrame))[1]
    params = propertynames(dFrame)
    deleteat!(params, numParams - 2)
    for i in 1:(numParams - 2)
        for j in (i+1):numParams - 1
            display(scatterParams(dFrame, params[i], params[j]))
        end
    end
end

pairwiseScatters(moreData)

evenMoreData = moreData
evenMoreData.NormalizedAmp = evenMoreData.amp ./ evenMoreData.L
