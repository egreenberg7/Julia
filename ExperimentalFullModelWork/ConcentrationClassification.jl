#File for functions based around classifying solutions over range of
#initial concentrations rather than rate constants.

#include("NoNegConstants.jl")
#include("EvaluationFunctions.jl")
using CSV

"""
Function to generate dictionary of concentration ranges to iterate over.
- `logStepSize` scaling factor between elements in the range
- `XRange` optional parameters to extract min and max from
"""
function getu0Ranges(logStepSize=0.2; Lrange=Lrange, Krange=Krange,Prange=Prange,Arange=Arange)
    myRanges = Dict()
    keys = [:L, :K, :P, :A]
    ranges = [Lrange, Prange, Krange, Arange]
    for i in zip(keys, ranges)
        logRange = range(log10(i[2][1]), log10(i[2][end]); step=logStepSize)
        myRanges[i[1]] = 10 .^ logRange
    end
    return myRanges
end
#TODO DOUBLE CHECK RANGES IN NoNegCONSTANTS.jl FOR CONCNETRATIONS!!!

"""
Generate a dataframe of classifications of solutions to various initial
conditions; classifications based on `singleSolve` in EvaluationFunctions.jl 
- `p` rate constants being used
- `u0RangeDict` Dictionary of initial conditions to iterate over from `getu0Ranges`
"""
function classifyConcentrations(p, u0RangeDict; u0 = zeros(16), savedata = true, filename="ConcentrationOutput")
    if savedata
        mkdir(filename)
        cd(filename)
        open("README", "a") do io
            println(io, "P = $p")
            println(io, "u0Ranges:")
            println(io, u0RangeDict)
        end
    end
    numU0 = 1
    for rng in values(u0RangeDict) numU0 *= length(rng) end
    data = Dict(:L=>zeros(numU0), :K=>zeros(numU0), :P=>zeros(numU0), :A=>zeros(numU0), 
                :retcode=>zeros(numU0), :per=>zeros(numU0), :amp=>zeros(numU0))
    count = 1
    for L in u0RangeDict[:L]
        u0[1] = L
        for K in u0RangeDict[:K]
            u0[2] = K
            for P in u0RangeDict[:P]
                u0[3] = P
                for A in u0RangeDict[:A]
                    u0[4] = A
                    retcode, per, amp = singleSolve(prob, u0, longSpan, p)
                    data[:L][count] = L
                    data[:K][count] = K
                    data[:P][count] = P
                    data[:A][count] = A
                    data[:retcode][count] = retcode
                    data[:per][count] = per
                    data[:amp][count] = amp
                    count += 1
                end
            end
        end
    end
    df=DataFrame(data)
    if savedata
        CSV.write("$(filename).csv", df)
        cd("..")
    end
    return df
end

function getConcentrationCombos(u0RangeDict; u0 = zeros(16), savedata = false, filename="ConcentrationOutput")
    numU0 = 1
    for rng in values(u0RangeDict) numU0 *= length(rng) end
    data = Dict(:L=>zeros(numU0), :K=>zeros(numU0), :P=>zeros(numU0), :A=>zeros(numU0))
    count = 1
    for L in u0RangeDict[:L]
        u0[1] = L
        for K in u0RangeDict[:K]
            u0[2] = K
            for P in u0RangeDict[:P]
                u0[3] = P
                for A in u0RangeDict[:A]
                    u0[4] = A
                    data[:L][count] = L
                    data[:K][count] = K
                    data[:P][count] = P
                    data[:A][count] = A
                    count += 1
                end
            end
        end
    end
    df=DataFrame(data)
    if savedata
        CSV.write("$(filename).csv", df)
        cd("..")
    end
    return df
end

