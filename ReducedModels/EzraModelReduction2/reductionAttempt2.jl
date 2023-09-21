include("reductionSetup.jl")

#Eliminate everything except the four variables that work but in a more computationally
#efficient manner

function setLHSto0!(equationDictionary, Key)
    equationDictionary[Key] = 0 ~ equationDictionary[Key].rhs
end

#Based off of simplification attempt 1
initiallyEliminated = [:dLK, :dLpAK, :dAK, :dAP, :dAKL, :dAPLp,
                        :dLpP, :dLpAP]

originalEqns = deepcopy(rateEquations)

for i in initiallyEliminated
    setLHSto0!(rateEquations, i)
end
rateEquations = deepcopy(originalEqns)

michaelisMenten!(rateEquations, dLK, LK, :dLK)
compareModels(100)

#setLHSto0 does not work as intended :(