module DeficiencyCalculator
    using Catalyst
    using LinearAlgebra
    using Graphs

    export deficiency

    function deficiency(rxn) 
        return  numComplexesMinusLinkageClasses(rxn) - rank(netstoichmat(rxn))
    end

    function numComplexesMinusLinkageClasses(rxn)
        complexIncMatrix = reactioncomplexes(rxn)[2]
        numComplexes = size(complexIncMatrix)[1]
        diGraph = SimpleDiGraph(complexIncMatrix)
        undiGraph = SimpleGraph(diGraph)
        linkageClasses = connected_components(undiGraph)
        return numComplexes - length(linkageClasses)
    end
end