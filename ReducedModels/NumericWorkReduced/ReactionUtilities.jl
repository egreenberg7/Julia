"Ranges of constants used in our reaction system"
module ReactionConstants
    export dfRange, kaRange, kbRange, kcatRange, Lrange, Krange, Prange, Arange, randInRange

    const dfRange = (0.1, 10000) #exponential range from 0.1 to 1000
    const kaRange = (0.001, 10)
    const kbRange = (0.001, 1000)
    const kcatRange = (0.001, 1000)
    const Lrange = (0.01,10) 
    const Krange = (0.001,10)
    const Prange = (0.001,1) 
    const Arange = (0.01,10) 

    """
        randInRange(range)
    Returns random number within range, where range is a two element tuple
    """
    function randInRange(range)
        return range[1] + rand() * (range[end] - range[1])
    end
end

"Reactions where only one enzyme can bind to membrane"
module SimpleReactions
    using Catalyst
    using DifferentialEquations
    using .Main.ReactionConstants
    export simplePRxn, simpleKRxn, randomP, randomU

    """
    Reaction with only kinase able to bind to the membrane and adaptor proteins
    only binding to enzymes on the membrane.
    """
    const simplePRxn = @reaction_network simplePRxn begin
        @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka7 kb7 kcat7 df
        @species L(t) K(t) P(t) A(t) Lp(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAKL(t) 
        # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
        # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
        (ka1,kb1), L + K <--> LK # L binding to kinase
        kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
        (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
        (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
        (ka1*df,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
        kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
        (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
        kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    end

    """
    Reaction with only phosphatase able to bind to the membrane and adaptor proteins
    only binding to enzymes on the membrane.
    """
    const simpleKRxn = @reaction_network simpleKRxn begin
        @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 df
        @species L(t) K(t) P(t) A(t) Lp(t) LpA(t) LK(t) LpP(t) LpAP(t) LpAPLp(t) 
        # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
        # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
        (ka1,kb1), L + K <--> LK # L binding to kinase
        kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
        (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
        (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
        kcat7, LpP --> L + P # L dephosphorylation by phosphatase
        (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
        (ka7*df,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
        kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality
    end

    

    """
        randomP()
    Get a random parameter vector for a simplified model
    """
    function randomP()
        kAs = [randInRange(kaRange) for i in 1:4]
        kBs = [randInRange(kbRange) for i in 1:4]
        kCats = [randInRange(kcatRange) for i in 1:2]
        df = randInRange(dfRange)
        return [kAs[1], kBs[1], kCats[1], kAs[2], kBs[2], kAs[3], kBs[3], kAs[4], kBs[4], kCats[2], df]
    end

    """
        randomU()
    Get a random concentration vector for a simplified model
    """
    function randomU()
        return vcat([randInRange(i) for i in [Lrange, Krange, Prange, Arange]], zeros(6))
    end
end

"Reactinos without membrane localization"
module NoMembraneReactions
    using Catalyst
    using DifferentialEquations
    using .Main.ReactionConstants
    export simplePRxn, simpleKRxn, randomP, randomU

    """
    React with neither kinase nor phosphatase able to bind to the membrane.
    """
    const noMembraneReactions = @reaction_network simplePRxn begin
        @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka7 kb7 kcat7 df
        @species L(t) K(t) P(t) Lp(t) LK(t) LpP(t)
        # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
        # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
        (ka1,kb1), L + K <--> LK # L binding to kinase
        kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
        (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
        kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    end

end

"Toy reactions for testing"
module ToyReactions
    using Catalyst
    const enzymeReaction = @reaction_network enzymeReaction begin
        @parameters ka kb kcat
        @species S(t) E(t) ES(t) P(t)
        (ka, kb), S + E <--> ES
        kcat, ES --> P
    end
end

"Module with full reaction and random concentration and parameter generators"
module FullReaction
    using Catalyst
    using .Main.ReactionConstants

    export fullRxn, randomP, randomU

    const fullRxn = @reaction_network fullrn begin
        @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
        @species L(t) K(t) P(t) A(t) Lp(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
        # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
        # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
        (ka1,kb1), L + K <--> LK # L binding to kinase
        kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
        (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
        (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
        (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
        kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
        (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
        kcat7, LpP --> L + P # L dephosphorylation by phosphatase
        (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
        (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
        kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

        #previously excluded reactions, all possible combinations possible in vitro
        (ka2,kb2), Lp + AK <--> LpAK
        (ka2*y,kb2), Lp + AKL <--> LpAKL
        (ka2,kb2), Lp + AP <--> LpAP
        (ka2*y,kb2), Lp + APLp <--> LpAPLp
        (ka3,kb3), A + K <--> AK
        (ka4,kb4), A + P <--> AP
        (ka3,kb3), A + LK <--> AKL
        (ka4,kb4), A + LpP <--> APLp
        (ka3*y,kb3), LpA + LK <--> LpAKL
        (ka4*y,kb4), LpA + LpP <--> LpAPLp
        (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
        kcat1, AKL --> Lp + AK #phosphorylation of lipid
        (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
        kcat7, APLp --> L + AP #dephosphorylation of lipid
    end

    """
        randomP()
    Get a random parameter vector for the full model
    """
    function randomP()
        kAs = [randInRange(kaRange) for i in 1:5]
        kBs = [randInRange(kbRange) for i in 1:5]
        kCats = [randInRange(kcatRange) for i in 1:2]
        df = randInRange(dfRange)
        return [kAs[1], kBs[1], kCats[1], kAs[2], kBs[2], kAs[3], kBs[3], kAs[4], kBs[4], kAs[5], kBs[5], kCats[2], df]
    end

    """
        randomU()
    Get a random concentration vector for full model
    """
    function randomU()
        return vcat([randInRange(i) for i in [Lrange, Krange, Prange, Arange]], zeros(12))
    end
end