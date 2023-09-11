using Symbolics
using DifferentialEquations
using Plots
using Catalyst
using CSV
using DataFrames

const fullrn = @reaction_network fullrn begin
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
Function representing differential equations of full model
"""
function fullmodel_ode!(du, u, p, t)
    L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
    ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y = p #parameters
    du[1] = kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK #* L
    du[2] = kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA #* K
    du[3] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P #* P
    du[4] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P #* A
    du[5] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP #* Lp
    du[6] = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*y*LK*LpA - ka4*y*LpA*LpP #* LpA
    du[7] = kb3*AKL + kb3*LpAKL + ka1*K*L - kb1*LK - kcat1*LK - ka3*A*LK - ka3*y*LK*LpA #* LK
    du[8] = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*y*LpA*LpP #* LpP
    du[9] = kb1*LpAKL + kcat1*LpAKL + ka2*AK*Lp + ka3*K*LpA - kb2*LpAK - kb3*LpAK - ka1*y*L*LpAK #* LpAK
    du[10] = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*y*Lp*LpAP #* LpAP
    du[11] = ka1*y*L*LpAK + ka2*y*AKL*Lp + ka3*y*LK*LpA - kb1*LpAKL - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL #* LpAKL
    du[12] = ka2*y*APLp*Lp + ka7*y*Lp*LpAP + ka4*y*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp #* LpAPLp
    du[13] = kb1*AKL + kb2*LpAK + kcat1*AKL + ka3*A*K - kb3*AK - ka1*AK*L - ka2*AK*Lp #* AK
    du[14] = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp #* AP
    du[15] = kb2*LpAKL + ka3*A*LK + ka1*AK*L - kb1*AKL - kb3*AKL - kcat1*AKL - ka2*y*AKL*Lp #* AKL
    du[16] = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*y*APLp*Lp #* APLp
    nothing
end

@variables ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 df
@variables L K P A Lp LpA LK LpP LpAK LpAP LpAKL LpAPLp AK AP AKL APLp
@variables dL dK dP dA dLp dLpA dLK dLpP dLpAK dLpAP dLpAKL dLpAPLp dAK dAP dAKL dAPLp
@variables Ltot, Atot, Ptot, Ktot

"""
Dictionary of the Julia symbolic objects representing our differential equations.
We will be manipulating this to reduce the model.
"""
rateEquations  = Dict([:dL => (dL ~ kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*df*L*LpAK) 
                :dK => (dK ~ kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA)
                :dP => (dP ~ kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P) 
                :dA => (dA ~ kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P) 
                :dLp => (dLp ~ kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*df*AKL*Lp - ka2*df*APLp*Lp - ka7*df*Lp*LpAP) 
                :dLpA => (dLpA ~ kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*df*LK*LpA - ka4*df*LpA*LpP)
                :dLK => (dLK ~ kb3*AKL + kb3*LpAKL + ka1*K*L - kb1*LK - kcat1*LK - ka3*A*LK - ka3*df*LK*LpA)
                :dLpP => (dLpP ~ kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*df*LpA*LpP)
                :dLpAK => (dLpAK ~ kb1*LpAKL + kcat1*LpAKL + ka2*AK*Lp + ka3*K*LpA - kb2*LpAK - kb3*LpAK - ka1*df*L*LpAK)
                :dLpAP => (dLpAP ~ kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*df*Lp*LpAP)
                :dLpAKL => (dLpAKL ~ ka1*df*L*LpAK + ka2*df*AKL*Lp + ka3*df*LK*LpA - kb1*LpAKL - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL)
                :dLpAPLp => (dLpAPLp ~ ka2*df*APLp*Lp + ka7*df*Lp*LpAP + ka4*df*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp)
                :dAK => (dAK ~ kb1*AKL + kb2*LpAK + kcat1*AKL + ka3*A*K - kb3*AK - ka1*AK*L - ka2*AK*Lp)
                :dAP => (dAP ~ kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp)
                :dAKL => (dAKL ~ kb2*LpAKL + ka3*A*LK + ka1*AK*L - kb1*AKL - kb3*AKL - kcat1*AKL - ka2*df*AKL*Lp)
                :dAPLp => (dAPLp ~ kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*df*APLp*Lp)]
)

"""
Function to allow us to evaluate our reduced models through the DifferentialEquations packages.
"""
function symbolicmodel_ode!(du, u, p, t; rateEquations = rateEquations)
    theL, theK, theP, theA, theLp, theLpA, theLK, theLpP, theLpAK, theLpAP, theLpAKL, theLpAPLp, theAK, theAP, theAKL, theAPLp = u #initial conditions
    theka1, thekb1, thekcat1, theka2, thekb2, theka3, thekb3, theka4, thekb4, theka7, thekb7, thekcat7, thedf = p #parameters

    function getDerivative(symbol, u = u, p = p; rateEquations = rateEquations)
        if typeof(rateEquations[symbol].lhs == 0) == Bool && rateEquations[symbol].lhs == 0
             return 0.0
        else
            expr = (rateEquations[symbol].rhs)
            val =  Symbolics.value(substitute(expr, 
                Dict(L=>theL, K=>theK , P=>theP , A=>theA, Lp=>theLp, LpA=>theLpA,LK=>u[7], LpP=>u[8], LpAK=>u[9], LpAP=>u[10],
                    LpAKL=>u[11],  LpAPLp=>u[12], AK=>u[13], AP=>u[14], AKL=>u[15], APLp=>u[16], 
                    ka1=>p[1], kb1=>p[2], kcat1=>p[3], ka2=>p[4], kb2=>p[5], ka3=>p[6], kb3=>p[7], ka4=>p[8], kb4=>p[9],
                    ka7=>p[10], kb7=>p[11], kcat7=>p[12], df=>p[13])))
            return val
        end
    end

    du[1] = getDerivative(:dL)
    du[2] = getDerivative(:dK)
    du[3] = getDerivative(:dP)
    du[4] = getDerivative(:dA)
    du[5] = getDerivative(:dLp)
    du[6] = getDerivative(:dLpA)
    du[7] = getDerivative(:dLK)
    du[8] = getDerivative(:dLpP)
    du[9] = getDerivative(:dLpAK)
    du[10] = getDerivative(:dLpAP)
    du[11] = getDerivative(:dLpAKL)
    du[12] = getDerivative(:dLpAPLp)
    du[13] = getDerivative(:dAK)
    du[14] = getDerivative(:dAP)
    du[15] = getDerivative(:dAKL)
    du[16] = getDerivative(:dAPLp)
    nothing
end

"""
Oscillatory datapoints from Jonathan that we can use to validate the reduced models
"""
datapoints = DataFrame(CSV.File("ReducedModels/EzraModelReduction2/ezra_optimized_params.csv"))

"""
Function to extract parameters from dataframe
"""
function getP(df, row)
    return [df.ka1[row], df.kb1[row], df.kcat1[row], df.ka2[row], df.kb2[row], 
            df.ka3[row], df.kb3[row], df.ka4[row], df.kb4[row], df.ka7[row], df.kb7[row],
            df.kcat7[row], df.DF[row]]
end

"""
Function to extract concentrations from dataframe
"""
function getU(df, row)
    u = zeros(16)
    u[1] = df.L[row]
    u[2] = df.K[row]
    u[3] = df.P[row]
    u[4] = df.A[row]
    return u
end

"""
Function to generate plots of the solved equations from both the full model and
the reduced model for comparison
"""
function compareModels(row, model = symbolicmodel_ode!)
    p = getP(datapoints, row)
    u = getU(datapoints, row)
    actualProb = ODEProblem(fullmodel_ode!, u, (0,500), p)
    reducedProb = ODEProblem(model, u, (0,500), p)
    actualSol = solve(actualProb, save_idxs=1)
    modelSol = solve(reducedProb, save_idxs=1)
    plot(actualSol, label="Full model")
    plot!(modelSol, label="Reduced model")
end

"""
Function to reduce our model by substituting in 0 for the derivative of one of our variables
"""
function michaelisMenten!(rateEquations, eliminatedDifferential, eliminatedVar, eliminatedVarSymbol)
    for i in 1:2
        substitute(rateEquations[eliminatedVarSymbol], eliminatedDifferential=>0)
        LKsub = Symbolics.solve_for(rateEquations[eliminatedVarSymbol], eliminatedVar)
        map!(eq->(substitute(eq, (eliminatedVar=>LKsub, eliminatedDifferential=>0))), values(rateEquations))
    end
end

#We now making a series of Michaelis-Menten approximations

#First we eliminate LK
michaelisMenten!(rateEquations, dLK, LK, :dLK)
#compareModels(30)

#Now let's get rid of LpAK
michaelisMenten!(rateEquations, dLpAK, LpAK, :dLpAK)
#compareModels(30)

#Getting rid of LpAKL next ends the oscillations :( 
#Getting rid of AK next leads to damped oscillations for row 30, not for 100
#Also for some reason, I have to call michaelisMenten! twice to avoid errors?

michaelisMenten!(rateEquations, dAK, AK, :dAK)

#Still working for row 100, no oscillations for row 30
michaelisMenten!(rateEquations, dAP, AP, :dAP)
#compareModels(30)

#Row 100 continues to succeed
michaelisMenten!(rateEquations, dAKL, AKL, :dAKL)
#compareModels(100)

#Row 100 continues to work, the period has now decreased
michaelisMenten!(rateEquations, dAPLp, APLp, :dAPLp)
#compareModels(100)

#We graph all the variables of the reduced model 
#to inform our choice of the next variable to eliminate
modelprob = ODEProblem(symbolicmodel_ode!, getU(datapoints, 100), (0,150), getP(datapoints, 100))
sol = solve(modelprob)
plot(sol)

#LpAKL kills the model
#michaelisMenten!(rateEquations, dLpAKL, LpAKL, :dLpAKL)
#compareModels(100)

#Still workin good after LpP
michaelisMenten!(rateEquations, dLpP, LpP, :dLpP)
compareModels(100)

#Still working after LpAP eliminated!!!
michaelisMenten!(rateEquations, dLpAP, LpAP, :dLpAP)
compareModels(100)

#I get an error if i now try to get rid of LpA, though I can solve for LpAPLp this way
michaelisMenten!(rateEquations, dLpA, LpA, :dLpA)

michaelisMenten!(rateEquations, dP, P, :dP)


odes = convert(ODESystem, fullrn)
symbolicsODEs = ModelingToolkit.get_eqs(odes)

substitute(symbolicsODEs, Differential(t)(LK(t))=>0)

#Random commented out code
begin
    #delete!(rateEquations, :dLK)
    #x = Symbolics.build_function(rateEquations[:dL].rhs)

    #L=>L, K=>K, P=>P, A=>A, Lp=>Lp, LpA=>LpA,LK=>LK, LpP=>LpP, LpAK=>LpAK, LpAP=>LpAP,
    #       LpAKL=>LpAKL,  LpAPLp=>LpAPLp, AK=>AK, AP=>AK, AP=>AP, AKL=>AKL, APLp=>APLp
    #    aL, aK, aP, aA, aLp, aLpA, aLK, aLpP, aLpAK, aLpAP, aLpAKL, aLpAPLp, aAK, aAP, aAKL, aAPLp = u 

    # Dict(L=>u[1], K=>u[2], P=>u[3], u[4]=>A, Lp=>u[5], LpA=>u[6],LK=>u[7], LpP=>u[8], LpAK=>u[9], LpAP=>u[10],
    #         LpAKL=>u[11],  LpAPLp=>u[12], AK=>u[13], AP=>u[14], AKL=>u[15], APLp=>u[16], 
    #         ka1=>ka1, kb1=>kb1, kcat1=>kcat1, ka2=>ka2, kb2=>kb2, ka3=>ka3, kb3=>kb3, ka4=>ka4, kb4=>kb4,
    #         ka7=>ka7, kb7=>kb7, kcat7=>kcat7, df=>df)))

    #Reaction System
    # D(L) = kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK #* L
    # D(K) = kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA #* K
    # D(P) = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P #* P
    # D(A) = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P #* A
    # D(Lp) = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP #* Lp
    # D(LpA) = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*y*LK*LpA - ka4*y*LpA*LpP #* LpA
    # D(LK) = kb3*AKL + kb3*LpAKL + ka1*K*L - kb1*LK - kcat1*LK - ka3*A*LK - ka3*y*LK*LpA #* LK
    # D(LpP) = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*y*LpA*LpP #* LpP
    # D(LpAK) = kb1*LpAKL + kcat1*LpAKL + ka2*AK*Lp + ka3*K*LpA - kb2*LpAK - kb3*LpAK - ka1*y*L*LpAK #* LpAK
    # D(LpAP) = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*y*Lp*LpAP #* LpAP
    # D(LpAKL) = ka1*y*L*LpAK + ka2*y*AKL*Lp + ka3*y*LK*LpA - kb1*LpAKL - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL #* LpAKL
    # D(LpAPLp) = ka2*y*APLp*Lp + ka7*y*Lp*LpAP + ka4*y*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp #* LpAPLp
    # D(AK) = kb1*AKL + kb2*LpAK + kcat1*AKL + ka3*A*K - kb3*AK - ka1*AK*L - ka2*AK*Lp #* AK
    # D(AP) = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp #* AP
    # D(AKL) = kb2*LpAKL + ka3*A*LK + ka1*AK*L - kb1*AKL - kb3*AKL - kcat1*AKL - ka2*y*AKL*Lp #* AKL
    # D(APLp) = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*y*APLp*Lp #* APLp


    #Mass Conservation
    # ConservationEquations = Dict(
    #     [:L=> (Ltot ~ L + Lp + LpA + LK + LpP + LpAK + LpAP + 2 * LpAKL +  2 * LpAPLp + AKL + APLp)
    #     :P=> (Ptot ~ P + LpP + LpAP + LpAPLp + AP + APLp)
    #     :K=> (Ktot ~ K + LK + LpAK + LpAKL + AK + AKL)
    #     :A=> (Atot ~ A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)]
    # )

    #We start by substituting in conservation formulas
    # LpAPLpsub = Symbolics.solve_for(ConservationEquations[:P], LpAPLp)
    # Asub = Symbolics.solve_for(ConservationEquations[:A], A)
    # Ksub = Symbolics.solve_for(ConservationEquations[:K], K)
    # LpAKLsub = Symbolics.solve_for(ConservationEquations[:L], LpAKL)

    #map!(eq->simplify(substitute(eq, (LpAPLp=>LpAPLpsub, A=>Asub, K=>Ksub, LpAKL=>LpAKLsub))), values(rateEquations))

    # for i in [:dLpAPLp, :dA, :dK, :dLpAKL]
    #     delete!(rateEquations, i)
    # end
end