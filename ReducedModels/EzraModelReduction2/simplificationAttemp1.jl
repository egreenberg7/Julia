include("reductionSetup.jl")

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

#I now save this intermediate reduced model so I can try out different things on it.
RateEq4Var = deepcopy(rateEquations)

#I no longer eliminate the variable analytically due to the solver, but I set another
#variable equal to 0. We try P. Successful! 
#rateEquations[:dP] = 0 ~ RateEq4Var[:dP].rhs
#This might not be accurate acatually, do let's do MM again.
michaelisMenten!(rateEquations, dP, P, :dP)
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