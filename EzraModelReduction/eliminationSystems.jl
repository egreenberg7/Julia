using DifferentialEquations
using Plots

#=THESE CONSTRAINTS WERE GENERATED BY MATHEMATICA; I isolated dL from the the first 3 by simple division
    Each of these was used in the corresponding elimination system
    1. dL = (AKL^2*kb3*Km1 - AK*AKL*kb3*L - AKL*K*kb3*L + AKL*K*kcat1*L - APLp*K*kcat7*L - A*AKL*ka3*Km1*LK + AKL*kb3*Km1*LK - AKL*kcat1*Km1*LK + APLp*kcat7*Km1*LK + A*AK*ka3*L*LK + A*K*ka3*L*LK + K*kcat1*L*LK - A*ka3*Km1*LK^2 - kcat1*Km1*LK^2 - AKL*DF*ka3*Km1*LK*LpA + AK*DF*ka3*L*LK*LpA + DF*K*ka3*L*LK*LpA - DF*ka3*Km1*LK^2*LpA - AKL*DF*kb3*L*LpAK + A*DF*ka3*L*LK*LpAK + DF^2*ka3*L*LK*LpA*LpAK + 2*AKL*kb3*Km1*LpAKL - AK*kb3*L*LpAKL - K*kb3*L*LpAKL + K*kcat1*L*LpAKL - A*ka3*Km1*LK*LpAKL + kb3*Km1*LK*LpAKL - kcat1*Km1*LK*LpAKL - DF*ka3*Km1*LK*LpA*LpAKL - DF*kb3*L*LpAK*LpAKL + kb3*Km1*LpAKL^2 - K*kcat7*L*LpAPLp + kcat7*Km1*LK*LpAPLp - K*kcat7*L*LpP + kcat7*Km1*LK*LpP )/(-(K*L) + Km1*LK)
    2. dL = (AKL^2*DF*K*ka2*Lp - AK*AKL*DF*ka2*LK*Lp + AKL*DF*K*ka3*LK*LpA - AK*DF*ka3*LK^2*LpA - AKL^2*DF*kb3*LpAK + A*AKL*DF*ka3*LK*LpAK - AKL*DF*kb3*LK*LpAK + AKL*DF*kcat1*LK*LpAK - APLp*DF*kcat7*LK*LpAK + A*DF*ka3*LK^2*LpAK + DF*kcat1*LK^2*LpAK - AKL*DF^2*ka2*LK*Lp*LpAK + AKL*DF^2*ka3*LK*LpA*LpAK - AKL*K*kb2*LpAKL + AK*AKL*kb3*LpAKL - AKL*K*kcat1*LpAKL + APLp*K*kcat7*LpAKL - A*AK*ka3*LK*LpAKL - A*K*ka3*LK*LpAKL + AK*kb2*LK*LpAKL + AK*kb3*LK*LpAKL - K*kcat1*LK*LpAKL + AKL*DF*K*ka2*Lp*LpAKL - AK*DF*ka3*LK*LpA*LpAKL - AKL*DF*kb3*LpAK*LpAKL + DF*kb2*LK*LpAK*LpAKL + DF*kcat1*LK*LpAK*LpAKL - K*kb2*LpAKL^2 + AK*kb3*LpAKL^2 - K*kcat1*LpAKL^2 - DF*kcat7*LK*LpAK*LpAPLp + K*kcat7*LpAKL*LpAPLp - DF*kcat7*LK*LpAK*LpP + K*kcat7*LpAKL*LpP)/(-(DF*LK*LpAK) + K*LpAKL)
    3. dL = (AKL^2*DF*ka2*Km1*Lp - AK*AKL*DF*ka2*L*Lp + AKL*DF*ka3*Km1*LK*LpA - AK*DF*ka3*L*LK*LpA - AKL*DF*kb3*L*LpAK + AKL*DF*kcat1*L*LpAK - APLp*DF*kcat7*L*LpAK + A*DF*ka3*L*LK*LpAK + DF*kcat1*L*LK*LpAK - AKL*DF^2*ka2*L*Lp*LpAK - AKL*kb2*Km1*LpAKL - AKL*kcat1*Km1*LpAKL + APLp*kcat7*Km1*LpAKL + AK*kb2*L*LpAKL + AK*kb3*L*LpAKL - A*ka3*Km1*LK*LpAKL - kcat1*Km1*LK*LpAKL + AKL*DF*ka2*Km1*Lp*LpAKL + DF*kb2*L*LpAK*LpAKL + DF*kcat1*L*LpAK*LpAKL - kb2*Km1*LpAKL^2 - kcat1*Km1*LpAKL^2 - DF*kcat7*L*LpAK*LpAPLp + kcat7*Km1*LpAKL*LpAPLp - DF*kcat7*L*LpAK*LpP + kcat7*Km1*LpAKL*LpP) / ((-(DF*L*LpAK) + Km1*LpAKL))
    #Km1*(-(AKL*DF*ka2*LK*Lp) - DF*ka3*LK^2*LpA + AKL*kb3*LpAKL - A*ka3*LK*LpAKL + kb2*LK*LpAKL + kb3*LK*LpAKL - DF*ka3*LK*LpA*LpAKL + kb3*LpAKL^2) = L*(-(AKL*DF*K*ka2*Lp) - DF*K*ka3*LK*LpA + AKL*DF*kb3*LpAK - A*DF*ka3*LK*LpAK - DF^2*ka3*LK*LpA*LpAK + K*kb2*LpAKL + K*kb3*LpAKL + DF*kb3*LpAK*LpAKL)
=#
function eliminationSystem1!(du, u, p, t)
     #Parameter and concetrations, concentrations, and differentials definitions
     Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
     L, K , P , A, Lp,  LpA, LK, LpP, LpAK, LpAP , LpAKL , LpAPLp , AK , AP , AKL , APLp = u
    
     #dL = du[1] will have to be found from these constraints...
     #dL, dK, dAK, dAKL, dLK, dLpAK, dLpAKL are all from Mathematica
    #dL FROM FIRST MATHEMATICA CONSTRAINT; SEE BELOW THE function
    du[1] = (AKL^2*kb3*Km1 - AK*AKL*kb3*L - AKL*K*kb3*L + AKL*K*kcat1*L - APLp*K*kcat7*L - A*AKL*ka3*Km1*LK + AKL*kb3*Km1*LK - AKL*kcat1*Km1*LK + APLp*kcat7*Km1*LK + A*AK*ka3*L*LK + A*K*ka3*L*LK + K*kcat1*L*LK - A*ka3*Km1*LK^2 - kcat1*Km1*LK^2 - AKL*DF*ka3*Km1*LK*LpA + AK*DF*ka3*L*LK*LpA + DF*K*ka3*L*LK*LpA - DF*ka3*Km1*LK^2*LpA - AKL*DF*kb3*L*LpAK + A*DF*ka3*L*LK*LpAK + DF^2*ka3*L*LK*LpA*LpAK + 2*AKL*kb3*Km1*LpAKL - AK*kb3*L*LpAKL - K*kb3*L*LpAKL + K*kcat1*L*LpAKL - A*ka3*Km1*LK*LpAKL + kb3*Km1*LK*LpAKL - kcat1*Km1*LK*LpAKL - DF*ka3*Km1*LK*LpA*LpAKL - DF*kb3*L*LpAK*LpAKL + kb3*Km1*LpAKL^2 - K*kcat7*L*LpAPLp + kcat7*Km1*LK*LpAPLp - K*kcat7*L*LpP + kcat7*Km1*LK*LpP ) / (-(K*L) + Km1*LK)
    #dK 
    du[2] = -(A*K*ka3) + AK*kb3 + AKL*kb3 - A*ka3*LK - K*ka3*LpA - DF*ka3*LK*LpA + kb3*LpAK + kb3*LpAKL 
    #dP 
    du[3] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P 
    #dA
    du[4] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P 
    #dLp
    du[5] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*DF*AKL*Lp - ka2*DF*APLp*Lp - ka7*DF*Lp*LpAP 
    #dLpA
    du[6] = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*DF*LK*LpA - ka4*DF*LpA*LpP 
    #dLK (MM assumption)
    du[7] = 0
    #dLpP
    du[8] = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*DF*LpA*LpP 
    #dLpAK
    du[9] = AK*ka2*Lp + AKL*DF*ka2*Lp + K*ka3*LpA + DF*ka3*LK*LpA - kb2*LpAK - kb3*LpAK - kb2*LpAKL - kb3*LpAKL 
    #dLpAP 
    du[10] = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*DF*Lp*LpAP 
    #dLpAKL (MM Assumption)
    du[11] = 0
    #dLpAPLp
    du[12] = ka2*DF*APLp*Lp + ka7*DF*Lp*LpAP + ka4*DF*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp 
    #dAK
    du[13] = du[1] + A*K*ka3 - AK*kb3 - AKL*kb3 + AKL*kcat1 - APLp*kcat7 + A*ka3*LK + kcat1*LK - AK*ka2*Lp - AKL*DF*ka2*Lp + kb2*LpAK + kb2*LpAKL + kcat1*LpAKL - kcat7*LpAPLp - kcat7*LpP 
    #dAP
    du[14] = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp 
    #dAKL
    du[15] = -du[1] - AKL*kcat1 + APLp*kcat7 - kcat1*LK - kcat1*LpAKL + kcat7*LpAPLp + kcat7*LpP 
    #dAPLp
    du[16] = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*DF*APLp*Lp
    nothing
end

function eliminationSystem2!(du, u, p, t)
     #Parameter and concetrations, concentrations, and differentials definitions
     Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
     L, K , P , A, Lp,  LpA, LK, LpP, LpAK, LpAP , LpAKL , LpAPLp , AK , AP , AKL , APLp = u
    
     #dL = du[1] will have to be found from these constraints...
     #dL, dK, dAK, dAKL, dLK, dLpAK, dLpAKL are all from Mathematica
    #dL FROM FIRST MATHEMATICA CONSTRAINT; SEE BELOW THE function
    du[1] = (AKL^2*DF*K*ka2*Lp - AK*AKL*DF*ka2*LK*Lp + AKL*DF*K*ka3*LK*LpA - AK*DF*ka3*LK^2*LpA - AKL^2*DF*kb3*LpAK + A*AKL*DF*ka3*LK*LpAK - AKL*DF*kb3*LK*LpAK + AKL*DF*kcat1*LK*LpAK - APLp*DF*kcat7*LK*LpAK + A*DF*ka3*LK^2*LpAK + DF*kcat1*LK^2*LpAK - AKL*DF^2*ka2*LK*Lp*LpAK + AKL*DF^2*ka3*LK*LpA*LpAK - AKL*K*kb2*LpAKL + AK*AKL*kb3*LpAKL - AKL*K*kcat1*LpAKL + APLp*K*kcat7*LpAKL - A*AK*ka3*LK*LpAKL - A*K*ka3*LK*LpAKL + AK*kb2*LK*LpAKL + AK*kb3*LK*LpAKL - K*kcat1*LK*LpAKL + AKL*DF*K*ka2*Lp*LpAKL - AK*DF*ka3*LK*LpA*LpAKL - AKL*DF*kb3*LpAK*LpAKL + DF*kb2*LK*LpAK*LpAKL + DF*kcat1*LK*LpAK*LpAKL - K*kb2*LpAKL^2 + AK*kb3*LpAKL^2 - K*kcat1*LpAKL^2 - DF*kcat7*LK*LpAK*LpAPLp + K*kcat7*LpAKL*LpAPLp - DF*kcat7*LK*LpAK*LpP + K*kcat7*LpAKL*LpP)/(-(DF*LK*LpAK) + K*LpAKL)
    #dK 
    du[2] = -(A*K*ka3) + AK*kb3 + AKL*kb3 - A*ka3*LK - K*ka3*LpA - DF*ka3*LK*LpA + kb3*LpAK + kb3*LpAKL 
    #dP 
    du[3] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P 
    #dA
    du[4] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P 
    #dLp
    du[5] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*DF*AKL*Lp - ka2*DF*APLp*Lp - ka7*DF*Lp*LpAP 
    #dLpA
    du[6] = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*DF*LK*LpA - ka4*DF*LpA*LpP 
    #dLK (MM assumption)
    du[7] = 0
    #dLpP
    du[8] = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*DF*LpA*LpP 
    #dLpAK
    du[9] = AK*ka2*Lp + AKL*DF*ka2*Lp + K*ka3*LpA + DF*ka3*LK*LpA - kb2*LpAK - kb3*LpAK - kb2*LpAKL - kb3*LpAKL 
    #dLpAP 
    du[10] = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*DF*Lp*LpAP 
    #dLpAKL (MM Assumption)
    du[11] = 0
    #dLpAPLp
    du[12] = ka2*DF*APLp*Lp + ka7*DF*Lp*LpAP + ka4*DF*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp 
    #dAK
    du[13] = du[1] + A*K*ka3 - AK*kb3 - AKL*kb3 + AKL*kcat1 - APLp*kcat7 + A*ka3*LK + kcat1*LK - AK*ka2*Lp - AKL*DF*ka2*Lp + kb2*LpAK + kb2*LpAKL + kcat1*LpAKL - kcat7*LpAPLp - kcat7*LpP 
    #dAP
    du[14] = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp 
    #dAKL
    du[15] = -du[1] - AKL*kcat1 + APLp*kcat7 - kcat1*LK - kcat1*LpAKL + kcat7*LpAPLp + kcat7*LpP 
    #dAPLp
    du[16] = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*DF*APLp*Lp
    nothing
end

function eliminationSystem3!(du, u, p, t)
    #Parameter and concetrations, concentrations, and differentials definitions
    Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
    L, K , P , A, Lp,  LpA, LK, LpP, LpAK, LpAP , LpAKL , LpAPLp , AK , AP , AKL , APLp = u
   
    #dL = du[1] will have to be found from these constraints...
    #dL, dK, dAK, dAKL, dLK, dLpAK, dLpAKL are all from Mathematica
    du[1] = (AKL^2*DF*ka2*Km1*Lp - AK*AKL*DF*ka2*L*Lp + AKL*DF*ka3*Km1*LK*LpA - AK*DF*ka3*L*LK*LpA - AKL*DF*kb3*L*LpAK + AKL*DF*kcat1*L*LpAK - APLp*DF*kcat7*L*LpAK + A*DF*ka3*L*LK*LpAK + DF*kcat1*L*LK*LpAK - AKL*DF^2*ka2*L*Lp*LpAK - AKL*kb2*Km1*LpAKL - AKL*kcat1*Km1*LpAKL + APLp*kcat7*Km1*LpAKL + AK*kb2*L*LpAKL + AK*kb3*L*LpAKL - A*ka3*Km1*LK*LpAKL - kcat1*Km1*LK*LpAKL + AKL*DF*ka2*Km1*Lp*LpAKL + DF*kb2*L*LpAK*LpAKL + DF*kcat1*L*LpAK*LpAKL - kb2*Km1*LpAKL^2 - kcat1*Km1*LpAKL^2 - DF*kcat7*L*LpAK*LpAPLp + kcat7*Km1*LpAKL*LpAPLp - DF*kcat7*L*LpAK*LpP + kcat7*Km1*LpAKL*LpP) / ((-(DF*L*LpAK) + Km1*LpAKL))
    #dK 
    du[2] = -(A*K*ka3) + AK*kb3 + AKL*kb3 - A*ka3*LK - K*ka3*LpA - DF*ka3*LK*LpA + kb3*LpAK + kb3*LpAKL 
    #dP 
    du[3] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P 
    #dA
    du[4] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P 
    #dLp
    du[5] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*DF*AKL*Lp - ka2*DF*APLp*Lp - ka7*DF*Lp*LpAP 
    #dLpA
    du[6] = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*DF*LK*LpA - ka4*DF*LpA*LpP 
    #dLK (MM assumption)
    du[7] = 0
    #dLpP
    du[8] = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*DF*LpA*LpP 
    #dLpAK
    du[9] = AK*ka2*Lp + AKL*DF*ka2*Lp + K*ka3*LpA + DF*ka3*LK*LpA - kb2*LpAK - kb3*LpAK - kb2*LpAKL - kb3*LpAKL 
    #dLpAP 
    du[10] = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*DF*Lp*LpAP 
    #dLpAKL (MM Assumption)
    du[11] = 0
    #dLpAPLp
    du[12] = ka2*DF*APLp*Lp + ka7*DF*Lp*LpAP + ka4*DF*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp 
    #dAK
    du[13] = du[1] + A*K*ka3 - AK*kb3 - AKL*kb3 + AKL*kcat1 - APLp*kcat7 + A*ka3*LK + kcat1*LK - AK*ka2*Lp - AKL*DF*ka2*Lp + kb2*LpAK + kb2*LpAKL + kcat1*LpAKL - kcat7*LpAPLp - kcat7*LpP 
    #dAP
    du[14] = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp 
    #dAKL
    du[15] = -du[1] - AKL*kcat1 + APLp*kcat7 - kcat1*LK - kcat1*LpAKL + kcat7*LpAPLp + kcat7*LpP 
    #dAPLp
    du[16] = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*DF*APLp*Lp
    nothing
end

include("../UTILITIES/ODESystems.jl")

pMM = zeros(12)
Km1 = sum(p[2:3])/p[1]
pMM[1] = Km1
pMM[2:12] = p[3:13]
usym = [:L => 10^1.6,  :K => 10^0.8, :P => 10^0.6, :A => 10^1.2,:Lp => 0.1, :LpA => 0.1, :LK => 0.1, 
        :LpP => 0.1, :LpAK => 0.1, :LpAP => 0.1, :LpAKL => 0.1, :LpAPLp => 0.1 , :AK => 0.1, :AP => 0.1, 
        :AKL => 0.1, :APLp => 0.1]
myu0=[i[2] for i in usym]
funcArr = [eliminationSystem1!, eliminationSystem2!, eliminationSystem3!]
ODEArr = [ODEProblem(i, myu0, tspan, pMM)  for i in funcArr]
solArr = [solve(i, Rosenbrock23(), save_idxs = 1, save_at = 0.1) for i in ODEArr]
solPlots = plot(solve(ODEProblem(fullmodel_ode!, myu0, tspan, p), Rosenbrock23(), save_idxs = 1, saveat = 0.1), label = "full")
for i in enumerate(solArr) plot!(i[2], label = "$(i[1])"); print(i[1]) end
solPlots