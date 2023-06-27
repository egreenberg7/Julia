using DifferentialEquations
#Full reaction model with kb1 eliminated using dLK = 0 and dLpAKL = 0
function eliminationSystem()
     #Parameter and concetrations, concentrations, and differentials definitions
     p = ka1, Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF
     #Equations from Mathematica
     dAK = A  *  K  *  ka3 - AK  *  kb3 + AKL  *  ka1  *  Km1 - AK  *  ka1  *  L - AK  *  ka2  *  Lp + kb2  *  LpAK  
     dAKL =  - AKL * kb3 - AKL * ka1 * Km1 + AK * ka1 * L + A * ka3 * LK - AKL * DF * ka2 * Lp + kb2 * LpAKL
     dK =  - A * K * ka3 + AK * kb3 + AKL * kb3 - A * ka3 * LK - K * ka3 * LpA - DF * ka3 * LK * LpA + kb3 * LpAK + kb3 * LpAKL
     dL = AKL * kb3 - AKL * kcat1 + APLp * kcat7 + AKL * ka1 * Km1 - AK * ka1 * L - A * ka3 * LK - kcat1 * LK + AKL * DF * ka2 * Lp - kb2 * LpAKL - kcat1 * LpAKL + kcat7 * LpAPLp + kcat7 * LpP
     dLK = 0
     dLpAK = AK * ka2 * Lp + AKL * DF * ka2 * Lp + K * ka3 * LpA + DF * ka3 * LK * LpA - kb2 * LpAK - kb3 * LpAK - kb2 * LpAKL - kb3 * LpAKL
     dLpAKL = 0
     Km1 * LK = (AKL * kb3)/ka1 + K * L - (A * ka3 * LK)/ka1 - (DF * ka3 * LK * LpA)/ka1 + (kb3 * LpAKL)/ka1
     Km1 * LpAKL = (AKL * DF * ka2 * Lp)/ka1 + (DF * ka3 * LK * LpA)/ka1 + DF * L * LpAK - (kb2 * LpAKL)/ka1 - (kb3 * LpAKL)/ka1
     kb2 * LK * LpAKL = AKL * DF * ka2 * LK * Lp + DF * ka3 * LK^2 * LpA + DF * ka1 * L * LK * LpAK - AKL * kb3 * LpAKL - K * ka1 * L * LpAKL + A * ka3 * LK * LpAKL - kb3 * LK * LpAKL + DF * ka3 * LK * LpA * LpAKL - kb3 * LpAKL^2
     
     #Original equations from reaction system
     dLp = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*DF*AKL*Lp - ka2*DF*APLp*Lp - ka7*DF*Lp*LpAP 
     dP = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P 
     dA = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P 
     dLpA = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*DF*LK*LpA - ka4*DF*LpA*LpP 
     dLpP = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*DF*LpA*LpP 
     dLpAP = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*DF*Lp*LpAP 
     dLpAPLp = ka2*DF*APLp*Lp + ka7*DF*Lp*LpAP + ka4*DF*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp 
     dAP = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp 
     dAPLp = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*DF*APLp*Lp
end