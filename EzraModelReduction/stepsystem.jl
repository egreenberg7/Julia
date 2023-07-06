using DifferentialEquations
using Plots

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

function mmLKmodel_ode!(du, u, p, t)
    L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
    Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y = p #parameters

    ka1 = (-(AKL*kb3) + A*ka3*LK + kcat1*LK + y*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + y*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK))/(K*L)
    kb1 = (-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + y*ka3*Km1*LK*LpA - kb3*Km1*LpAKL)/(K*L - Km1*LK)

    du[1] = kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK #* L
    du[2] = kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA #* K
    du[3] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P #* P
    du[4] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P #* A
    du[5] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP #* Lp
    du[6] = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*y*LK*LpA - ka4*y*LpA*LpP #* LpA
    du[7] = 0 #* LK
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

function mmLKsubmodel_ode!(du,u, p, t)
    L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
    Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p #parameters


      du[1] = AKL*kb3 + APLp*kcat7 - A*ka3*LK - kcat1*LK - DF*ka3*LK*LpA + kb3*LpAKL + (AKL*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK) + (LpAKL*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK) - (AK*(-(AKL*kb3) + A*ka3*LK + kcat1*LK + DF*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK)))/K - (DF*LpAK*(-(AKL*kb3) + A*ka3*LK + kcat1*LK + DF*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK)))/K + kcat7*LpAPLp + kcat7*LpP #*L
      du[5] = APLp*kb7 + AKL*kcat1 + kcat1*LK - A*ka2*Lp - AK*ka2*Lp - AP*ka2*Lp - AKL*DF*ka2*Lp - APLp*DF*ka2*Lp - AP*ka7*Lp + kb2*LpA + kb2*LpAK + kb2*LpAKL + kcat1*LpAKL + kb2*LpAP - DF*ka7*Lp*LpAP + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka7*Lp*P #*Lp 
      du[2] = -(A*K*ka3) + AK*kb3 + AKL*kb3 - A*ka3*LK - K*ka3*LpA - DF*ka3*LK*LpA + kb3*LpAK + kb3*LpAKL #*K
      du[3] = AP*kb4 + kb4*LpAP + kb7*LpP + kcat7*LpP - A*ka4*P - ka7*Lp*P - ka4*LpA*P #*P
      du[4] = -(A*K*ka3) + AK*kb3 + AKL*kb3 + AP*kb4 + APLp*kb4 - A*ka3*LK - A*ka2*Lp + kb2*LpA - A*ka4*LpP - A*ka4*P #*A 
      du[6] = A*ka2*Lp - K*ka3*LpA - kb2*LpA - DF*ka3*LK*LpA + kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp - DF*ka4*LpA*LpP - ka4*LpA*P #*LpA 
      du[7] = 0 #*LK 
      du[8] = APLp*kb4 + kb4*LpAPLp - A*ka4*LpP - kb7*LpP - kcat7*LpP - DF*ka4*LpA*LpP + ka7*Lp*P #*LpP 
      du[9] = AK*ka2*Lp + K*ka3*LpA - kb2*LpAK - kb3*LpAK + kcat1*LpAKL + (LpAKL*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK) - (DF*LpAK*(-(AKL*kb3) + A*ka3*LK + kcat1*LK + DF*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK)))/K #*LpAK 
      du[10] = AP*ka2*Lp - kb2*LpAP - kb4*LpAP - DF*ka7*Lp*LpAP + kb7*LpAPLp + kcat7*LpAPLp + ka4*LpA*P #*LpAP
      du[11] = AKL*DF*ka2*Lp + DF*ka3*LK*LpA - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL - (LpAKL*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK) + (DF*LpAK*(-(AKL*kb3) + A*ka3*LK + kcat1*LK + DF*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK)))/K #*LpAKL 
      du[12] = APLp*DF*ka2*Lp + DF*ka7*Lp*LpAP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp + DF*ka4*LpA*LpP #*LpAPLp
      du[13] = A*K*ka3 - AK*kb3 + AKL*kcat1 - AK*ka2*Lp + kb2*LpAK + (AKL*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK) - (AK*(-(AKL*kb3) + A*ka3*LK + kcat1*LK + DF*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK)))/K #*AK 
      du[14] = -(AP*kb4) + APLp*kb7 + APLp*kcat7 - AP*ka2*Lp - AP*ka7*Lp + kb2*LpAP + A*ka4*P #*AP 
      du[15] = -(AKL*kb3) - AKL*kcat1 + A*ka3*LK - AKL*DF*ka2*Lp + kb2*LpAKL - (AKL*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK) + (AK*(-(AKL*kb3) + A*ka3*LK + kcat1*LK + DF*ka3*LK*LpA - kb3*LpAKL + (LK*(-(AKL*kb3*Km1) - K*kcat1*L + A*ka3*Km1*LK + kcat1*Km1*LK + DF*ka3*Km1*LK*LpA - kb3*Km1*LpAKL))/(K*L - Km1*LK)))/K #*AKL 
      du[16] = -(APLp*kb4) - APLp*kb7 - APLp*kcat7 - APLp*DF*ka2*Lp + AP*ka7*Lp + kb2*LpAPLp + A*ka4*LpP #*APLp
end

function mmwithMassModel_ode!(du,u, p, t)
    L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
    Km1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p #parameters
    #Ltot
    Ltot = L + Lp + LpA + LpP + LpAK + LpAP + 2 * LpAKL +  2 * LpAPLp + AKL + APLp + LK
    #dL
    du[1] = AKL*kb3 + APLp*kcat7 + kb3*LpAKL + kcat7*LpAPLp + kcat7*LpP - A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + (AKL*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA* (-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/ (K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)) + (LpAKL*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA* (-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/ (K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)) - (AK*(-(AKL*kb3) - kb3*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + ((-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2* LpAPLp - LpP + Ltot)))/(K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot))))/K - (DF*LpAK*(-(AKL*kb3) - kb3*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + ((-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2* LpAPLp - LpP + Ltot)))/(K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot))))/K
    #dK
    du[2] = -(A*K*ka3) + AK*kb3 + AKL*kb3 - K*ka3*LpA + kb3*LpAK + kb3*LpAKL - A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)
    #dP
    du[3] = AP*kb4 + kb4*LpAP + kb7*LpP + kcat7*LpP - A*ka4*P - ka7*Lp*P - ka4*LpA*P
    #dA
    du[4] = -(A*K*ka3) + AK*kb3 + AKL*kb3 + AP*kb4 + APLp*kb4 - A*ka2*Lp + kb2*LpA - A*ka4*LpP - A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - A*ka4*P 
    #dLp
    du[5] = APLp*kb7 + AKL*kcat1 - A*ka2*Lp - AK*ka2*Lp - AP*ka2*Lp - AKL*DF*ka2*Lp - APLp*DF*ka2*Lp - AP*ka7*Lp + kb2*LpA + kb2*LpAK + kb2*LpAKL + kcat1*LpAKL + kb2*LpAP - DF*ka7*Lp*LpAP + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - ka7*Lp*P 
    #dLpA
    du[6] = A*ka2*Lp - K*ka3*LpA - kb2*LpA + kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp - DF*ka4*LpA*LpP - DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - ka4*LpA*P
    #dLK 
    du[7] = 0
    #dLpP
    du[8] = APLp*kb4 + kb4*LpAPLp - A*ka4*LpP - kb7*LpP - kcat7*LpP - DF*ka4*LpA*LpP + ka7*Lp*P
    #dLpAK
    du[9] = AK*ka2*Lp + K*ka3*LpA - kb2*LpAK - kb3*LpAK + kcat1*LpAKL + (LpAKL*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA* (-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/ (K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)) - (DF*LpAK*(-(AKL*kb3) - kb3*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + ((-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2* LpAPLp - LpP + Ltot)))/(K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot))))/K 
    #dLpAP
    du[10] = AP*ka2*Lp - kb2*LpAP - kb4*LpAP - DF*ka7*Lp*LpAP + kb7*LpAPLp + kcat7*LpAPLp + ka4*LpA*P
    #dLpAKL
    du[11] = AKL*DF*ka2*Lp - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - (LpAKL*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA* (-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/ (K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)) + (DF*LpAK*(-(AKL*kb3) - kb3*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + ((-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2* LpAPLp - LpP + Ltot)))/(K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot))))/K
    #dLpAPLp
    du[12] = APLp*DF*ka2*Lp + DF*ka7*Lp*LpAP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp + DF*ka4*LpA*LpP
    #dAK
    du[13] = A*K*ka3 - AK*kb3 + AKL*kcat1 - AK*ka2*Lp + kb2*LpAK + (AKL*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA* (-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/ (K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)) - (AK*(-(AKL*kb3) - kb3*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + ((-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2* LpAPLp - LpP + Ltot)))/(K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot))))/K
    #dAP
    du[14] = -(AP*kb4) + APLp*kb7 + APLp*kcat7 - AP*ka2*Lp - AP*ka7*Lp + kb2*LpAP + A*ka4*P
    #dAKL
    du[15] = -(AKL*kb3) - AKL*kcat1 - AKL*DF*ka2*Lp + kb2*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) - (AKL*(-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/(K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)) + (AK*(-(AKL*kb3) - kb3*LpAKL + A*ka3*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*LpA*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + ((-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)* (-(AKL*kb3*Km1) - K*kcat1*L - kb3*Km1*LpAKL + A*ka3*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + kcat1*Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot) + DF*ka3*Km1*LpA* (-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot)))/ (K*L - Km1*(-AKL - APLp - L - Lp - LpA - LpAK - 2*LpAKL - LpAP - 2*LpAPLp - LpP + Ltot))))/K
    #dAPLp
    du[16] = -(APLp*kb4) - APLp*kb7 - APLp*kcat7 - APLp*DF*ka2*Lp + AP*ka7*Lp + kb2*LpAPLp + A*ka4*LpP
end

begin
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => 0.055, :kb1 => 19.8, :kcat1 => 241, :ka2 => 1, :kb2 => 0.95,
        :ka3 => 41, :kb3 => 193, :ka4 => 0.19, :kb4 => 0.13, :ka7 => 0.62, 
        :kb7 => 3.39, :kcat7 => 4.6, :y => 750]
p = [x[2] for x in psym]
    
#initial condition list
usym = [:L => 1.5, :K => 0.2, :P => 0.3, :A => 0.6, :Lp => 3.0, :LpA => 0.0, :LK => 0, 
        :LpP => 0., :LpAK => 0.1, :LpAP => 0.01, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0., 
        :AP => 0., :AKL => 0, :APLp => 0.]
u0 = [x[2] for x in usym]
pMM = zeros(12)
Km1 = sum(p[2:3])/p[1]
pMM[1] = Km1
pMM[2:12] = p[3:13]
tspan = (0, 100)
#
usym = [:L => 10^1.6,  :K => 10^0.8, :P => 10^0.6, :A => 10^1.2,:Lp => 0, :LpA => 0, :LK => 0, 
        :LpP => 0, :LpAK => 0, :LpAP => 0, :LpAKL => 0, :LpAPLp => 0 , :AK => 0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0=[i[2] for i in usym]

probfull = ODEProblem(fullmodel_ode!, u0, tspan, p)
probMM = ODEProblem(mmwithMassModel_ode!,u0,tspan,pMM)
solfull = solve(probfull, Rosenbrock23(), save_idxs = [1,5])
solMM = solve(probMM, Rosenbrock23(), save_idxs = [1,5])
full = Plots.plot(solfull)
end
Plots.plot!(solMM)