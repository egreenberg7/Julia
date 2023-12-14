using DifferentialEquations
using Plots
using CSV
using DataFrames

function twoVarSystemKdToZero(du, u, p, t)
    L, LpA = u
    n1DF, n7DF, ka2, kb2, ka3, ka4, Ltot, Ktot, Ptot, Atot = p
    du[1] = -n1DF*Ktot*L - n7DF*Ptot*(Ltot-(LpA+Ktot+Ptot)-L)
    du[2] = -LpA*(Ktot*ka3+Ptot*ka4)+ka2*(Ltot-L-(LpA+Ktot+Ptot))*(Atot-(LpA+Ktot+Ptot))-kb2*LpA
    if L < 0 && LpA < 0
        throw("Negative concentrations")
    end
end

WorkingParams = CSV.read("ReducedModels/VeryReducedModel/ViablePoints_12varModelCleaned.csv", DataFrame)

function calculateSelectivity(ka, kb, kcat) 
    return ka * kcat / (kb + kcat)
end

function getP(row)
    r = WorkingParams[row, :]
    n1DF = calculateSelectivity(r.ka1 * r.VA,r.kb1,r.kcat1)
    n7DF = calculateSelectivity(r.ka7 * r.VA, r.kb7,r.kcat7)
    return [n1DF,n7DF,r.ka2,r.kb2,r.ka3 * 10,r.ka4 * 10,r.Lp,r.K,r.P,r.A]
end

function getU(row)
    c = 1
    r = WorkingParams[row,:]
    return [r.Lp *c, r.A]
end

function makeProb(row; tspan = (0,500))
    println(getP(row))
    return ODEProblem(twoVarSystemKdToZero, getU(row), tspan, getP(row))
end

for i in 1:50
    prob = makeProb(i; tspan=(0,1))
    solve(prob)
    println("Made it!")
    try
        display(plot(prob))
    catch y
        if isa(y, TypeError)
            println{"Instability detected"}
        end
    end
end

for i in 1:100 
    #ap = 0.01.*[1,100, 3,3,40,30,8,1,1,5]
    #au = 0.01.*[5,2]
    ap = [rand() *100, rand()*100,rand()*10,rand()*10,rand()*100,rand()*100,rand()*20,rand()*5,rand()*5,rand()*8]
    au = [rand() * ap[7], rand()]
    try
        sol=solve(ODEProblem(twoVarSystemKdToZero,au,(0,10),ap), maxtime=30) 
    if sol.retcode == SciMLBase.ReturnCode.Success
        display(plot(sol))
        println("p: " * string(ap))
        println("u: " * string(au))
    end
    catch y
    end
end

