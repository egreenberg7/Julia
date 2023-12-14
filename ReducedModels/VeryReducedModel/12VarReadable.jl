using DifferentialEquations
using DataFrames
using CSV
using Plots
include("../../ExperimentalFullModelWork/OutputHandling.jl")

WorkingParams = CSV.read("ReducedModels/VeryReducedModel/ViablePoints_12varModelCleaned.csv", DataFrame)
function getKValues(row)
    p = WorkingParams[row, :]
    return [p.ka1,p.kb1,p.kcat1,p.ka2,p.ka3,p.kb3,p.kb4,p.ka7,p.kb7,p.kcat7, p.VA, p.Lp,p.K,p.P,p.A]
end

getKValues(1)
oscdata = DataFrame(CSV.File("/Users/ezragreenberg/JLab/Julia/ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))
function getSelectivities(row)
    p = WorkingParams[row, :]
    n1 = (p.ka1 * p.kcat1) / (p.kb1 + p.kcat1)
    n1DF = (p.ka1 * p.VA * p.kcat1) / (p.kb1 + p.VA * p.kcat1)
    n7 = (p.ka7 * p.kcat7) / (p.kb7 + p.kcat7)
    n7DF = (p.ka7 * p.VA* p.kcat7) / (p.kb7 + p.VA * p.kcat7)
    return [n1,n1DF,n7,n7DF]
end

function getp(row)
    p = WorkingParams[row, :]
    return vcat(getSelectivities(row), p.ka2,p.kb2,p.ka3,p.kb3,p.ka4,p.kb4,p.Lp + (p.Lp + p.A)/1.5,p.K,p.P,p.A + (p.Lp+p.A)/1.5)
end

function getSelectivitiesExperimental(row)
    p = getP(oscdata, row)
    ka1 = p[1]
    kb1 = p[2] * 10
    kcat1 = p[3]
    VA=p[end]
    ka7 = p[10]
    kb7 = p[11]
    kcat7 = p[12]
    n1 = (ka1 * kcat1) / (kb1 + kcat1)
    n1DF = (ka1 * VA * kcat1) / (kb1 + VA * kcat1)
    n7 = (ka7 * kcat7) / (kb7 + kcat7)
    n7DF = (ka7 * VA* kcat7) / (kb7 + VA * kcat7)
    return [n1,n1DF,n7,n7DF]
end

function getpExperimental(row)
    p = getP(oscdata, row)
    uDF = oscdata[row, :]
    return vcat(getSelectivitiesExperimental(row), ka2exp ,kb2exp,ka3exp,kb3exp, p[8] ,p[9], uDF.L/2 +uDF.A*0.2,uDF.K,uDF.P,uDF.A *1.2)
end

function getu2(row) 
    p = oscdata[row, :]
    return [p.L/2.0, p.A * 0.2]
end

function getU(row)
    p = WorkingParams[row,:]
    return [p.Lp, 0]
end

myP = getp(2)
myU = getU(2)


function twoVarSystem(du, u, p, t)
    L, LpA = u
    n1, n1DF, n7, n7DF, ka2, kb2, ka3, kb3, ka4, kb4, Ltot, Ktot, Ptot, Atot = p
    du[1] = L*((-n1*Ktot)/(1.0+LpA/(kb3/ka3)) - n1DF*LpA*Ktot/(kb3/ka3)+LpA) + 
            (Ltot - L - LpA * (1.0 + (Ktot/(kb3/ka3+LpA) + Ptot/(kb4/ka4 + LpA))) * Ptot * 
            ((n7DF*LpA)/(kb4/ka4+LpA) + n7/(1.0+LpA/kb4/ka4)))
    du[2] = (LpA * Ktot * (kb3 - ka3*LpA)/(kb3/ka3 + LpA)) + (LpA * Ptot * (kb4 - ka4*LpA)/(kb4/ka4 + LpA)) +
            ka2*(Ltot-L-LpA*(1+(Ktot/(kb3/ka3+LpA))+(Ptot/(kb4/ka4+LpA))))*
            (Atot-LpA*(1+(Ktot/(kb3/ka3+LpA))+(Ptot/(kb4/ka4+LpA)))) - kb2 * LpA
end

myprob = ODEProblem(twoVarSystem,myU,(0,500),myP)
mysol = solve(myprob)
myplot = plot(mysol)

function generatePlot(row)
    myP=getp(row)
    myU=getU(row)
    myprob = ODEProblem(twoVarSystem,myU,(0,3),myP)
    mysol = solve(myprob)
    plot(mysol)
end

generatePlot(25)

function generatePlot2(row)
    myP=getpExperimental(row)
    myU=getu2(row)
    myprob = ODEProblem(twoVarSystem, myU, (0,500),myP)
    mysol = solve(myprob,Rodas4P())
    plot(mysol)
end

for i in 400:500
    display(generatePlot2(i))
end

function traceAtFixedPoint(Ltot, Ktot, Ptot, Atot, ka2, kb2, ka3, kb3, ka4, kb4, n1DF, nDF7) 
    function Sqrt(x) return x^0.5 end
    return -(Atot*ka2) - kb2 - Ktot*(ka3 + n1DF) - (ka4 + nDF7)*Ptot + ka2*(2*Ktot - Ltot + 2*Ptot) + 
    (2*Atot*ka2*Ktot^2*n1DF^2 + 2*kb2*Ktot^2*n1DF^2 - 4*ka2*Ktot^3*n1DF^2 + 2*ka3*Ktot^3*n1DF^2 + 
      2*ka2*Ktot^2*Ltot*n1DF^2 - 4*ka2*Ktot^2*n1DF^2*Ptot + 2*ka4*Ktot^2*n1DF^2*Ptot + 
      Atot*ka2*Ktot*n1DF*nDF7*Ptot + 3*kb2*Ktot*n1DF*nDF7*Ptot - 4*ka2*Ktot^2*n1DF*nDF7*Ptot + 
      3*ka3*Ktot^2*n1DF*nDF7*Ptot + 3*ka2*Ktot*Ltot*n1DF*nDF7*Ptot - 4*ka2*Ktot*n1DF*nDF7*Ptot^2 + 
      3*ka4*Ktot*n1DF*nDF7*Ptot^2 + kb2*nDF7^2*Ptot^2 + ka3*Ktot*nDF7^2*Ptot^2 + ka4*nDF7^2*Ptot^3 - 
      (2*Ktot*n1DF + nDF7*Ptot)*Sqrt(Atot^2*ka2^2*Ktot^2*n1DF^2 + 
         Ktot^2*(kb2^2 + ka3^2*Ktot^2 + ka2^2*Ltot^2 + 2*ka2*Ktot*(-2*(kb3 + ka3*Ktot) + ka3*Ltot) + 
           2*kb2*(ka3*Ktot + ka2*(-2*Ktot + Ltot)))*n1DF^2 + 
         2*Ktot*n1DF*(ka2*Ktot*(-2*(kb2 + kb4 + (ka3 + ka4)*Ktot) + ka4*Ltot)*n1DF + 
           ka2*(-2*Ktot*(kb2 + kb3 + ka3*Ktot) + (kb2 + ka3*Ktot)*Ltot)*nDF7 + 
           (kb2 + ka3*Ktot)*(ka4*Ktot*n1DF + (kb2 + ka3*Ktot)*nDF7))*Ptot + 
         (ka4*(-4*ka2 + ka4)*Ktot^2*n1DF^2 + 2*Ktot*(2*ka4*(kb2 + ka3*Ktot) - 
             2*ka2*(kb2 + kb4 + (ka3 + ka4)*Ktot) + ka2*ka4*Ltot)*n1DF*nDF7 + (kb2 + ka3*Ktot)^2*nDF7^2)*
          Ptot^2 + 2*ka4*nDF7*(-2*ka2*Ktot*n1DF + ka4*Ktot*n1DF + (kb2 + ka3*Ktot)*nDF7)*Ptot^3 + 
         ka4^2*nDF7^2*Ptot^4 + 2*Atot*ka2*Ktot*n1DF*(Ktot*(kb2 + ka3*Ktot - ka2*Ltot)*n1DF + 
           ka4*Ktot*n1DF*Ptot + (kb2 + ka3*Ktot)*nDF7*Ptot + ka4*nDF7*Ptot^2)))/
     (2*Ktot*n1DF*(Ktot*n1DF + nDF7*Ptot))
end

nRange = 1:100
count = 0
for Ltot in Lrange
    for Ktot in Krange
        for Ptot in Prange
            for Atot in Arange
                for ka2 in kaRange
                    for kb2 in kbRange
                        #for ka3 in kaRange
                            for kb3 in kbRange
                                #for ka4 in kaRange
                                    for kb4 in kbRange
                                        for n1DF in nRange
                                            for nDF7 in nRange
                                                if traceAtFixedPoint(Ltot, Ktot, Ptot, Atot, ka2, kb2, 100 * kb3, kb3, 100*kb4, kb4, n1DF, nDF7) < 0
                                                    count=count+1
                                                    #println(Ltot, Ktot, Ptot, Atot, ka2, kb2, ka3, kb3, ka4, kb4, n1DF, nDF7)
                                                    for i in [Ltot, Ktot, Ptot, Atot, ka2, kb2, 100*kb3, kb3, 100*kb4, kb4, n1DF, nDF7]
                                                    if (count == 1000)
                                                        println(i)
                                                        throw(Exception)
                                                    end
                                                end

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
myvalues = [0.010000000000000002,0.001,0.001,0.010000000000000002,0.001,0.001,0.1,0.001,0.1,0.001,1.0,1.0]
myprob = ODEProblem(twoVarSystem, [0.4,0.001], (0,500),vcat([myvalues[end-1]/2,myvalues[end-1]/2,myvalues[end],myvalues[end]], myvalues[5:end], myvalues[1:4]))
plot(solve(myprob))
