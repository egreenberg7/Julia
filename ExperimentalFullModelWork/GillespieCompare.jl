include("OutputHandling.jl")
include("../UTILITIES/GillespieConverter.jl")

oscdata = DataFrame(CSV.File("/Users/ezragreenberg/Julia/ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))

for i in 1:size(oscdata)[1]
    sol = entryToSol(oscdata, i, tspan=600.0)
    exSolPlt = Plots.plot(sol, title="Representative Solution",
        label="Numerical Solution", xlabel = "Time (s)", ylabel = "PIP (μM)", dpi=300, size = (600,400)) 
    volume = 0.5
    p=getP(oscdata, i)
    currow = oscdata[i, :]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    jumpU0 = GillespieConverter.convertU0(u0, volume)
    jumpP = GillespieConverter.convertP(p, volume)
    jumpProb = GillespieConverter.getJumpProb(fullrn, jumpU0, jumpP, (0.0, 600.0))
    jumpSol = solve(jumpProb, SSAStepper())
    GillespieConcSol = [j[1] for j in jumpSol.u] ./ (GillespieConverter.Nₐ * volume * 1e-21)
    Plots.plot!(exSolPlt, jumpSol.t[begin:100:end], GillespieConcSol[begin:100:end], label = "Gillespie Simulation", legend=:topright)
    savefig(exSolPlt, "/Users/ezragreenberg/Julia/ExperimentalFullModelWork/graphStorage/gillespiecomps/sol$i.png")
end
