include("OutputHandling.jl")
include("../UTILITIES/GillespieConverter.jl")

#oscdata = DataFrame(CSV.File("/Users/ezragreenberg/Julia/ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))
oscdata = DataFrame(CSV.File("ExperimentalFullModelWork/TrulyOscillatoryData.csv"))

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

for i in 1:size(oscdata)[1]
    sol = entryToSol(oscdata, i, tspan=5000.0)
    exSolPlt = Plots.plot(sol, title="Representative Solution",
        label="Numerical Solution", xlabel = "Time (s)", ylabel = "PIP (μM)", dpi=300, size = (600,400)) 
    savefig(exSolPlt, "/Users/ezragreenberg/JLab/Julia/graphStorage/numSolutions/sol$i.png")
end

    i = 700
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
    x=Plots.plot!(exSolPlt, jumpSol.t[begin:100:end], GillespieConcSol[begin:100:end], label = "Gillespie Simulation", legend=:topright)
    #savefig(exSolPlt, "/Users/ezragreenberg/Julia/ExperimentalFullModelWork/graphStorage/gillespiecomps/sol$i.png")
plot(jumpSol)

plotL = plot(jumpSol.L)

function solutionsGroupedBySpecies(sol; df = oscdata, tspan = 600, xlims = (0, tspan), title="Example Oscillator")
    speciesIndices = Dict([:L => 1, :K => 2, :P => 3, :A => 4, :Lp => 5, :LpA => 6, :LK => 7, 
        :LpP => 8, :LpAK => 9, :LpAP => 10, :LpAKL => 11, :LpAPLp => 12, :AK => 13, :AP => 14, 
        :AKL => 15, :APLp => 16])
    #L, Lp
    lipids = [1, 5]
    #A, LpA
    adaptors = [4, 6]
    kinases = [speciesIndices[i] for i in [:K, :LK, :LpAK, :LpAKL, :AK, :AKL]]
    phosphatases = [speciesIndices[i] for i in [:P, :LpP, :LpAP, :LpAPLp, :AP, :APLp]]

    Llabels = ["L" "Lp"]
    Alabels = ["A" "LpA"]
    Klabels = ["K" "LK" "LpAK" "LpAKL" "AK" "AKL"]
    Plabels = ["P" "LpP" "LpAP" "LpAPLp" "AP" "APLp"]
    
    speciesLabels = [Llabels, Alabels, Klabels, Plabels]

    speciesTypes = [lipids, adaptors, kinases, phosphatases]

    solutions = [solutions[i] for i in [1,3,2,4]]
    speciesLabels = [speciesLabels[i] for i in [1,3,2,4]]

    plots = [plot(i[1], labels = i[2]) for i in zip(solutions, speciesLabels)]
    for i in ["ka1", "kb1", "ka4", "ka7", "per"]
        println(i * ": " * string(round(oscdata[row, i], sigdigits=2)))
    end
    plot(plots..., layout=(2,2), plot_title=title,
    xlims = xlims,
    #xaxis="Time (s)", yaxis="Concentration (μM)", 
    dpi = 600,
    xaxis="")
    #, xlims=(100, tspan), ylims=(0,oscdata[:L, row]) 
end