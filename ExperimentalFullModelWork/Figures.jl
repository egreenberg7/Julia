include("OutputHandling.jl")
include("../UTILITIES/GillespieConverter.jl")

oscdata = DataFrame(CSV.File("ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))
cd("ExperimentalFullModelWork/NoNegRunResultsCSVs")
alldata = getAllCSVs()
cd("../..")

oscParams = unique(oscdata[:, [:ka1, :ka4, :ka7, :kb1]])
allParams = unique(removeBadValues(alldata)[:, [:ka1, :ka4, :ka7, :kb1]])


groupedDF = groupby(oscParams, :ka4)

function addPlt(plt, df; markershape=:circle)
    scatter!(plt, log10.(df.ka1), log10.(df.kb1), log10.(df.ka7), label="$(log10(df.ka4[1]))", markershape = markershape, alpha = 0.5)
    plt
end
function addProj(plt, df; markershape=:circle)
    numElems = size(df)[1]
    scatter!(plt, log10.(df.ka1), log10.(df.kb1), fill(zlims(plt)[1], numElems), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    scatter!(plt, log10.(df.ka1), fill(ylims(plt)[2], numElems), log10.(df.ka7), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    scatter!(plt, fill(xlims(plt)[1], numElems), log10.(df.kb1), log10.(df.ka7), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    plt
end

plt = Plots.plot(legendtitle="log10(ka4)", 
    title="Dimensionality Factor = 10,000",
    legend=:outerright)
shapes = [:circle, :diamond, :utriangle]

for (df,shape) in zip(groupedDF, shapes)
    addPlt(plt, df, markershape=shape)
end
xbounds = xlims(plt)
ybounds = ylims(plt)
zbounds = zlims(plt)
for (df,shape) in zip(groupedDF, shapes)
    addProj(plt, df, markershape=shape)
    xlims!(plt, xbounds)
    ylims!(plt, ybounds)
    zlims!(plt, zbounds)
end

xlabel!(plt, "log10(ka1)")
ylabel!(plt, "log10(kb1)")
zlabel!(plt, "log10(ka7)")
title!(plt, "Oscillatory Parameters\n Dimensionality Factor = 10,000")
plt
Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/oscillatoryParams.png")

searchPLT = Plots.scatter3d([],[],[],
    xlims=(log10(kaRange[1]), log10(kaRange[end])), 
    ylims = (log10(kbRange[1]), log10(kbRange[end])), 
    zlims = (log10(ka7Range[1]), log10(ka7Range[end])), 
    legendtitle="log10(ka4)", 
    label=:none,
    legend=:topright;
    xlabel = "log10(ka1)",
    ylabel = "log10(kb1)",
    zlabel = "log10(ka7)",
    zticks = range(0.4,1.0,4))
for i in log10.(kaRange)
    scatter3d!(searchPLT, [],[],[],label = "$i", shape = :circle, color=:black)
end
addProj(searchPLT, allParams)
searchPLT
Plots.savefig(searchPLT,"ExperimentalFullModelWork/graphStorage/searchspace.png")

#Representative parameter combination 
paramset = oscdata[oscdata[:, :ka1].==0.1 .&& oscdata[:,:kb1] .== 0.1 .&& oscdata[:, :ka7] .== 10^0.6 .&& oscdata[:, :ka4] .== 10 ^-2.5,:]
u0Graph = make3DAmpGraph(paramset)
Plots.savefig(u0Graph, "ExperimentalFullModelWork/graphStorage/u0graph.png")

sol = entryToSol(paramset, 4, tspan=600)
exSolPlt = Plots.plot(sol)
p=getP(paramset, 4)
    currow = paramset[4, :]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    jumpU0 = GillespieConverter.convertU0(u0, 3010)
    jumpP = GillespieConverter.convertP(p, 3010)
    jumpProb = GillespieConverter.getJumpProb(fullrn, jumpU0, jumpP, (0.0,600.0))
    jumpSol = solve(jumpProb, SSAStepper())#, saveat=0.1)
    display(Plots.plot(jumpSol))
solplts = Plots.plot(exSolPlt, Plots.plot(jumpSol, idxs=1))