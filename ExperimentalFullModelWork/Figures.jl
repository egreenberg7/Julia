include("OutputHandling.jl")
include("../UTILITIES/GillespieConverter.jl")

oscdata = DataFrame(CSV.File("ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))
cd("ExperimentalFullModelWork/NoNegRunResultsCSVs")
alldata = getAllCSVs()
cd("../..")

oscParams = unique(oscdata[:, [:ka1, :ka4, :ka7, :kb1]])
allParams = unique(removeBadValues(alldata)[:, [:ka1, :ka4, :ka7, :kb1]])

"""
Check if I had invalid parameters...
"""
begin
    for i in size(oscdata)[1]
        curP = getP(oscdata, i)
        if curP[3] > 10.0^3 || curP[3] < 10.0^-3.0
            println("Bad kcat1")
        elseif curP[8] < 10.0 ^ -3.0 || curP[3] > 10.0^3.0 
            println("Bad kb") 
        elseif curP[10] <10.0^-3.0 || curP[3] > 10.0^3.0
            println("Another bad kb")
        end
    end
end


groupedDF = groupby(oscParams, :ka4)

function addPlt(plt, df; markershape=:circle)
    scatter!(plt, log10.(df.ka1), log10.(df.kb1), log10.(df.ka7), label="10^{$(log10(df.ka4[1]))}", markershape = markershape, alpha = 0.5)
    plt
end
function addProj(plt, df; markershape=:circle)
    numElems = size(df)[1]
    scatter!(plt, log10.(df.ka1), log10.(df.kb1), fill(zlims(plt)[1], numElems), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    scatter!(plt, log10.(df.ka1), fill(ylims(plt)[2], numElems), log10.(df.ka7), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    scatter!(plt, fill(xlims(plt)[1], numElems), log10.(df.kb1), log10.(df.ka7), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    plt
end

function addConcPlt(plt, df; markershape=:circle)
    scatter!(plt, log10.(df.K), log10.(df.P), log10.(df.A), label="10^{$(log10(df.L[1]))}", markershape = markershape, alpha = 0.5)
    plt
end
function addConcProj(plt, df; markershape=:circle)
    numElems = size(df)[1]
    scatter!(plt, log10.(df.K), log10.(df.P), fill(zlims(plt)[1], numElems), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    scatter!(plt, log10.(df.K), fill(ylims(plt)[2], numElems), log10.(df.A), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    scatter!(plt, fill(xlims(plt)[1], numElems), log10.(df.P), log10.(df.A), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
    plt
end


#Code to create zoomed in graph of oscillatory rate constants
begin 
    plt = Plots.plot(legendtitle="ka4 (μM/s)", 
        title="Dimensionality Factor = 10,000",
        legend=:outerright, dpi = 300,
        formatter=x->"10^{$(round(x, digits=1))}")
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

    xlabel!(plt, "ka1 (μM/s)")
    ylabel!(plt, "kb1 (μM/s)")
    zlabel!(plt, "ka7 (μM/s)")
    title!(plt, "Oscillatory Parameters:\n Dimensionality Factor = 10,000")
    plt
    Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/oscillatoryParams.png")
end 

#Code to create zoomed out graph of oscillatory params
begin
    plt = Plots.plot(legendtitle="ka4 (μM/s)", 
        title="Dimensionality Factor = 10,000",
        legend=:outerright, dpi = 300,
        formatter=x->"10^{$(round(x, digits=1))}",
        xlims=(log10(kaRange[1]), log10(kaRange[end])), 
        ylims = (log10(kbRange[1]), log10(kbRange[end])), 
        zlims = (log10(ka7Range[1]), log10(ka7Range[end])), 
        xtickfontsize = 6,
        ytickfontsize = 6,
        ztickfontsize = 6)
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

    xlabel!(plt, "ka1 (μM/s)")
    ylabel!(plt, "kb1 (μM/s)")
    zlabel!(plt, "ka7 (μM/s)")
    title!(plt, "Oscillatory Parameters:\n Dimensionality Factor = 10,000")
    plt
    Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/oscillatoryParamsOut.png")
end

#Code to create plot of search space 
begin
    searchPLT = Plots.scatter3d([],[],[],
        xlims=(log10(kaRange[1]), log10(kaRange[end])), 
        ylims = (log10(kbRange[1]), log10(kbRange[end])), 
        zlims = (log10(ka7Range[1]), log10(ka7Range[end])), 
        legendtitle="ka4 (μM/s)", 
        label=:none,
        legend=:topright,
        title = "Parameter Search Space: \n Dimensionality Factor = [0.1, 10,000]",
        xlabel = "ka1 (μM/s)",
        ylabel = "kb1 (μM/s)",
        zlabel = "ka7 (μM/s)",
        xtickfontsize = 6,
        ytickfontsize = 6,
        ztickfontsize = 6,
        zticks = range(0.4,1.0,4),
        dpi = 300,
        formatter = x->"10^{$x}")
    for i in log10.(kaRange)
        scatter3d!(searchPLT, [],[],[],label = "10^{$i}", shape = :circle, color=:black)
    end
    addProj(searchPLT, allParams)
    searchPLT
    Plots.savefig(searchPLT,"ExperimentalFullModelWork/graphStorage/searchspace.png")
end


#Representative parameter combination 
paramset = oscdata[oscdata[:, :ka1].==0.1 .&& oscdata[:,:kb1] .== 0.1 .&& oscdata[:, :ka7] .== 10^0.6 .&& oscdata[:, :ka4] .== 10 ^-2.5,:]
u0Graph = make3DAmpGraph(paramset)
Plots.savefig(u0Graph, "ExperimentalFullModelWork/graphStorage/u0graph.png")

sol = entryToSol(paramset, 5, tspan=600)
exSolPlt = Plots.plot(sol, title="Representative Solution",
    label="Numerical Solution", xlabel = "Time (s)", ylabel = "PIP (μM)", 
    dpi=300, size=(1800,400), legendfontsize=12)
volume = 0.5
p=getP(paramset, 5)
    currow = paramset[5, :]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    jumpU0 = GillespieConverter.convertU0(u0, volume)
    jumpP = GillespieConverter.convertP(p, volume)
    jumpProb = GillespieConverter.getJumpProb(fullrn, jumpU0, jumpP, (0.0,600.0))
    jumpSol = solve(jumpProb, SSAStepper())#, saveat=0.1)
    GillespieConcSol = [i[1] for i in jumpSol.u] ./ (GillespieConverter.Nₐ * volume * 1e-21)
    plot!(exSolPlt, jumpSol.t, GillespieConcSol, label = "Gillespie Simulation", legend=:topright)
    xlabel!(exSolPlt, "Time (s)")
savefig(exSolPlt, "ExperimentalFullModelWork/graphStorage/exampleSol5.png")
#Code to create a concentration plot
begin 
    groupedDFC = groupby(paramset, :L)
    paramset = oscdata[oscdata[:, :ka1].==0.1 .&& oscdata[:,:kb1] .== 0.1 .&& oscdata[:, :ka7] .== 10^0.6 .&& oscdata[:, :ka4] .== 10 ^-2.5,:]
    plt = Plots.plot(legendtitle="PIP (μM)", 
        title="Oscillatory Concentrations for \nRepresentative Rate Constant Set",
        legend=:topright, 
        dpi = 300,
        formatter=x->"10^{$(round(x, digits=1))}",
        xlims = log10.((Krange[1],Krange[end])),
        ylims = log10.((Prange[1], Prange[end])),
        zlims = log10.((Arange[1],Arange[end])),
        xtickfontsize = 7,
        ytickfontsize = 7,
        ztickfontsize = 7,
        xguidefontsize = 10,
        yguidefontsize = 10,
        zguidefontsize = 10
        )
    shapes = [:circle, :diamond, :utriangle]

    for (df,shape) in zip(groupedDFC, shapes)
        addConcPlt(plt, df, markershape=shape)
    end
    xbounds = xlims(plt)
    ybounds = ylims(plt)
    zbounds = zlims(plt)
    for (df,shape) in zip(groupedDFC, shapes)
        addConcProj(plt, df, markershape=shape)
        xlims!(plt, xbounds)
        ylims!(plt, ybounds)
        zlims!(plt, zbounds)
    end

    xlabel!(plt, "PIP5K (μM)")
    #ylabel!(plt, "Synaptojanin 1 (μM)")
    zlabel!(plt, "AP2 (μM)")
    plt
    Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/oscillatoryParams1.png")
end 