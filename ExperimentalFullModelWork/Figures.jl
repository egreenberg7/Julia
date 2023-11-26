include("OutputHandling.jl")
include("../UTILITIES/GillespieConverter.jl")
using StatsPlots
using GLM
using LaTeXStrings

PLOTS_DEFAULTS=(Dict(:dpi=>600))

#oscdata = DataFrame(CSV.File("/Users/ezragreenberg/JLab/Julia/ExperimentalFullModelWork/MaybeOscValuesAnalysis/AllExpOsc.csv"))
oscdata = CSV.read("/Users/ezragreenberg/Documents/GitHub/Julia/ExperimentalFullModelWork/paramsWithMinDF.csv", DataFrame)

# cd("/Users/ezragreenberg/JLab/Julia/ExperimentalFullModelWork/NoNegRunResultsCSVs")
# alldata = getAllCSVs()
# cd("../..")

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

#TODO Add correct legend labels, move to 1 plot with 4 subplots
function solutionsGroupedBySpecies(row; df = oscdata, tspan = 600)
    speciesIndices = Dict([:L => 1, :K => 2, :P => 3, :A => 4, :Lp => 5, :LpA => 6, :LK => 7, 
        :LpP => 8, :LpAK => 9, :LpAP => 10, :LpAKL => 11, :LpAPLp => 12, :AK => 13, :AP => 14, 
        :AKL => 15, :APLp => 16])
    #L, Lp
    lipids = [1, 2]
    #A, LpA
    adaptors = [4, 6]
    kinases = [speciesIndices[i] for i in [:K, :LK, :LpAK, :LpAKL, :AK, :AKL]]
    phosphatases = [speciesIndices[i] for i in [:P, :LpP, :LpAP, :LpAPLp, :AP, :APLp]]
    speciesTypes = [lipids, adaptors, kinases, phosphatases]
    solutions = [entryToSol(oscdata, row; tspan = tspan, save_idxs = i) for i in speciesTypes]
    plots = [plot(i) for i in solutions]
    #plot(plots, layout=4) 
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
    Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/trulyOscillatoryParams.png")
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
    Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/trulyOscillatoryParamsOut.png")
end

#Code to create plot of search space of parameters
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

#Code to create plot of primary search space for concentrations
begin
    include("ConcentrationClassification.jl")
    r = getu0Ranges()
    searchPLT = Plots.scatter3d([],[],[],
        xlims=(log10(Krange[1]), log10(Krange[end])), 
        ylims = (log10(Prange[1]), log10(Prange[end])), 
        zlims = (log10(Arange[1]), log10(Arange[end])), 
        legendtitle="PIP (μM)", 
        label=:none,
        legend=:topright,
        title = "Concentration Primary Search Space",
        xlabel = "PIP5K (μM)",
        #ylabel = " (μM)",
        zlabel = "AP2 (μM)",
        xtickfontsize = 6,
        ytickfontsize = 6,
        ztickfontsize = 6,
        dpi = 300,
        formatter = x->"10^{$x}")
    #for i in log10.(r[:L])
        scatter3d!(searchPLT, [],[],[],label = "10^{-2.0} to 10.^{2.0} by Increments of √10", shape = :circle, color=:black)
    #end
    cRangeDict = Dict([:L => Lrange, :K=>Krange,:A=>Arange,:P=>Prange])
    addConcProj(searchPLT, DataFrame(getConcentrationCombos(cRangeDict)))
    searchPLT
    Plots.savefig(searchPLT,"ExperimentalFullModelWork/graphStorage/Concentrationsearchspace1.png")
end

#Code to create plot of secondary search space of concentrations
begin
    include("ConcentrationClassification.jl")
    r = getu0Ranges()
    searchPLT = Plots.scatter3d([],[],[],
        xlims=(log10(r[:K][1]), log10(r[:K][end])), 
        ylims = (log10(r[:P][1]), log10(r[:P][end])), 
        zlims = (log10(r[:A][1]), log10(r[:A][end])), 
        legendtitle="PIP (μM)", 
        label=:none,
        legend=:topright,
        title = "Concentration Secondary Search Space",
        xlabel = "PIP5K (μM)",
        #ylabel = " (μM)",
        zlabel = "AP2 (μM)",
        xtickfontsize = 6,
        ytickfontsize = 6,
        ztickfontsize = 6,
        dpi = 300,
        formatter = x->"10^{$x}")
    #for i in log10.(r[:L])
        scatter3d!(searchPLT, [],[],[],label = "10^{-2.0} to 10.^{2.0} by Increments of 10^{0.2}", shape = :circle, color=:black)
    #end
    addConcProj(searchPLT, getConcentrationCombos(getu0Ranges()))
    searchPLT
    Plots.savefig(searchPLT,"ExperimentalFullModelWork/graphStorage/Concentrationsearchspace.png")
end 

#Code to create K vs A very nicely
begin
    regLine = lm(@formula(log10(K)~log10(A)), oscdata)
    scatter(log10.(oscdata.K), log10.(oscdata.A), 
        xaxis="AP2 (μM)", yaxis="PIP5K (μM)", smooth = true,
        formatter=x->"10^{$x}",
        label=:none,
        color=:black,
        legendfontsize=10,
        title="Initial Kinase and Adaptor Concentrations")
    plot!([],[],ms=0,color="white",label=L"PIP5K = 1.46 * AP2^{0.955}")
    plot!([],[],ms=0,color="white",label=L"r^2 = 0.933", 
    dpi = 600)
    savefig("ExperimentalFullModelWork/graphStorage/KvsA.png")
end

#Code to create L vs A nicely NOT REALLY CLEAR CORRELATION
begin
    regLine = lm(@formula(log10(L)~log10(A)), oscdata)
    scatter(log10.(oscdata.L), log10.(oscdata.A), 
        xaxis="AP2 (μM)", yaxis="PIP2 (μM)", smooth = true,
        formatter=x->"10^{$x}",
        label=:none,
        color=:black,
        legendfontsize=10,
        title="Initial Lipid and Adaptor Protein Concentrations")
    plot!([],[],ms=0,color="white",label=L"PIP5K = 1.46 * AP2^{0.955}")
    plot!([],[],ms=0,color="white",label=L"r^2 = 0.933", 
    dpi = 600)
    #savefig("ExperimentalFullModelWork/graphStorage/LvsA.png")
end

#Code to plot P vs Period very nicely
begin
    regLine = lm(@formula(log10(per)~log10(P)), oscdata)
    constantMultiplier = round(10^coef(regLine)[1],digits=3)
    exponent = round(coef(regLine)[2], digits = 3)
    correlation = round(r2(regLine), digits = 3)
    residualPlot=scatter(log10.(oscdata.P), predict(regLine, oscdata).-log10.(oscdata.per),ylims=(-0.4,0.4))
    scatter(log10.(oscdata.per), log10.(oscdata.P), 
        xaxis="Synaptojanin (μM)", yaxis="Period (s)", smooth = true,
        formatter=x->"10^{$x}",
        label=:none,
        linecolor=:black,
        mc=:blue,
        legendfontsize=10,
        title="Period vs Initial Phosphatase Concentration",
        linewidth = 3,
        linestyle=:dash,
        dpi=600)
    plot!([],[],ms=0,color="white",label=L"Period = %$constantMultiplier * Synaptojanin^{%$exponent}")
    plot!([],[],ms=0,color="white",label=L"r^2 = %$correlation")
    savefig("ExperimentalFullModelWork/graphStorage/PerVsP.png")
end

begin
    regLine = lm(@formula(log10(per)~log10(ka4)), oscdata)
    constantMultiplier = round(10^coef(regLine)[1],digits=3)
    exponent = round(coef(regLine)[2], digits = 3)
    correlation = round(r2(regLine), digits = 3)
    residualPlot=scatter(log10.(oscdata.ka4), predict(regLine, oscdata).-log10.(oscdata.per),ylims=(-0.4,0.4))
    scatter(log10.(oscdata.per), log10.(oscdata.ka4), 
        xaxis="Rate Constant for Binding of Synaptojanin to AP2 1/(μM*̇s)", yaxis="Period (s)", smooth = true,
        formatter=x->"10^{$(round(x,digits=2))}",
        label=:none,
        linecolor=:black,
        mc=:blue,
        legendfontsize=10,
        title="Period vs Binding Rate of Phosphatase to Adaptor",
        linewidth = 3,
        linestyle=:dash,
        dpi = 600)
    plot!([],[],ms=0,color="white",label=L"Period = %$constantMultiplier * k_{on}^{%$exponent}")
    plot!([],[],ms=0,color="white",label=L"r^2 = %$correlation")
    savefig("ExperimentalFullModelWork/graphStorage/PerVska4.png")
end

#Code to create plot of oscillatory concentrations
x = makeConcentrationGraph(oscdata)
savefig(x, "ExperimentalFullModelWork/graphStorage/OscillatoryConcentrations.png")

#Representative parameter combination 
paramset = oscdata[oscdata[:, :ka1].==0.1 .&& oscdata[:,:kb1] .== 0.1 .&& oscdata[:, :ka7] .== 10^0.6 .&& oscdata[:, :ka4] .== 10 ^-2.5,:]
#u0Graph = make3DAmpGraph(paramset)
#Plots.savefig(u0Graph, "ExperimentalFullModelWork/graphStorage/u0graph.png")

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
    #Plots.savefig(plt, "ExperimentalFullModelWork/graphStorage/oscillatoryParams1.png")
end 


#Miscelanneous statistical plots
scatter(log10.((oscdata.ka1 .* Km1exp .- oscdata.kb1)),(oscdata.L))
scatter((oscdata.K), (oscdata.P))

@df oscdata scatter(:ka1, :ka4)
using PlotlyJS
PlotlyJS.plot(oscdata,y=:A,kind="box",Layout(yaxis_type="log"))

num = 900
Plots.plot(entryToSol(oscdata,num;tspan=5000,save_idxs=[10,11]),labels=["LpAP" "LpAKL"])
Plots.plot(entryToSol(oscdata, num; tspan=5000, save_idxs=[1,2,5,6,10,11]),labels=["L" "K" "Lp" "LpA" "LpAP" "LpAKL"], xlims=(100,1000), ylims=(0,1))
Plots.plot(entryToSol(oscdata,num;tspan=1000, save_idxs=[1,2,5,6,10,11]),labels=["L" "K" "Lp" "LpA" "LpAP" "LpAKL"])
# for i in 1:1008
#     x = Plots.plot(entryToSol(oscdata, i; tspan=5000, save_idxs=[1,2,3]))
#     Plots.savefig(x,"graphstorage2/OscdataSolution$i.png")
# end

