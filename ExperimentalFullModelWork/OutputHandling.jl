using Plots
using DataFrames
using CSV
using Interact
using FFTW
using LaTeXStrings
include("NoNegConstants.jl")
include("EvaluationFunctions.jl")

#osc_data = DataFrame(CSV.File("ExperimentalFullModelWork/OscCsvs/Osc_df_1.0.csv"))

"""
Converts row in oscillatory solution dataframe to an ODESolution
- `df` Dataframe of oscillatory solutions
- `row` Row of dataframe to solve
"""
function entryToSol(df, row; tspan = shortSpan, save_idxs = 1)
    currow = df[row,:]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    ka1est = currow[:ka1]
    kb1est = currow[:kb1]
    ka4est = currow[:ka4]
    ka7est = currow[:ka7]
    dfest = currow[:df]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
    p = [x[2] for x in psym]
    return solve(remake(prob, u0=u0, tspan=(0.0, tspan), p=p), Rodas4(), abstol=1e-8, reltol=1e-12, saveat=0.1, save_idxs=save_idxs, maxiters=200 * tspan, verbose=false, callback=TerminateSteadyState(1e-8, 1e-12))
end

"""
Converts row in dataframe of oscillatory concentrations for a given p 
into a solution.
- `df` Dataframe of oscillatory solutions
- `row` Row of dataframe to solve
- `p` Rate constants for the given dataframe
"""
function entryToSol(df, row, p; tspan = shortSpan)
    currow = df[row,:]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A] 
    return solve(remake(prob, u0=u0, tspan=(0.0, tspan), p=p), Rodas4(), abstol=1e-8, reltol=1e-12, saveat=0.1, save_idxs=1, maxiters=200 * tspan, verbose=false, callback=TerminateSteadyState(1e-8, 1e-12))
end

"""
Creates plots of all solutions and their Fourier transforms for a dataframe of oscillatory points
in directory sum_plots
- `osc_df` Dataframe of oscillatory solutions
"""
function PlotSolutions(osc_df, numsols = size(osc_df)[1])
    for i in 1:numsols
        cursol = entryToSol(osc_df, i)
        solPlot = plot(cursol, title="L vs t")
        fftSolPlot = plot((2:length(cursol.t)/2 + 1),(abs.(rfft(cursol.u)))[2:end], title="Fourier Transform")
        plot(solPlot,fftSolPlot, layout=(2,1))
        savefig("Someplots/$i.png")
    end
end

"""
Rewrites CSV to only have one df value; useful for initial implementation where CSVs just appended new values as DF increased
Also good to quickly read in one of the later files.
- `dfval` df value that dataframe should have 
"""
function isolateDF(dfval)
    my_df = DataFrame(CSV.File("ExperimentalFullModelWork/OscCsvs/Osc_df_$(dfval).csv"))
    booldf =  my_df[:,:DF].==dfval
    isolated_df = my_df[booldf, :]
    return isolated_df
end

""""
Overwrites alllcsv files with isolatedDF files.
"""
function isolateDFfiles(dfRange=dfRange)
    for df in dfRange
        CSV.write("ExperimentalFullModelWork/OscCsvs/Osc_df_$(df).csv",isolateDF(df))
    end
end

"""
Makes 3d plots of points with oscillatory solutions by looping over ka7,
maybe eventually with colors given by period and amplitude; returns dictionary of ka7 values mapped to plots
- `df` Dataframe from ExperimentalOscillators output 
- `kaRange=kaRange` kaRange used to generate data being analyzed
- `toSave=false` Set to true if you want to save plots to file
"""
function PlotOscillatoryRegime(df, kaRange = kaRange, toSave=false)
    my_plots = Dict{Float64, Plots.Plot{Plots.GRBackend}}()#(undef, length(kaRange))
    for ka7 in kaRange
        boolka7 = df[:,:ka7].==ka7
        roundedka7 = round(log10(ka7),sigdigits=1)
        curdf = df[boolka7, :]
        curplot = PlotNoka7WithShadow(curdf, "log(ka7) = $(roundedka7)")
        my_plots[ka7] = curplot
        if toSave
            savefig(curplot, "Users/ezragreenberg/Julia/ExperimentalFullModelWork/Someplots/ka7=$(roundedka7).png")
        end
    end
    return my_plots
end

"""
Plots a 3d graph of dataframe entries with log parameter values for ka1, kb1, and ka4 as the axes.
    Values are projected onto the coordinate planes
- `osc_df` Dataframe with values to be plotted
- `title` Title to be given to plot
"""
function PlotNoka7WithShadow(df, title)
    #Set coordinate values and plot them
    x = log10.(df[:,:ka1])
    y = log10.(df[:,:kb1])
    z = log10.(df[:,:ka4])

    #Plot coordinates to find limits for projections
    myplot = Plots.scatter(x,y,z)

    #Make it so there are set axis sizes based on the ranges, can be commented out
    xlims!((log10.((kaRange[1], kaRange[end]))))
    ylims!(log10.((kbRange[1], kbRange[end])))
    zlims!(log10.((kaRange[1], kaRange[end])))

    #Create series from the graph lims to create projections
    numElements = length(x)
    xyPlane = fill(zlims(myplot)[1], numElements)
    xzPlane = fill(ylims(myplot)[2], numElements)
    yzPlane = fill(xlims(myplot)[1], numElements)

    final_plot = Plots.scatter(x, y, xyPlane, markershape=:xcross, mc=:blue, msa=0.5, ma=0.5, ms = 1)
    Plots.scatter!(x, xzPlane, z, markershape=:xcross,mc=:blue,msa=0.5,ma=0.5, ms = 1)
    Plots.scatter!(yzPlane, y, z, markershape=:xcross,mc=:blue,msa=0.5,ma=0.5, ms = 1)

    #Add labels
    Plots.scatter!(xlabel="log(ka1)",ylabel="log(kb1)",zlabel="log(ka4)", legend=:none, title=title)

    #Replot coordinates to make sure they are on top of projected points
    Plots.scatter!(x,y,z,mc=:black, ms = 2)
    #Make it so there are set axis sizes based on the ranges, can be commented out
    xlims!((log10.((kaRange[1], kaRange[end]))))
    ylims!(log10.((kbRange[1], kbRange[end])))
    zlims!(log10.((kaRange[1], kaRange[end])))

    return final_plot
end

"""
Uploads all csv files of oscillatory parameters generated by ExperimentalOscillators into a dictionary of df values mapped to DataFrames
- `dfRange` dfRange of ExperimentalOscillators run for CSV files being extracted, default is current dfRange in Constants.jl
"""
function uploadOscCsvs(dfRange=dfRange)
    osc_dfs = Dict{Float64, DataFrame}()
    for i in dfRange
        fileName = "/Users/ezragreenberg/Julia/ExperimentalFullModelWork/OscCsvs/Osc_df_$i.csv"
        osc_dfs[i] = DataFrame(CSV.File(fileName))
    end
    return osc_dfs
end

"""
Uploads all csv files of not oscillatory parameters generated by ExperimentalOscillators into an array of DataFrames
- `dfRange` dfRange of ExperimentalOscillators run for CSV files being extracted, default is current dfRange in Constants.jl
"""
function uploadNonOscCsvs(dfRange=dfRange)
    osc_dfs = Array{DataFrame}(undef, length(dfRange))
    for i in enumerate(dfRange)
        fileName = "/Users/ezragreenberg/Julia/ExperimentalFullModelWork/OscCsvs/No_osc_df_$(i[2]).csv"
        osc_dfs[i[1]] = DataFrame(CSV.File(fileName))
    end
    return osc_dfs
end

"""
Uploads all csv files of oscillatory parameters and generates a dictionary with df and ka7 values keys to plots of oscillatory regimes
"""
function makeOscCsvPlots(dfRange=dfRange, kaRange=kaRange)
    my_dfs = uploadOscCsvs(dfRange)
    plot_arr = Dict{Float64, Dict{Float64, Plots.Plot{Plots.GRBackend}}}()
    for key in keys(my_dfs)
        ka7_dict = PlotOscillatoryRegime(my_dfs[key], kaRange)
        plot_arr[key] = ka7_dict
    end
    return plot_arr
end

"""
Makes interactive graph with slider to visualize all CSV data
"""
function oscPlotWithSlider(plot_dict; dfRange = dfRange, kaRange = kaRange)
    @manipulate for df=dfRange, ka7=kaRange
        Plots.plot(plot_dict[df][ka7])
    end
end


"""
Upload all csvs in the working directory into 1 dataframe. The working
directory must contain only CSV files of output from running
NoNegConstants.jl
"""
function getAllCSVs()
    mydf = DataFrame(Dict(:df => Float64[], :ka7 => Float64[],
    :ka4 => Float64[], :ka1 => Float64[], :kb1 => Float64[],
    :fit => Float64[], :per => Float64[], :amp => Float64[],
    :L => Float64[], :K => Float64[], :P => Float64[], :A => Float64[],
    :oscFound => Int64[]))          
    filepaths = readdir(;join=true)
    for i in filepaths
        if occursin("csv", i)
            x = DataFrame(CSV.File(i))
            append!(mydf, x)
        end
    end
    return mydf
end

"""
Remove all evaluated values that would result in a negative rate constant
Also remove empty rows
"""
function removeBadValues(mydf; minka7 = (kcat7exp / Km7exp), Km1 = Km1exp)
    lowka7df = mydf[:, :ka7]  .> minka7
    newdf = mydf[lowka7df,:]
    lowkcat1df = newdf[:, :ka1] .* Km1 .> newdf[:, :kb1]
    newdf = newdf[lowkcat1df, :]
    emptyrowremoval = newdf[:,:df] .!= 0.0
    newdf = newdf[emptyrowremoval, :]
    return newdf
end
"""
Get dataframe holding only the values with at least minOscSols oscillatory
solutions found
"""
function getOscValues(mydf; minOscSols = 3)
    booldf = mydf[:,:oscFound ] .>= minOscSols
    return mydf[booldf,:]
end
"""
Extracts and calculates the array of rate constants from a dataframe
"""
function getP(mydf, row)
    ka1est = mydf[row, :ka1]
    kb1est = mydf[row, :kb1]
    ka7est = mydf[row, :ka7]
    ka4est = mydf[row, :ka4]
    dfest = mydf[row, :df]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
        :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
        :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
    p = [x[2] for x in psym]
end

"""
Takes CSV file and constructs an array of all the p values contained inside. Useful
for parsing CSV with osc values to do more indepth concentration search.
"""
function getAllUniqueP(mydf)
    #Extract dataframe containing only unique p combinations for non-experimental values
    pDF = unique(mydf[:,[:ka1, :kb1, :ka7,:ka4, :df]])
    numRows = size(pDF,1)
    pArr = zeros(numRows, 13)
    #Construct the rest of the p values based on the experimental values
    for i in 1:size(pDF, 1)
        pArr[i, :] = getP(pDF, i)
    end
    return pArr    
end

"""
Take dataframe of concentrations from concentration classifier,
the parameters, and a lipid concentration and generate graph
"""
function make3DAmpGraph(mydf,L)
    booldf = mydf[:,:L] .== L
    df = mydf[booldf, :]
    x = df.K #K
    y = df.P #P 
    z = df.A #A
    constantTwos=fill(-2,length(x))
    constantThrees=fill(-3,length(x))
            
    log_x = log10.(x)
    log_y = log10.(y)
    log_z = log10.(z)

    colorVals = df.amp ./ L # NormalizedAmplitude

    graph = Plots.plot(Plots.scatter(log_x,log_y,constantTwos,marker_z=colorVals,markershape= :xcross,alpha=0.5,ms=2))
    Plots.scatter!(graph, log_x,-constantTwos,log_z,marker_z=colorVals,markershape= :xcross, alpha=0.5,ms=2)
    Plots.scatter!(graph, constantThrees,log_y,log_z,markershape= :xcross,marker_z=colorVals,alpha=0.5,ms=2)
    Plots.scatter!(graph, constantTwos,log_y,log_z,marker_z=colorVals, markershape= :xcross,alpha=0.5,ms=2)
    Plots.scatter!(graph,
        log_x,
        log_y,
        log_z,
        marker_z=colorVals,
        title=L"$10^{%$(log10(L))}\textrm{ μM PIP}$",
        titlefontsize = 14,
        xlims=(-2,2),
        ylims=(-2,2),
        zlims=(-2,2),
        legend=:none,
        markerstrokealpha=0,
        markersize = 3,
        markerstrokewidth = 0.2,
        xaxis=(L"\textrm{\log(PIP5K) (μM)}"),
        yaxis=(L"\textrm{\log(Synaptojanin) (μM)}"),
        zaxis=(L"\textrm{\log(AP2) (μM)}"),
        xguidefontsize=12,
        yguidefontsize = 12,
        zguidefontsize=12)#,
        #xticks=(-2:2,[L"$10^{-2}$",L"$10^{-1}$",L"$10^{0}$",L"$10^{1}$",L"$10^{2}$"]), 
        #zticks=(-2:2,[L"$10^{-2}$",L"$10^{-1}$",L"$10^{0}$",L"$10^{1}$",L"$10^{2}$"]),
        #yticks=(-2:2,[L"$10^{-2}$",L"$10^{-1}$",L"$10^{0}$",L"$10^{1}$",L"$10^{2}$"]))
end

"""
Take dataframe of concentrations and make graph where colors represent lipids.
"""
function makeConcentrationGraph(mydf)
    df = mydf
    x = df.K #K
    y = df.P #P 
    z = df.A #A
    constantTwos=fill(-2,length(x))
    constantThrees=fill(-3,length(x))
            
    log_x = log10.(x)
    log_y = log10.(y)
    log_z = log10.(z)

    colorVals = log10.(df.L)
    
    graph = Plots.plot(Plots.scatter(log_x,log_y,constantTwos,marker_z=colorVals,markershape= :xcross,alpha=0.5,ms=2))
    Plots.scatter!(graph, log_x,-constantTwos,log_z,marker_z=colorVals,markershape= :xcross, alpha=0.5,ms=2)
    Plots.scatter!(graph, constantThrees,log_y,log_z,markershape= :xcross,marker_z=colorVals,alpha=0.5,ms=2)
    Plots.scatter!(graph, constantThrees,log_y,log_z,marker_z=colorVals, markershape= :xcross,alpha=0.5,ms=2)
    Plots.scatter!(graph,
        log_x,
        log_y,
        log_z,
        marker_z=colorVals,
        formatter=x->"10^{$x}",
        colorbar_formatter=x->"10^{$x}",
        #colorbar_discrete_values = unique(colorVals), 
        #color = palette(:thermal, length(unique(colorVals))),
        title="Oscillatory Initial Concentrations",
        titlefontsize = 14,
        xlims=(-3,2),
        ylims=(-3,2),
        zlims=(-2,2),
        legend=:none,
        markerstrokealpha=0,
        markersize = 3,
        markerstrokewidth = 0.2,
        xaxis="PIP5K (μM)",
        #yaxis="Synaptojanin (μM)",
        zaxis="AP2 (μM)",
        dpi = 600,
        xguidefontsize=12,
        yguidefontsize = 12,
        zguidefontsize=12)#,
        #colorbar_title=("PIP (μM)")),
        #xticks=(-2:2,[L"$10^{-2}$",L"$10^{-1}$",L"$10^{0}$",L"$10^{1}$",L"$10^{2}$"]), 
        #zticks=(-2:2,[L"$10^{-2}$",L"$10^{-1}$",L"$10^{0}$",L"$10^{1}$",L"$10^{2}$"]),
        #yticks=(-2:2,[L"$10^{-2}$",L"$10^{-1}$",L"$10^{0}$",L"$10^{1}$",L"$10^{2}$"]))
end

"""
Function to plot 3d group graphs with legend of labels??? Currently Only for L
"""
function makeGroupGraph(df, labels)
    plt = Plots.scatter3d(xlabel="PIP5K (μM)", ylabel="Synaptojanin (μM)", zlabel="AP2 (μM)", 
                            title="Oscillatory Concentrations", legendtitle="PIP2 (μM)",
                            xlims = log10.((Krange[1], Krange[end])), ylims = log10.((Prange[1], Prange[end])),
                            zlims = log10.((Arange[1], Arange[end])), clims = log10.((Lrange[1], Lrange[end])),
                            formatter=x->"10^{$(round(x, digits=1))}",
                            colorbar_formatter=x->"10^{$(round(x, digits=1))}", 
                            cbartitle="PIP2 Concentrations (μM)",
                            dpi=300)
    groupedDF = groupby(df, labels)
    for i in groupedDF
        #Lval = log10(i.L[1])
        logL = log10.(i.L)
        logx = log10.(i.K)
        logy = log10.(i.P)
        logz = log10.(i.A)
        Plots.scatter3d!(plt, logx, logy, logz, marker_z = logL, labels = :none)
    end
    function addConcProj(plt, df; markershape=:circle)
        numElems = size(df)[1]
        scatter!(plt, log10.(df.K), log10.(df.P), fill(zlims(plt)[1], numElems), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
        scatter!(plt, log10.(df.K), fill(ylims(plt)[2], numElems), log10.(df.A), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
        scatter!(plt, fill(xlims(plt)[1], numElems), log10.(df.P), log10.(df.A), label=:none, markershape = markershape, alpha = 0.5, color = :black, ms = 2)
        plt
    end
    addConcProj(plt, oscdata)
    return plt
end

"""
Same as above but with Period not L
"""
function makePerGraph(df)
    plt = Plots.scatter3d(xlabel="PIP5K (μM)", ylabel="Synaptojanin (μM)", zlabel="AP2 (μM)", 
                            title="Periods of Oscillatory Solutions",
                            xlims = log10.((Krange[1], Krange[end])), ylims = log10.((Prange[1], Prange[end])),
                            zlims = log10.((Arange[1], Arange[end])), clims = (0,2),
                            formatter=x->"10^{$(round(x, digits=1))}",
                            colorbar_formatter=x->"10^{$(round(x, digits=1))}", 
                            #cbartitle="log(Period (s))",
                            dpi=300)
        #Lval = log10(i.L[1])
        logL = log10.(df.per)
        logx = log10.(df.K)
        logy = log10.(df.P)
        logz = log10.(df.A)
        Plots.scatter3d!(plt, logx, logy, logz, marker_z = logL, labels = :none)
    function addConcProj(plt, df; markershape=:circle)
        numElems = size(df)[1]
        scatter!(plt, log10.(df.K), log10.(df.P), fill(zlims(plt)[1], numElems), label=:none, markershape = markershape, alpha = 0.5, marker_z = logL, ms = 2)
        scatter!(plt, log10.(df.K), fill(ylims(plt)[2], numElems), log10.(df.A), label=:none, markershape = markershape, alpha = 0.5, marker_z = logL, ms = 2)
        scatter!(plt, fill(xlims(plt)[1], numElems), log10.(df.P), log10.(df.A), label=:none, markershape = markershape, alpha = 0.5, marker_z = logL, ms = 2)
        plt
    end
    addConcProj(plt, oscdata)
    Plots.gr_update_colorbar!
    return plt
end