using Plots
using DataFrames
using CSV
using Interact
include("Constants.jl")
include("EvaluationFunctions.jl")

#osc_data = DataFrame(CSV.File("ExperimentalFullModelWork/OscCsvs/Osc_df_1.0.csv"))

"""
Converts row in oscillatory solution dataframe to an ODESolution
- `df` Dataframe of oscillatory solutions
- `row` Row of dataframe to solve
"""
function entryToSol(df, row; tspan = shortSpan)
    currow = df[row,:]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    ka1est = currow[:ka1]
    kb1est = currow[:kb1]
    ka4est = currow[:ka4]
    ka7est = currow[:ka7]
    dfest = currow[:DF]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
    p = [x[2] for x in psym]
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
    #Plots.scatter!(x,y,z,mc=:black, ms = 2)
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