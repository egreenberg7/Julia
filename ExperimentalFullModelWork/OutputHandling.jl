using DifferentialEquations
using Plots
using Catalyst
using DataFrames
using CSV
using Statistics
using Peaks
using FFTW
include("Constants.jl")
include("EvaluationFunctions.jl")

osc_data = DataFrame(CSV.File("ExperimentalFullModelWork/OscCsvs/Osc_df_1.0.csv"))

"""
Converts row in oscillatory solution dataframe to an ODESolution
- `df` Dataframe of oscillatory solutions
- `row` Row of dataframe to solve
"""
function entryToSol(df, row)
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
    return solve(remake(prob, u0=u0, p=p), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
end

"""
Creates plots of all solutions and their Fourier transforms for a dataframe of oscillatory points
in directory sum_plots
- `osc_df` Dataframe of oscillatory solutions
"""
function PlotSolutions(osc_df)
    for i in 1:size(osc_df)[1]
        cursol = entryToSol(osc_df, i)
        solPlot = plot(cursol, title="L vs t")
        fftSolPlot = plot((2:length(cursol.t)/2 + 1),(abs.(rfft(cursol.u)))[2:end], title="Fourier Transform")
        plot(solPlot,fftSolPlot, layout=(2,1))
        savefig("Someplots/$i.png")
    end
end

"""
Makes 3d plots of points with oscillatory solutions, maybe eventually with colors given by period and amplitude
"""
function PlotOscillatoryRegime(df)
    for ka7 in kaRange
        boolka7 = df[:,:ka7].==ka7
        roundedka7 = round(log10(ka7),sigdigits=1)
        curdf = df[boolka7, :]
        curplot = PlotNoka7WithShadow(curdf, "log(ka7) = $(roundedka7)")
        savefig(curplot, "Someplots/ka7=$(roundedka7).png")
    end
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

    #Create series from the graph lims to create projections
    #TODO Make sure that we do not have overlapping points, set colors to be good, set marker size also
    numElements = length(x)
    xyPlane = fill(zlims(myplot)[1], numElements)
    xzPlane = fill(ylims(myplot)[2], numElements)
    yzPlane = fill(xlims(myplot)[1], numElements)

    final_plot = Plots.scatter(x, y, xyPlane, markershape=:xcross, mc=:blue, msa=0.5, ma=0.5)
    Plots.scatter!(x, xzPlane, z, markershape=:xcross,mc=:blue,msa=0.5,ma=0.5)
    Plots.scatter!(yzPlane, y, z, markershape=:xcross,mc=:blue,msa=0.5,ma=0.5)

    #Add labels
    Plots.scatter!(xlabel="log(ka1)",ylabel="log(kb1)",zlabel="log(ka4)", legend=:none, title=title)

    #Replot coordinates to make sure they are on top of projected points
    Plots.scatter!(x,y,z,mc=:black)
    return final_plot
end