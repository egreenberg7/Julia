include("../ExperimentalFullModelWork/EvaluationFunctions.jl")
using DataFrames
using CSV
using Catalyst
using Random
using Distributions
using Plots

"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) K(t) P(t) A(t) Lp(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
    # ALIASES: L = PIP, Lp = PIP2, K = Kinase (PIP5K), P = Phosphatase (Synaptojanin), A = AP2 
    # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    (ka1,kb1), L + K <--> LK # L binding to kinase
    kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
    kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

    #previously excluded reactions, all possible combinations possible in vitro
    (ka2,kb2), Lp + AK <--> LpAK
    (ka2*y,kb2), Lp + AKL <--> LpAKL
    (ka2,kb2), Lp + AP <--> LpAP
    (ka2*y,kb2), Lp + APLp <--> LpAPLp
    (ka3,kb3), A + K <--> AK
    (ka4,kb4), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3*y,kb3), LpA + LK <--> LpAKL
    (ka4*y,kb4), LpA + LpP <--> LpAPLp
    (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end  

#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
        :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :y => 3631.050539219606]
p = [x[2] for x in psym]
    
#initial condition list
usym = [:L => 0, :K => 10^-0.2895987, :P => 0.820348, :A => 10^0.42483, :Lp => 0.0, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]

#timespan for integration
shortspan = 100.0
longspan = 600.0
#solve the reduced ODEs
prob = ODEProblem(fullrn, u0, shortspan, p)

function testClassifier(numIterations; saveCSV = true, outputdirectory = "Users/ezragreenberg/Julia/Someplots/",filename="mytest.csv")
    #df = DataFrame(u0=Vector{Float64}[], p=Vector{Float64}[], retcode=Float64[], per=Float64[], amp=Float64[])
    #CSV of array gets read in as string, will have to generalize this later
    df = DataFrame(L=Float64[],K=Float64[],P=Float64[],A=Float64[], retcode=Float64[], per=Float64[], amp=Float64[])
    for i in 1:numIterations
        u0[1] = rand(Random.seed!(i),Distributions.LogUniform(0.01, 100)) #Lp #1 for L, 2 for Lp
        u0[2] = rand(Random.seed!(numIterations + i),Distributions.LogUniform(0.001, 100)) #K
        u0[3] = rand(Random.seed!(2 * numIterations + i),Distributions.LogUniform(0.01, 100)) #P
        u0[4] = rand(Random.seed!(3 * numIterations + i),Distributions.LogUniform(0.001, 100)) #A
        retcode = adaptiveSolve(prob, u0, shortspan, longspan, p)
        push!(df, Dict(:L=>u0[1],:K=>u0[2],:P=>u0[3],:A=>u0[4], :retcode => retcode[1], :per => retcode[2], :amp =>retcode[3]))
        #push!(df, Dict(:u0 => u0, :p => p, :retcode => retcode[1], :per => retcode[2], :amp =>retcode[3]))
        if i % 500 == 0 
            println("$i iterations completed")
        end
    end
    return df
end

mydf = testClassifier(10000)
booldf =  mydf[:,:retcode].!=1.1
non_terminated_df = mydf[booldf, :]
osc_df = mydf[mydf[:,:retcode] .< 0, :]

function entryToSol(df, row)
    currow = df[row,:]
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    return solve(remake(prob, u0=u0, tspan=(0,longspan)), Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
end

function PlotSolutions(osc_df, numsols = size(osc_df)[1])
    for i in 1:numsols
        cursol = entryToSol(osc_df, i)
        solPlot = plot(cursol, title="L vs t")
        fftSolPlot = plot((2:length(cursol.t)/2 + 1),(abs.(rfft(cursol.u)))[2:end], title="Fourier Transform")
        plot(solPlot,fftSolPlot, layout=(2,1))
        savefig("Someplots/$i.png")
    end
end

PlotSolutions(non_terminated_df)
    
OscilatorsInNonTerminated = [100, 102, 104, 105, 108, 109, 111, 112, 113, 115, 119, 12, 122, 123, 126, 129, 130, 135, 14, 141, 143, 144, 147, 148, 151, 152, 153, 156, 160, 161, 168, 169, 170, 173, 177, 178, 179, 18, 181, 183, 185, 19, 22, 24, 25, 36, 37, 39, 40, 41, 42, 45, 47, 48, 50, 51, 56, 57, 58, 6, 62, 63, 67, 70, 75, 77, 79, 80, 81, 84, 88, 93, 94]


CSV.write("/Users/ezragreenberg/Julia/Someplots/"*"mytest.csv", mydf)

newdf = DataFrame(CSV.File("/Users/ezragreenberg/Julia/Someplots/"*"mytest.csv"))