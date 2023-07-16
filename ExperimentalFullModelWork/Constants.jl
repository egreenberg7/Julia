using Catalyst

#Known parameters
const Km1exp = 15.0
const ka2exp = 0.7*10^-3
const kb2exp = 2.0 * 10^-3
const ka3exp = 0.118
const kb3exp = 0.0609
const Kd4exp =29.0
const Km7exp = 39.0
const kcat7exp = 85.3

#Estimates for unknown parameters
ka1est = 0.9 * 10^-3 #Slightly less than ka2
kb1est = 0.2 * 10^-3
ka4est = 0.12 #Estimating LpA bindking to K (ka3) similar to LpA binding to L
ka7est = 0.9 #Pretty much guess, around ka1
dfest = 750

#parameter list Changed this around
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
        :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
        :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
p = [x[2] for x in psym]

    
#initial condition list
usym = [:L => 5, :K => 0.1, :P => 0.1, :A => 0.1, :Lp => 0., :LpA => 0.0, :LK => 0., 
        :LpP => 0., :LpAK => 0., :LpAP => 0., :LpAKL => 0., :LpAPLp => 0., :AK => 0., :AP => 0., 
        :AKL => 0., :APLp => 0.]
u0 = [x[2] for x in usym]

#Reaction network
const fullrn = @reaction_network fullrn begin
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

#Time bounds for ODE problem #// TODO make adaptive integration range
const shortSpan = (0.0, 100.0)
const longSpan = (0.0, 600.0)

#Set up ODE problem
const prob = ODEProblem(fullrn, u0, shortSpan, p)


#RANGES THAT VARY BY POWER OF 10
#// TODO DIMENSIONALITY CAN GO TO 10,000
dfRange = 10.0 .^ ((-2:8) ./2) #exponential range from 0.1 to 1000
kaRange = 10.0 .^ ((-6:2) ./2)
kbRange = 10.0 .^ ((-6:6) ./2)

#TODO Double resolution
Lrange = 10.0 .^ (-2:2)
Krange = 10.0 .^ (-3:2)
Prange = 10.0 .^ (-3:2)
Arange = 10.0 .^ (-2:2)

#=
#3 times resolution as above
dfRange = 10.0 .^ ((-3:9) ./3) #exponential range from 0.1 to 1000
kaRange = 10.0 .^ ((-3:9) ./3)
kbRange = 10.0 .^ ((-3:9) ./3)

Lrange = 10.0 .^ ((-6:6)./3)
Krange = 10.0 .^ ((-9:6)./3)
Prange = 10.0 .^ ((-9:6)./3)
Arange = 10.0 .^ ((-6:6)./3)
=#

#Set up array of initial conditions that will be tested in advance for ease of later code
function makeCombos(Lrange, Krange, Prange, Arange)
    arr = Array{Float64}(undef, length(Lrange)*length(Krange)*length(Arange)*length(Prange),4)
    count = 1
    for L in Lrange
        for K in Krange
            for P in Prange
                for A in Arange
                    arr[count, 1] = L
                    arr[count, 2] = K
                    arr[count, 3] = P
                    arr[count, 4] = A
                    count += 1
                end
            end
        end
    end
    return arr
end

u0combos = makeCombos(Lrange, Krange, Prange, Arange)

numU0combos = length(u0combos[:,1])

#TODO Set to const once ranges finalized

