#TODO Currently, we have it so that LpAK_0=K_0=sum(both). Need to find more realistic initial concentrations.
#TODO Also, we have to find what regimes the MM will hold on the 2D in addition to just the regular set-up.

"""
Michaelis-Menten approximation of the original oscillator model, without all possible pairs of reactions
Note that V_max=kcat*E. Because the velocity of the reaction is given by V_max[S]/(Km+S),
we can rewrite this as kcat*E0*S/(Km+S), and in the code for the rate equations,
we then divide this by Eß*S and use mass conservation to find E0 since our equation is interpreted as
 a typical rate equation.
"""
mm_rn = @reaction_network mm_rn begin
    @parameters kcat1 Km1 ka2 kb2 ka3 kb3 ka4 kb4 kcat7 Km7 DF
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LpAK(t) LpAP(t) 

    kcat1*(K + LpAK)/(K*(Km1+L)), L + K --> Lp + K # L phosphorylation by kinase into Lp using Michaelis-Menten
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (kcat1*(K + LpAK) * DF)/(LpAK*(Km1+L)), LpAK + L --> Lp + LpAK # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A) using Michaelis-Menten
    (kcat7*(P+LpAP))/((Km7+Lp) * P), Lp + P --> L + P # Lp dephosphorylation by phosphatase using Michaelis-Menten
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    DF * (kcat7*(P+LpAP))/((Km7+Lp) * LpAP), Lp + LpAP --> L + LpAP # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A) using Michaelis-Menten
end

psym=[0, 15, 0.0007, 0.002, 0.118, 0.0609, 0, 0, 39, 85.3, 1500] #0s are unknown, 1500 put in as nominal value for df based on Jonathan

#=
p = rand(11)
u = rand(8)
tspan = (0,100)
using DifferentialEquations
mmprob = ODEProblem(mm_rn,u,tspan,p)
x = convert(ODESystem, mm_rn)
for i in x.eqs println(i) end

one_rn = @reaction_network one_rn begin
    @parameters kcat1 Km1
    @species L(t) K(t)

    kcat1*(K)/(K*(Km1+L)), L + K --> Lp + K # L phosphorylation by kinase into Lp using Michaelis-Menten
end
test = convert(ODESystem, one_rn)
for i in test.eqs println(i) end
=#