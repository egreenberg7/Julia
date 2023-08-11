module GillespieConverter
export convertU0, convertP, getJumpProb, Nₐ

using Catalyst
using DifferentialEquations
using JumpProcesses

#Avogadro's number
const Nₐ = 6.02214076 * 1e23

"""
Converts array of concentrations to array of copy numbers by multiplying
by a given volume and Avogadro's constant and rounding to the nearest int.
- `u0` Vector of concentrations in μM
- `volume` Volume of reaction vessel in μm³
"""
function convertU0(u0, volume)
    return round.(Int, ((volume *  Nₐ * 10.0^-21) .* u0))
end

"""
Given a list of rate constants for fullrn, converts them to transition rates by
dividing by volume or multiplying by other factors as appropriate for the reaction order.
- `p` Vector of rate constants for the full reaction
- `volume` Volume of reaction vessel in μm³
"""
function convertP(p, volume)
    gillespieP = zeros(13)
    gillespieP[1] = p[1] / (volume *  Nₐ * 10.0^-21) #ka1, 2nd order
    gillespieP[2] = p[2] #kb1, 1st order
    gillespieP[3] = p[3] #kcat1, 1st order
    gillespieP[4] = p[4] / (volume *  Nₐ * 10.0^-21) #ka2, 2nd order
    gillespieP[5] = p[5] #kb2, 1st order
    gillespieP[6] = p[6] / (volume *  Nₐ * 10.0^-21) #ka3, 2nd order
    gillespieP[7] = p[7] #kb3, 1st order
    gillespieP[8] = p[8] / volume #ka4, 2nd order
    gillespieP[9] = p[9] #kb4, 1st order
    gillespieP[10] = p[10] / (volume *  Nₐ * 10.0^-21) #ka7, 2nd order
    gillespieP[11] = p[11] #kb7, 1st order
    gillespieP[12] = p[12] #kcat7, 1st order
    gillespieP[13] = p[13] #df, not dependent on order???
    return gillespieP
end

"""
Given a reaction, jumpU0, and jumpP, creates the jump problem for the reaction. Note that convertP and
convertU0 have to be called prior to using this.
- `rn` The reaction system
- `jumpU0` The initial copy numbers
- `tspan` The timespan to evaluate the numbers over.
- `jumpP` The transition rates
"""
function getJumpProb(rn, jumpU0, jumpP, tspan)
    dProb = DiscreteProblem(rn, jumpU0, tspan, jumpP)
    jProb = JumpProblem(rn, dProb, Direct()) 
    return jProb
end
end