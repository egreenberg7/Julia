#! Solve model for arbitrary oscillatory parameters and initial conditions
function make_ODEProb(rn = make_fullrn(); psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
                                                        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
                                                        :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606], 
                                        usym = [:L => 0.0, :K => 0.5, :P => 0.3, :A => 2.0, :Lp => 3.0, :LpA => 0.0, :LK => 0.0, 
                                                :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0],
                                        tspan = (0., 100.) 
                                                )
    #? Parameter list
    p= [x[2] for x in psym]
 
    #? Initial condition list
    u0 = [x[2] for x in usym]

    #? Create ODE problem
    return ODEProblem(rn, u0, tspan, p) #TODO fix type stability
end



function make_nerdssODEProb(rn; p::Vector{Float64} = [1.4353350626245021e-5, 0.20705634097197367, 461.9447701768986, 0.6642156268695386, 15.902417897116093, 7.971443130885463e-5, 0.7016715934086192, 
                                                1.5736952400326606e-5, 0.10411397662131976, 0.0007148707925785353, 0.05915700148314913, 0.8831745113213687, 1500.0], 
                                        u0::Vector{Float64} = [0.0, 301.15000000000003, 240.92000000000004, 602.3000000000001, 903.45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                        tspan = (0., 100.) 
                                                )
    #? Create ODE problem
    return ODEProblem(rn, u0, tspan, p) #TODO fix type stability
end
# const fullprob = make_ODEProb()