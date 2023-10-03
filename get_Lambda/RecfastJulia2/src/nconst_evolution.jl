####We want to use the Recfast machinery to do non cosmology evolutions

# Collisional ionization
function sigmaciv(m, alpha, T)
   if T < 0
    return 0
end
    T = T * 8.617333262145e-5 # convert ev
    y = sqrt(alpha^2 * m/(2*T))
    integrand(u) = u * exp(-u^2)/(1+ (2 * y^2/u^2)) * ( 1 - y^2/u^2 - 0.5 * (1 - y^4/u^4) * log(y^2/u^2)
        + y^2 * log(y^2/u^2)/(u^2+y^2))
    #uevals = 10 .^ range(log10(y), stop = log10(1000 * y) , length = 50)
    #return sqrt(2^7 * pi/(m^3 * T)) * sum(integrand.(uevals[1:end-1]) .* diff(uevals))
    RF_hb_ev^2 * RF_cLight^3* sqrt(2^7 * pi/(m^3 * T)) * quadgk(u -> integrand(u), y, Inf )[1]
end

#And the corresponding, Case A reocmbination rate
function sigmarecv(m, alpha, T)
    if T < 0
        return 0
    end
    T = T * 8.617333262145e-5 # convert ev
    y = sqrt(alpha^2 * m/(2*T))
    summand(u, n) = u * exp(-u^2)/(u^2 * n^3 + y^2 * n)
    integrand(u) = sum(summand.(u, 1:20))
    return RF_hb_ev^2 * RF_cLight^3*alpha^5 * sqrt(2^11 * pi/(3^3 * m * T^3)) * quadgk(u -> integrand(u), 0, Inf)[1]
end


function Evaluate_chemical_network_nconst(p, nHTot, t, y)
    cosmoVals = calc_abundances_nconst(p, nHTot, y)
    calc_dxe_nconst(p, t, cosmoVals)
end



function calc_abundances_nconst(p, nHTot, y)
    xp = y[2]
    TM = y[3]
    x_Hm = y[4]
    x_H2p = y[5]
    x_H2 = y[6]

    x_H3p = y[7]
    xe = max(0, xp -  x_Hm + x_H2p + x_H3p)
    xH1 = max(0,(1 - xp - x_Hm - 2 * (x_H2p + x_H2) - 3 * x_H3p))

    nH1 =   xH1 * nHTot
    np = xp * nHTot
    nH2 = x_H2 * nHTot
    nH3p = x_H3p * nHTot
    return @SVector [nHTot, xp, TM, x_Hm, x_H2p, x_H2,
    xe, xH1, nH1, np, nH2,  x_H3p,  nH3p]
end

function calc_dxe_nconst(p, t, cosmoVals)
    nHTot,  xp, TM, x_Hm, x_H2p, x_H2, xe, xH1, nH1, np, nH2, x_H3p, nH3p = cosmoVals



    #volumetric rate, J s^-1 m^-3

    cooling_rate = (atomic_cooling(p, TM, xe, xH1, nHTot)
        + ionization_cooling(p, TM, xe, xH1, nHTot)
        + ff_cooling(p, TM, xe, xp, nHTot)
        + hh2cooling(p, TM, nH1, nH2)
        + h2h2cooling(p, TM, nH2))
    nTot = (xH1 + xp + xe + x_H2p + x_Hm + x_H2 + x_H3p) * ( nHTot)
    #println("cooling", cooling_rate)
    #println(atomic_cooling(p, TM, xe, xH1, nHTot))

   # println(t, [xp, TM, x_Hm, x_H2p, x_H2, xH1])
    gamma = (5 * (xe + xH1) + 7 * (x_H2))/(3 * (xe + xH1) + 5 * (x_H2))
    dTM = -(cooling_rate/RF_kBoltz) * (gamma-1)/nTot


    dxp  =  (sigmaciv(RF_me_eV* p.rescaling.RF_mScale, RF_alpha* p.rescaling.RF_aScale, TM) * nH1 * xe
            - sigmarecv(RF_me_eV* p.rescaling.RF_mScale, RF_alpha * p.rescaling.RF_aScale, TM) * xe * xp * nHTot
            -ch7(p,TM)* np * x_Hm -ch8(p, TM) * xp * nH1+  ch10(p) * x_H2p * nH1 -  ch15(p, TM) * x_H2 * np )


    dHm = (ch3(p, TM) * xe * nH1 - ch7(p, TM) * np * x_Hm
                -  ch5(p, TM) * x_Hm * nH1 )

    dH2p =  (ch8(p, TM) * xp * nH1  - ch10(p) * x_H2p * nH1 -ch13(p, TM) * nH2 * x_H2p
                + ch15(p, TM) * x_H2 * np )


    dH2 =  ( ch5(p, TM) * x_Hm * nH1  + ch10(p) * x_H2p *nH1 -ch13(p, TM) * nH2 * x_H2p
            -  ch15(p, TM) * x_H2 * np  + ch20(p,TM)*nH3p *xe- h2hdis(p, nH1, TM) *x_H2 * nH1) #- ch17(p, TM) * nH2 * xe)

    dH3p = ( ch13(p, TM) * nH2 * x_H2p - ch20(p,TM)*nH3p *xe )

   # println([0, dxp, dTM, dHm, dH2p, dH2])
  # println("[0, dxp, dTM, dHm, dH2p, dH2]",[0, dxp, dTM, dHm, dH2p, dH2])
    return @SVector [0, dxp, dTM, dHm, dH2p, dH2, dH3p]
end
"""
$(TYPEDSIGNATURES)
Evolve the chemical network from a starting time to an ending time (in seconds),
at a fixed specified density.
"""
function Evaluate_abundances_nconst(p, nHTot, initial_conditions, tstart, tend; dt = 1e-5,
    reltol = 1e-10, abstol = 1e-8, dtmax = 1e12, dtmin = 1e-10, maxiters=1e3)
    y0 = initial_conditions
    function evalODE(y, k, t)
        dy = Evaluate_chemical_network_nconst(p, nHTot, t, y)
    end
    trange = (tstart, tend)
    prob = ODEProblem(evalODE, y0, trange)
    sol = solve(prob, Rodas4P2(autodiff=false), isoutofdomain=(u,p,t) -> any(x -> x < 0, u), dt = dt, dtmax = dtmax, reltol = reltol, abstol = abstol, dtmin = dtmin, maxiters=maxiters)

    return sol

end
