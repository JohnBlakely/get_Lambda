
#=
function Evaluate_Recfast_System_h2(p, logz, y)
    #calc h2 here
    cosmoVals = calc_cosmoVals_h2(p, 10^logz, y)
    calc_dxe_h2(p, logz, cosmoVals)

end

function calc_cosmoVals_h2(p, z, y)
    nHTot = NH(p, z)
    fHe = p.fHe
    xHep = y[1]
    xp = y[2]
    TM = y[3]
    x_Hm = y[4]
    x_H2p = y[5]
    x_H2 = y[6]

    TR = TCMB(p, z)
    xe = max(0, xp + xHep -  x_Hm + x_H2p)
    xH1 = max(0,(1 - xp - x_Hm - 2 * (x_H2p + x_H2)))
    nH1 =   xH1 * nHTot
    np = xp * nHTot
    nH2 = x_H2 * nHTot
    nHe1 = nHTot * fHe - nHTot * xHep
    Hz = H_z_loc(p, z)
    nTot = (xH1 + xp + xe + x_H2p + x_Hm + x_H2 + fHe) * ( nHTot)

    return [nTot, nHTot, xHep, xp, TM, x_Hm, x_H2p, x_H2,  TR,
    xe, xH1, nH1, np, nH2, nHe1, Hz]
end
function calc_dxe_h2(p, logz, cosmoVals)
    nTot, nHTot, xHep, xp, TM, x_Hm, x_H2p, x_H2,  TR,
    xe, xH1, nH1, np, nH2, nHe1, Hz = cosmoVals
    z = 10^logz

            dHm = -(z * log(10)) * (ch3(p, TM) * xe * nH1 -  ch4(p, TR) * x_Hm - ch7(p, TM) * np * x_Hm
                        -  ch5(p, TM) * x_Hm * nH1) / (Hz * (1 + z))
            dH2p = - (z * log(10)) * (ch8(p, TM) * xp * nH1 - ch9(p, TR) * x_H2p -  ch10(p) * x_H2p * nH1
                        + ch15(p, TM) * x_H2 * np + ch18(p, TR) *x_H2 + h_h_hp_to_h2p(p, TM) * nH1^2 * xp) / (Hz * (1 + z))
            dH2 = -  (z * log(10)) * ( ch5(p, TM) * x_Hm * nH1  + ch10(p) * x_H2p *nH1
                    -  ch15(p, TM) * x_H2 * np -ch18(p, TR) *x_H2 - h2hdis(p, nH1, TM) *x_H2 * nH1
                    + threeh_association(p, TM) * nH1^2 * xH1 + h_h_h2_association(p, TM) * nH1^2 * x_H2 + h_h_hp_to_h2(p,TM) * nH1^2 * xp) / ((1 + z) * Hz)

            #compute the adjustment to the recombination rate from these reactions
            dHp =  -  (z * log(10)) * (-ch7(p,TM)* np * x_Hm -ch8(p, TM) * xp * nH1 + ch9(p, TR) * x_H2p +  ch10(p) * x_H2p * nH1 -  ch15(p, TM) * x_H2 * np
                        - h_h_hp_to_h2(p,TM) * nH1^2 * xp)/((1 + z) * Hz)
            recombination = calc_dxe(p, logz, [nTot, nHTot, xHep, xp, TM, TR, xe, xH1, nH1, nHe1, Hz])

    return @SVector [recombination[1], recombination[2] + dHp, recombination[3], dHm, dH2p, dH2]

end

"""
$(TYPEDSIGNATURES)
Evaluate the recombination history, including molecular hydrogen. Note that we assume everything is completely
ionized at `zstart`. If this is not the case, we will get wrong results. Returns
a DifferentialEquations.jl solution object. If you have problems, try adjusting
the solver in the source code.
"""
function Evaluate_recombination_h2(p; logzstart = 5., logzend = 1, dt = 1e-5,
    reltol = 1e-10, abstol = 1e-8, dtmax = 0.5, dtmin = 1e-20)
    fHe = p.fHe
    Xe_He0 = fHe
    Xe_H0 = 1
    TM0 = p.T0 * (1 + 10^logzstart)


    X_Hm = 0
    X_H2p = 0
    X_H2 = 0
    y0 = @SVector [Xe_He0, Xe_H0, TM0, X_Hm, X_H2p, X_H2]
    function evalODE(y, k, logz)
        dy = Evaluate_Recfast_System_h2(p, logz, y)

    end
    zrange = (logzstart, logzend)
    prob = ODEProblem(evalODE, y0, zrange)
    cf(integrator) = choice_function(integrator, p)
    alg_switch = CompositeAlgorithm((Rodas4P2(), Rodas5(autodiff=false)), cf)
    sol = solve(prob, alg_switch, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), dt = dt, dtmax = dtmax, reltol = reltol, abstol = abstol, dtmin = dtmin, maxiters=1e6)

    t =  10 .^sol.t
    u = sol.u

    function mySol(z)
        t
        u
        if z < t[1]
            return sol(log10(z))
        else
            return [u[1][1], u[1][2], u[1][3] * (1 +z)/(1 + t[1]), u[1][4], u[1][5], u[1][6]]
        end
    end


    return mySol

end
