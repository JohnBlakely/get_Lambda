using PyPlot
include("src/recfast.jl")
using .Recfast
p = Recfast.Params(Yp = 0.0,  T0 = 0.0546, Omega_M =0.3, Omega_B = 0.0156, h100 = 0.7, n_eff = 4.66e7)
Recfast.Set_VFC_params!(p, aS = 1.37*5, mS = .274, mxS = 16)
function testh2dbh(dtmax, zend)
    a = Recfast.Evaluate_recombination_h2(p, zstart = 4.64e6, dtmax = dtmax, zend = zend)
    tpoints = a.t[1:1000:end]
    #plot(tpoints, [g[1] for g in a(tpoints)], label = L"x_{He}")

    loglog(tpoints, [g[2] for g in a(tpoints)], label = L"x_H^+")

    #loglog(tpoints, [g[3] for g in a(tpoints)])

    loglog(tpoints, [g[4] for g in a(tpoints)], label = L"x_{H^-}")

    loglog(tpoints, [g[5] for g in a(tpoints)], label = L"x_{H_2^{+}}")

    loglog(tpoints, [g[6] for g in a(tpoints)], label = L"x_{H_2}")
    #scatter([250, 78], fittingh2()[2:2:4], label = "fitting formula")
    #loglog(a.t[1:100:end], [u[6] for u in a.u[1:100:end]])
    ##loglog(a.t[1:100:end], [u[4] for u in a.u[1:100:end]])
    #loglog(a.t[1:100:end], [u[5] for u in a.u[1:100:end]])
    xlim(1e5, 1)
    ylim(1e-25, 10)
    xlabel("z")
    legend()
    #savefig("h2sm.png")
    return a
end

#a fitting formula for the abundance from the H+ channel
#function fittingh2()
#    a = Recfast.Evaluate_recombination_h2(p, zstart = 1e4, dtmax = 100)
    #This is 2.30 in First Light to Reionization: x_H2 = x_HP R n_HI(zeff) t_H(zeff)
    #motivated by laziness, we just take the values from table 2.1 where possible
#    return (a(250)[end], 2.08e-18*1.54e14 *3.15 * a(250)[2],#
#    a(78)[end] - a(250)[end], 1.08e-16 * 8.49e14 * 0.098 * a(78)[2])

#end
