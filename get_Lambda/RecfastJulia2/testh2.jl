using PyPlot
include("src/recfast.jl")
using .Recfast
p = Recfast.Params()
function testh2(dtmax)

    p = Recfast.Params()
    a = Recfast.Evaluate_recombination_h2(p, logzstart = 4., dtmax = dtmax)
    tpoints = a.t
    fig, ax = subplots()
    plot(tpoints, [g[1] for g in a.(tpoints)], label = L"x_{He}")

    semilogy(tpoints, [g[2] for g in a.(tpoints)], label = L"x_H")

    #loglog(tpoints, [g[3] for g in a(tpoints)])

    semilogy(tpoints, [g[4] for g in a.(tpoints)], label = L"x_{H^-}")

    semilogy(tpoints, [g[5] for g in a.(tpoints)], label = L"x_{H_2^{+}}")

    semilogy(tpoints, [g[6] for g in a.(tpoints)], label = L"x_{H_2}")
    scatter([log10(250), log10(78)], fittingh2()[2:2:4], label = "fitting formula")



    #loglog(a.t[1:100:end], [u[6] for u in a.u[1:100:end]])
    ##loglog(a.t[1:100:end], [u[4] for u in a.u[1:100:end]])
    #loglog(a.t[1:100:end], [u[5] for u in a.u[1:100:end]])
    xlim(4, 1)
    ylim(1e-25, 10)
    xlabel("z")
    legend()
    savefig("h2sm.png")
    return a
end

#a fitting formula for the abundance from the H+ channel
function fittingh2()
    a = Recfast.Evaluate_recombination_h2(p, logzstart = 4., dtmax = .001)
    #This is 2.30 in First Light to Reionization: x_H2 = x_HP R n_HI(zeff) t_H(zeff)
    #motivated by laziness, we just take the values from table 2.1 where possible
    return (a((250))[end], 2.08e-18*1.54e14 *3.15 * a((250))[2],
    a((78))[end] - a((250))[end], 1.08e-16 * 8.49e14 * 0.098 * a((78))[2])

end
