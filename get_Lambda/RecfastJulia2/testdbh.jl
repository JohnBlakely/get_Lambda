include("src/recfast.jl")
#include("calc_visibility.jl")
using .Recfast

 p = Recfast.Params(Yp = 0.0,  T0 = 0.0546, Omega_M =0.3, Omega_B = 0.0156, h100 = 0.7, n_eff = 4.66e7)
 Recfast.Set_VFC_params!(p, aS = 1.37, mS = 0.274)
 sol = Recfast.Evaluate_recombination(p, zstart = 2.5e5, dtmax = 5.)
 az, aXe,aTM = read_Recfast_data("Recfast++.vDJ\\outputs\\Xe_Recfast.alp.me.all.dat")
 loglog(sol.t, [u[2] for u in sol.u], "k-", label = "julia code")
 loglog(az, aXe, "r--", label = "c code")
 xlabel("z")
 ylabel(L"x_e(z)")
 maximum(([u[2] for u in sol(az).u] .- aXe) ./(aXe))
