#### For now, not going to do anything with filtering the background. Might possibly speed things up, but I need to have a better picture of the full code before I do that ####


include("RFjl_ChemCnsts.jl")
#### For now, not going to do anything with filtering the background. Might possibly speed things up, but I need to have a better picture of the full code before I do that ####


include("RFjl_ChemCnsts.jl")
include("RFjl_react.jl")
include("RFjl_cool.jl")

include("src/recfast.jl")
using .Recfast


function set_aDM_vars(r_m, r_MH, r_alpha, xi, epsilon; Yp = 0.0, h100=0.6774, Omega_B = 0.044, Omega_M = 0.3089, Omega_K = 0.0)
    
    Omega_DM = Omega_M - Omega_B
    Omega_aDM = (epsilon/r_MH) * Omega_DM
    N_eff = 7.449 / (xi^4)
    T_DCMB_0 = xi*2.725
    z_rec_DM = ((r_m * (r_alpha^2))/xi) * (1100 + 1) - 1
    z_form = 100 #### Not correct, need to put PS treatement in 

    cps = Recfast.Params(Yp = Yp, T0 = T_DCMB_0, Omega_M = Omega_M, Omega_B = Omega_aDM, Omega_K = Omega_K, h100 = h100, n_eff = N_eff)
    
    Set_VFC_params!(cps, mS= r_m, mxS = 14, aS = 1.37)

    return z_rec_DM, z_form, cps
end


### Stolen from ChatGPT
function redshift_to_time_steps(z_interval, cps)
    # Constants for cosmological parameters (adjust as needed)

    H0 = 100*cps.h100
    Omega_m = cps.Omega_M
    Omega_lambda = cps.Omega_L

    # Initialize an empty array to store time steps
    t_steps = zeros(length(z_interval))

    for z in z_interval
        # Define the scale factor at redshift z
        a_z = 1.0 / (1.0 + z)

        # Calculate the age of the universe at redshift z
        integrand(a) = 1.0 / (a * sqrt(Omega_m / a^3 + Omega_lambda))
        age_universe, _ = quadgk(integrand, a_z, 1.0)
        age_universe /= H0  # Convert to Gyr using Hubble constant

        push!(t_steps, age_universe)
    end

    return zip(t_steps,z_interval)
end


function get_cooling(n, T, z_t_list, r_m, r_MH, r_alpha, input_abunds::Vector)
	i=0
        list_new = [12-6, 1e-6, 1, 0, 0, 0]
        cool = 0
	n = 1
	T = 1e3
	r_m = r_MH = r_alpha = 1
	dt = 1e-4
	for i in z_t_list
	   dt, z = i
	   list = list_new
	   list_new = abundances(list,n,T,r_m, r_MH, r_alpha, dt)
           cool += total_cooling(n, T, z, r_m, r_MH, r_alpha, list_new[1], list_new[2], list_new[3], list_new[4], list_new[5])
	end
        return cool
end

function run_cooling(n, T, z::Vector, r_m, r_MH, r_alpha, xi, epsilon; Yp = 0.0, h100=0.6774, Omega_B = 0.044, Omega_M = 0.3089, Omega_K = 0.0)
    z_rec_DM, z_form, cps = set_aDM_vars(r_m, r_MH, r_alpha, xi, epsilon; Yp = 0.0, h100=0.6774, Omega_B = 0.044, Omega_M = 0.3089, Omega_K = 0.0)
    zt_list = redshift_to_time_steps(z, cps)
    sol = Evaluate_recombination_h2(cps, logzstart = log10(z_rec_DM), logzend = log10(z_form), dt = 1e-10, dtmin = 1e-30, dtmax= 0.01)
    sol_last = sol.u[end]
    abund_list = [sol_last[2], sol_last[6], 1-sol_last[2]-2*sol_last[6], sol_last[4], sol_last[5], 0.0] # 0 for H3p for now
    coolin = get_cooling(n, T, zt_list, r_m, r_MH, r_alpha, abund_list)
    return coolin
end

z_test = 10 .^ range(2, stop=0, length=50)
## Currently there's a bug where you can start at a time before where the initial abunds are calculated

println(run_cooling(1, 1e3, z_test, 1, 1, 1, 1, 1))
