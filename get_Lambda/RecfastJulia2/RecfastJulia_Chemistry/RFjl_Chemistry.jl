#### For now, not going to do anything with filtering the background. Might possibly speed things up, but I need to have a better picture of the full code before I do that ####


include("RFjl_ChemCnsts.jl")
include("RFjl_react.jl")
include("RFjl_cool.jl")

include("../src/recfast.jl")
using .Recfast


function set_aDM_vars(r_m, r_MH, r_alpha, xi, epsilon; Yp = 0.0, h100=0.6774, Omega_B = 0.044, Omega_M = 0.3089, Omega_K = 0.0)
    
    Omega_DM = Omega_M - Omega_B
    Omega_aDM = (epsilon/r_MH) * Omega_DM
    N_eff = 7.449 / (xi^4)
    T_DCMB_0 = xi*2.725
    z_rec_DM = ((r_m * (r_alpha^2))/xi) * (z_rec_sm + 1) - 1
    z_form = 100 #### Not correct, need to put PS treatement in 

    cps = Recfast.Params(Yp = Yp, T_0 = T_DCMB_0, Omega_M = Omega_M, Omega_B = Omega_aDM, Omega_K = Omega_K, h100 = h100, n_eff = N_eff)
    
    Set_VFC_params!(cps, mS= r_m, mxS = 14, aS = 1.37)

    return z_rec_DM, z_form, cps
end


function


function redshift_to_time_steps(z_interval, cps)
    # Constants for cosmological parameters (adjust as needed)

    H0 = 100*cps.h100
    Omega_m = cps.Omega_M
    Omega_lambda = cps.Omega_Lambda

    # Initialize an empty array to store time steps
    t_steps = zeros(z_interval)

    for z in z_interval
        # Define the scale factor at redshift z
        a_z = 1.0 / (1.0 + z)

        # Calculate the age of the universe at redshift z
        integrand(a) = 1.0 / (a * sqrt(Omega_m / a^3 + Omega_lambda))
        age_universe, _ = quadgk(integrand, a_z, 1.0)
        age_universe /= H0  # Convert to Gyr using Hubble constant

        push!(t_steps, age_universe)
    end

    return t_steps
end



# Example usage:
z_interval = log10.(range(1000, 999.9, length=10))  # Logarithmic redshift interval
t_steps = redshift_to_time_steps(z_interval)

function get_cooling()
	i=0
        list_new = [12-6, 1e-6, 1, 0, 0, 0]
        cool = 0
	n = 1
	T = 1e3
        z = 10
	r_m = r_MH = r_alpha = 1
	dt = 1e-4
	while i < 100
	   list = list_new
	   list_new = abundances(list,n,T,r_m, r_MH, r_alpha, dt)
           cool += total_cooling(n, T, z, r_m, r_MH, r_alpha, list_new[1], list_new[2], list_new[3], list_new[4], list_new[5])
	   i = i+1
	end
        println(cool)
end

get_cooling()
