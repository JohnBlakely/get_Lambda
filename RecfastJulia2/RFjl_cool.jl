


include("RFjl_react.jl")
include("RFjl_ChemCnsts.jl")

using SpecialFunctions
using QuadGK













function atomic_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH)


    m = m_SM * r_m
    MH = MH_SM * r_MH
    alpha = alpha_SM * r_alpha
    T_DarkCMB = 2.73 * xi * (1 + z) 

    y = ((m * ((C * alpha)^2))/(2 * K_B * T))^(1/2)
    y2 = y^2
    log10y2 = log10(y2)
    lny2 = log(y2)

    prefac = ((HBAR^2) * (C^3))

    T5 = T/1e5
    Th = T/(r_m * (r_alpha^ 2))
   
    xH = (1 - x_ion - 2*f)

    atomBindE(T,r_m, r_alpha) = (((m_SM * r_m * (C^2)) * ((r_alpha * alpha_SM)^2)) / (2 * K_B * T))^(1 / 2)

    rrot = ((r_alpha*r_m)^2)/r_MH
    rvib = (((r_alpha^4)*(r_m^3))/(r_MH))^(1/2)

    #### Compton ####
    cmptn = ((K_B^5)/(HBAR*(C^6))) * (4 * (T-T_DarkCMB)/m) * ((8*pi)/3) * ((alpha/m)^2) * ((pi^2)/15) * (T_DarkCMB^4)

    #### Collisional Ionization #### (Using Binary-Encounter)
    f_approx = (exp(-y2)/2) + ((y2*expinti(-y2))/2)
    collion = prefac * ((((2^3)*pi*(alpha^4))/(m*K_B*T))^(1/2)) * f_approx

    #### Recombination ####
    bigy = prefac * (((2^9)*pi*(alpha^6)*K_B*T)/((3^3)*(C^6)*(m^3))) * (0.74 + log10y2 + (1/(3*y2)))
    lily = prefac * (((2^5)*pi*(alpha^10))/((3^3)*m*T*(C^2)*K_B)) * (5 + y2 * (2.860 + (14*log10y2)/3))
    if y>0.25
        recomb = bigy
    else
        recomb = lily
    end

    
    #### Collisional Excitation ####
    gapprox_arg(u) = ((log10(4*u/y) * u * exp(-1 * u^2))/(1 + (7*y2/(4 * (u^2)))))
    gapprox = quadgk(x -> gapprox_arg(x), ((0.75*y2)^0.5), Inf)[1]
    collexci = prefac * ((((2^16)/(3^9)) * ((2 * pi * (alpha^4))/(m * (C^2) * K_B * T)))^0.5) * gapprox


    #### Bremsstrahlung ####
    brem = prefac * (((2^9)*pi*(alpha^6)*K_B*T*(g_ff^2))/((3*m*(C^2))^3))


    #### Total #### (erg * cm^3 * s^-1)
    return x_ion*((x_ion)*(brem + recomb) + (x_ion) * (collion + collexci) + cmptn)

end






function molecular_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH)



#### Low density ####

    TscaleR = (r_alpha ^ -2) * (r_m ^ -2) * (r_MH)
    TscaleV = (r_alpha ^ -2) * (r_m ^ -1.5) * (r_MH ^ 0.5)
    TscaleA = (r_alpha ^ -2) * (r_m ^ -1)

    T_r = TscaleR * T
    T_v = TscaleV * T
    rrot = ((r_alpha*r_m)^2)/r_MH
    rvib = (((r_alpha^4)*(r_m^3))/(r_MH))^(1/2)

    

    function H_H2_coll(T, winL, winH)
        if T <= T_H_lim_1
            coefficient = a_i_H_10_100
        elseif T > T_H_lim_1 && T <= T_H_lim_2
            coefficient = a_i_H_100_1000
        elseif T > T_H_lim_2 && T <= T_H_lim_3
            coefficient = a_i_H_1000_6000
        end
        if @isdefined(coefficient) == true
            return 10^fitting_function(coefficient, T, 0)
        else
            return 1.862314467912518e-22 * wCool(log10(T), 1 + winL, log10(6000) + winH)
        end
    end
   
    function H2_H2_coll(T, winL, winH)
        H2_H2 = (10^fitting_function(a_i_h2_100_6000, T, 0)) * wCool(log10(T), 2 + winL, 4 + winH)
        if isfinite(H2_H2) == true
            return H2_H2
        else
            return 0
        end
    end
    
    function e_H2_coll(T, winL, winH)
        if T <= T_e_lim_1
            coefficient = a_i_e_10_200
        else
            coefficient = a_i_e_200_10000
        end
        return (10 ^ fitting_function(coefficient, T, 1)) * wCool(log10(T), 2 + winL, 4 + winH)
    end
  
    function p_H2_coll(T, winL, winH)
        if T<=1e4
            return 10^fitting_function(a_i_p_10_10000, T, 0)
        else
            return (1.182509139382060e-21) * wCool(log10(T), 1 + winL, 4 + winH)
        end
    end




    function H2_H2_DM_coll(T)
        T_0_SM = 5.40244e3
        T_0_DM = (((r_alpha) ^ 2) * ((r_m) ^ 1.5882) * ((r_MH) ^ -0.5861)) * T_0_SM
        T_0_r  = T_0_DM / TscaleR
        T_0_v  = T_0_DM / TscaleV

  
        H2_H2_DM_coll_r = H2_H2_coll(T_r, 0, log10((TscaleR/TscaleA)))
        H2_H2_DM_coll_v = H2_H2_coll(T_v, log10((TscaleV/TscaleR)), log10((TscaleV/TscaleA)))
   
        if r_MH <= r_m
            if T <= T_0_r && T < T_0_v
                return ((r_alpha*r_m)/(r_MH))*H2_H2_DM_coll_r
            elseif T >= T_0_v && T > T_0_r
                return (((r_alpha*r_m)/(r_MH^3))^0.25)*H2_H2_DM_coll_v
            end
        else 
            if T <= T_0_r
                return ((r_alpha*r_m)/(r_MH^2))*H2_H2_DM_coll_r
            elseif T > T_0_r && T <= T_0_DM
                x2 = T_0_r
                x1 = 0.98*x2
                y2 = H2_H2_coll(x2, 0,0)
                y1 = H2_H2_coll(x1, 0,0)
                slope = (y2 - y1)/(x2 - x1)
                return ((r_alpha*r_m)/(r_MH^2)) * (slope * (T_r - T_0_r) + y2)
            elseif T > T_0_DM && T <= T_0_v
                x2 = T_0_v
                x1 = 1.01*x2
                y2 = H2_H2_coll(x2, 0,0)
                y1 = H2_H2_coll(x1, 0,0)
                slope = (y2 - y1)/(x2 - x1)
                return (((r_alpha)*(r_m^0.25))/(r_MH^1.25)) * (slope * (T_v - T_0_v) + y2)
            else T > T_0_v
                return (((r_alpha)*(r_m^0.25))/(r_MH^1.25)) * H2_H2_DM_coll_v
            end
        end
    end

    function H_H2_DM_coll(T)
        T_0_SM = 855.833
        T_0_DM = ((r_alpha ^ -2) * ((r_m) ^ 1.5) * ((r_MH) ^ -0.5)) * T_0_SM
        T_0_r  = T_0_DM / TscaleR
        T_0_v  = T_0_DM / TscaleV


        H_H2_DM_coll_r = H_H2_coll(T_r, 0, log10((TscaleR/TscaleA)))
        H_H2_DM_coll_v = H_H2_coll(T_v, log10((TscaleV/TscaleR)), log10((TscaleV/TscaleA)))

        if r_MH <= r_m
            if T <= T_0_r && T < T_0_v
                return ((r_alpha*r_m)/(r_MH))*H_H2_DM_coll_r
            elseif T >= T_0_v && T > T_0_r
                return (((r_alpha*r_m)/(r_MH^3))^0.25)*H_H2_DM_coll_v
            end
        else 
            if T <= T_0_r
                return ((r_alpha*r_m)/(r_MH^2))*H_H2_DM_coll_r
            elseif T > T_0_r && T <= T_0_DM
                x2 = T_0_r
                x1 = 0.98*x2
                y2 = H_H2_coll(x2, 0,0)
                y1 = H_H2_coll(x1, 0,0)
                slope = (y2 - y1)/(x2 - x1)
                return ((r_alpha*r_m)/(r_MH^2)) * (slope * (T_r - T_0_r) + y2)
            elseif T > T_0_DM && T <= T_0_v
                x2 = T_0_v
                x1 = 1.01*x2
                y2 = H_H2_coll(x2, 0,0)
                y1 = H_H2_coll(x1, 0,0)
                slope = (y2 - y1)/(x2 - x1)
                return (((r_alpha)*(r_m^0.25))/(r_MH^1.25)) * (slope * (T_v - T_0_v) + y2)
            else T > T_0_v
                return (((r_alpha)*(r_m^0.25))/(r_MH^1.25)) * H_H2_DM_coll_v
            end
        end
    end

    winL = 0
    winH = log10((TscaleR/TscaleA))

    p_H2_DM_coll(T) = ((((r_alpha^2)*(r_m))/(r_MH^3))^(1/2)) * p_H2_coll(T_r, winL, winH)
    e_H2_DM_coll(T) = (r_alpha/r_MH) * e_H2_coll(T_r, winL, winH)


    total_low_density_H2(T) = n * ( (f * H2_H2_DM_coll(T)) + (xH * H_H2_DM_coll(T)) + (x_ion * (e_H2_DM_coll(T) + p_H2_DM_coll(T))))

    T3 = T/1e3

    R_SM(T) = ((9.5e-22 * (T3^3.76))/(1 + 0.12 * (T3^2.1))) * (exp((-0.13/T3)^3)) + (3e-24*exp(-0.51/T3))
    R_DM(T) = (((r_m^8)*(r_alpha^9))/(r_MH^6)) * R_SM(TscaleR * T)

    V_SM(T) = 6.7e-19 * exp(-5.86/T3) + 1.6e-18 * exp(-11.7/T3)
    V_DM(T) = (((r_m^5)*(r_alpha^9))/(r_MH^3)) * V_SM(TscaleV * T)

    total_high_density_H2(T) = R_DM(T) + V_DM(T)

    if total_high_density_H2(T) == 0 || total_low_density_H2(T) == 0
        return 0
    else 
        return (n*f)/((1/total_high_density_H2(T)) + (1/total_low_density_H2(T)))
    end
end


#chemical_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH, xHm) = (((r_alpha^2)/r_m)*ev_to_erg) * ( 4.48 * f * xH * k8(n, T, r_m, r_MH, r_alpha) + xHm * xH * k4(T, r_m, r_MH, r_alpha))



total_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH, xHm) = atomic_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH) + molecular_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH) # + chemical_cooling(n, T, z, r_m, r_MH, r_alpha, xi, x_ion, f, xH, xHm)





































    

