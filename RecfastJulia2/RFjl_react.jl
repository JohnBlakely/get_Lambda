


using QuadGK

include("RFjl_ChemCnsts.jl")

#### k1 and k6 are the ones which use QuadGK ####

#### Reaction Network ####
#
# SPECIES
# p: Dark Proton
# e: Dark Electron
# H: Dark Atomic Hydrogen
# H2: Dark Molecular Hydrogen
# Everything else is just the aDM equivalent SM abundance
#
# ALL REACTIONS ARE FROM Ryan et al (2022) 
# NOTE: John check the references -John
#___________________________________________________
# NAME |          Reaction              | Reference
#------|--------------------------------|-----------
# k1   | p + e -> H                     | Rosenberg & Fan (2017) 
# k2   | H + e -> H^-                   | Galli & Palla (1998) 
# k3   | H^- + H -> H2 + e              | Kreckel et al (2010)
# k4   | H + p -> H2^+                  | Ramaker & Peek (1976) & Copolla et al (2011)
# k5   | H2^+ + H -> H2 + p             | Galli & Palla (1998)
# k6   | H + e -> p + 2e                | Rosenberg & Fan (2017)
# k7   | H2 + p -> H2^+ + H             | Savin et al (2004) & Glover et al (2010)
# k8   | H2 + H -> 3H                   | Martin et al (1996)
# k9   | H^- + p -> 2H                  | Stenrup et al (2009)
# k10  | H2^+ + e -> 2H                 | Copolla et al (2010)
# k11  | 3H -> H2 + H                   | Forrey (2013)
# k12  | 2H + H2 -> 2H2                 | Glover & Abel (2008)
# k13  | 2H + p -> H2 + p               | Yoshida (2006)
# k14  | 2H + p -> H2^+ + H             | Yoshida (2006)
# k15  | H2^+ + H2 -> H3^+ + H          | Galli & Palla (1998)
# k16  | H3^+ + e -> H + H2             | Galli & Palla (1998)
# _____|________________________________|___________
#
# m & MH are in grams
Th(T, r_m, r_alpha) = (T / (r_m * (r_alpha^2)))

function k1(T, r_m, r_MH, r_alpha)
    m = m_SM * r_m
    alpha = alpha_SM * r_alpha
    summation_interval = 10 .^ range(0, stop=100, length=25)
    argument(u, n) = (u * exp(-(u^2))) / (((u^2) * (n^3)) + ((atomBindE(T, r_m, r_alpha)^2)*n))
    summedArg(u) = sum(argument.(u, summation_interval))
    return ((HBAR^2) * (C^2)) * (alpha^5) * (( ((2^11) * pi) / ((3^3) * m * ((K_B * T)^3)) )^(1/2)) * quadgk(x -> summedArg(x), 0, Inf)[1]
end

k2(T, r_m, r_MH, r_alpha) = ((r_alpha/r_m)^2) * (1.4e-18) * ((T/ (r_m * (r_alpha^2)))^0.928) * exp(-T/ (16200 * r_m * (r_alpha^2)))


function k4(T, r_m, r_MH, r_alpha)
    Tcrit = 30 #(30 / (r_m * (r_alpha^2)))
    if T < Tcrit
        return ((r_alpha^2)/(r_m * r_MH)) * 2.1e-20 * ((((T/30)/ (r_m * (r_alpha^2))))^(-0.15))
    elseif T >= Tcrit
        return ((r_alpha^2)/(r_m * r_MH)) * 10^(-18.20 - 3.194*log10((T / (r_m * (r_alpha^2)))) + 1.786*(log10((T / (r_m * (r_alpha^2))))^2) - 0.2072*(log10((T / (r_m * (r_alpha^2))))^3))
    end
end

function k5(T, r_m, r_MH, r_alpha)
    return (6.0e-10) / (r_alpha * (r_m^1.5) * (r_MH^0.5))
end

function k3(T, r_m, r_MH, r_alpha)
        Th = T/(r_m * (r_alpha ^ 2))
        a1 = 1.35e-9
        a2 = 9.8493e-2
        a3 = 3.2852e-1
        a4 = 5.5610e-1
        a5 = 2.771e-7
        a6 = 2.1826
        a7 = 6.191e-3
        a8 = 1.0461
        a9 = 8.9712e-11
        a10 = 3.0424
        a11 = 3.2576e-14
        a12 = 3.7741
        return (((r_alpha) * (r_m ^ 1.5) * (r_MH ^ 0.5)) ^ -1) * a1 * ((Th ^ a2) + a3 * (Th ^ a4) + a5 * (Th ^ a6))/(1 + a7 * (Th ^ a8) + a9 * (Th ^ a10) + a11 * (Th ^ a12)) 
end

function k6(T, r_m, r_MH, r_alpha) 
    return ((HBAR ^ 2) * (C ^ 0)) * ((((2 ^ 7) * pi) / (((m_SM * r_m) ^ 3) * (K_B * T))) ^ (1/2)) * quadgk(x ->  ((1 - ((atomBindE(T, r_m, r_alpha)^2)/(x^2))) * (x * (exp(-(x^2))))), atomBindE(T, r_m, r_alpha), Inf)[1]
end

function k7(T, r_m, r_MH, r_alpha)
    Th = (T / (r_m * (r_alpha^2)))
    Tcrit1 = (100/T)*Th
    Tcrit2 = (30000/T)*Th
    lTh = log(Th)
    if T >= Tcrit1 && T<= Tcrit2
        return (((r_alpha^2)*(r_m^3)*(r_MH))^(-1/2)) * exp(-2.1237150e4/Th) * (-3.3232183e-7 + 3.3735382e-7*lTh - 1.4491368e-7*(lTh^2) + 3.4172805e-8*(lTh^3) - 4.7813728e-9*(lTh^4) + 3.9731542e-10*(lTh^5) - 1.8171411e-11*(lTh^6) + 3.5311932e-13*(lTh^7))
    else
        return 0
    end
end

function k8(n, T, r_m, r_MH, r_alpha)
    rrot = ((r_alpha*r_m)^2)/r_MH
    rvib = (((r_alpha^4)*(r_m^3))/(r_MH))^(1/2)
    Th = min(T / (r_m * (r_alpha^2)), 4.5e4)
    iTh = 1 / Th
    logTh = log10(Th)
    ilogTh = 1 / logTh
    logTh2 = logTh^2
    logTh3 = logTh^3
    Tv = Th * (r_m * (r_alpha^2)) * (rvib^-1)
    iTv = 1 / Tv
    logTv = log10(Tv)
    ilogTv = 1 / logTv
    logTv2 = logTv^2
    logTv3 = logTv^3

    coeff_CD = [-178.4239, -68.42243, 43.20243, -4.633167, 69.70086, 40870.38, -23705.70, 128.8953, -53.91334, 5.315517, -19.73427, 16780.95, -25786.11, 14.82123, -4.890915, 0.4749030, -133.8283, -1.164408, 0.8227443, 0.5864073, -2.056313]
    coeff_DT = [-142.7664, 42.70741, -2.027365, -0.2582097, 21.36094, 27535.31, -21467.79, 60.34928, -27.43096, 2.676150, -11.28215, 14254.55, -23125.20, 9.305564, -2.464009, 0.1985955, 743.0600, -1.174242, 0.7502286, 0.2358848, 2.937507]
    n_c_rescale = ((r_alpha^8) * (r_m^4.75) * (r_MH^-1.75))

    ## Collisional Dissociative ##
    CD_kh1 = coeff_CD[1] + coeff_CD[2] * logTh + coeff_CD[3] * logTh2 + coeff_CD[4] * logTh3 + coeff_CD[5] * log10(1 + coeff_CD[6] * iTh)
    CD_kh2 = coeff_CD[7] * iTh
    CD_kl1 = coeff_CD[8] + coeff_CD[9] * logTh + coeff_CD[10] * logTh2 + coeff_CD[11] * log10(1 + coeff_CD[12] * iTh)
    CD_kl2 = coeff_CD[13] * iTh

    CD_nc1 = coeff_CD[14] + coeff_CD[15] * logTv + coeff_CD[16] * logTv2 + coeff_CD[17] * iTv
    CD_nc2 = CD_nc1 + coeff_CD[18]

    p = coeff_CD[19] + coeff_CD[20] * exp(-Th / 1.85e3) + coeff_CD[21] * exp(-Th / 4.4e2)

    CD_nc_1 = (10 ^ CD_nc1 * n_c_rescale)
    CD_nc_2 = (10 ^ CD_nc2 * n_c_rescale)

    k_CD = CD_kh1 - ((CD_kh1 - CD_kl1) / (1 + ((n / CD_nc_1) ^ p))) + CD_kh2 - ((CD_kh2 - CD_kl2) / (1 + ((n / CD_nc_2) ^ p)))

    ## Dissociative Tunneling ##
    DT_kh1 = coeff_DT[1] + coeff_DT[2] * logTh + coeff_DT[3] * logTh2 + coeff_DT[4] * logTh3 + coeff_DT[5] * log10(1 + coeff_DT[6] * iTh)
    DT_kh2 = coeff_DT[7] * iTh
    DT_kl1 = coeff_DT[8] + coeff_DT[9] * logTh + coeff_DT[10] * logTh2 + coeff_DT[11] * log10(1 + coeff_DT[12] * iTh)
    DT_kl2 = coeff_DT[13] * iTh

    DT_nc1 = coeff_DT[14] + coeff_DT[15] * logTv + coeff_DT[16] * logTv2 + coeff_DT[17] * iTv
    DT_nc2 = DT_nc1 + coeff_DT[18]

    p = coeff_DT[19] + coeff_DT[20] * exp(-Th / 1.85e3) + coeff_DT[21] * exp(-Th / 4.4e2)

    DT_nc_1 = (10 ^ DT_nc1 * n_c_rescale)
    DT_nc_2 = (10 ^ DT_nc2 * n_c_rescale)

    k_DT = DT_kh1 - ((DT_kh1 - DT_kl1) / (1 + ((n / DT_nc_1) ^ p))) + DT_kh2 - ((DT_kh2 - DT_kl2) / (1 + ((n / DT_nc_2) ^ p)))

    k8 = (10^k_CD) + (10^k_DT)

    if isfinite(k8) == true && isnan(k8) != true
        return ((r_alpha^-1) * (r_m^-1.5) * (r_MH^-0.5)) * k8
    else
        return 0
    end
end

function k9(T, r_m, r_MH, r_alpha)
    if T <= (1e5 * r_m * (r_alpha^2)) && T <= (10 * r_m * (r_alpha^2))
        return ((r_alpha*r_m)^-3)*(((2.96e-6)/(Th^0.5)) - 1.73e-9 + (2.5e-10 * (Th^0.5)) - 7.77e-13 * Th)
    else
        return 0
    end
end


k10(T, r_m, r_MH, r_alpha) = 0

k11(T, r_m, r_MH, r_alpha) = ((r_alpha*r_m*(r_MH^0.25))^-4) * (6e-32*((T / (r_m * (r_alpha^2)))^-0.25) + (2e-31*((T / (r_m * (r_alpha^2)))^-0.5)))

k12(T, r_m, r_MH, r_alpha) = k11(T, r_m, r_MH, r_alpha) / 8

k13(T, r_m, r_MH, r_alpha) = k11(T, r_m, r_MH, r_alpha)

k14(T, r_m, r_MH, r_alpha) = k11(T, r_m, r_MH, r_alpha) / 8

k15(T, r_m, r_MH, r_alpha) = 2e-9 * (((r_alpha^2)*(r_m^3)*(r_MH))^(-1/2))

k16(T, r_m, r_MH, r_alpha) = 4.6e-6*((T / (r_m * (r_alpha^2)))^-0.65)*(((r_alpha^2)*(r_m^3)*(r_MH))^(-1/2))


function dx_dt(x,f,n,T,r_m, r_MH, r_alpha, xH, xHm, xH2p, xH3p)
    return n*((xHm*xH*k3(T, r_m, r_MH, r_alpha)) + (xH2p*xH*k5(T, r_m, r_MH, r_alpha)) + 2*(xH*x*k6(T, r_m, r_MH, r_alpha)) - (x*x*k1(T, r_m, r_MH, r_alpha)) - (x*xH*k2(T, r_m, r_MH, r_alpha)) - (x*xH*k4(T, r_m, r_MH, r_alpha)) - (f*x*k7(T, r_m, r_MH, r_alpha)) - (xHm*x*k9(T, r_m, r_MH, r_alpha)) - (xH2p*x*k10(T, r_m, r_MH, r_alpha)) - (n*xH*xH*x*k14(T, r_m, r_MH, r_alpha)) - (xH3p * x * k16(T, r_m, r_MH, r_alpha)))
end

function df_dt(x,f,n,T,r_m, r_MH, r_alpha, xH, xHm, xH2p, xH3p)
    return n*((xHm*xH*k3(T, r_m, r_MH, r_alpha)) + (f*xH*k5(T, r_m, r_MH, r_alpha)) + (n*xH*xH*f*(k12(T, r_m, r_MH, r_alpha) + k13(T, r_m, r_MH, r_alpha))) + (xH3p*x*k16(T, r_m, r_MH, r_alpha)) - (f*x*k7(T, r_m, r_MH, r_alpha)) - (f*xH*k8(n, T, r_m, r_MH, r_alpha)) - (xH2p*x*k16(T, r_m, r_MH, r_alpha)))

end

#### Where Tegmark 97 comes in ####

function dxH_dt(dx_dt, df_dt)
    return (-1 * dx_dt - 2 * df_dt)
end

function dOtherSpecies_dt(f, df_dt, xHm, xH2p, xH3p)
    xHm_frac = xHm/f
    xH2p_frac = xH2p/f
    xH3p_frac = 1-xH2p_frac - xHm_frac
    return xHm_frac * df_dt, xH2p_frac * df_dt, xH3p_frac * df_dt
end


#### Turn x,f,xH,xHm,xH2p,xH3p to a vector ####
function abundances(abunds,n,T,r_m, r_MH, r_alpha, dt)
    x,f,xH,xHm,xH2p,xH3p = abunds

    total = x + f + xH + xHm + xH2p + xH3p
    if isapprox(1, total, atol=a_tol, rtol=r_tol) == false
        ### Add error to errorlog in future ###
        diff = 1 - total
        
        ### Assuming that the problem is xH ###
        xH = xH + diff
    end
    

    dfdt = df_dt(x,f,n,T,r_m, r_MH, r_alpha, xH, xHm, xH2p, xH3p)
    df = dfdt * dt
    dxdt = dx_dt(x,f,n,T,r_m, r_MH, r_alpha, xH, xHm, xH2p, xH3p)
    dx = dxdt * dt
    dxHdt = dxH_dt(dxdt, dfdt)
    dxH = dxHdt * dt
    dxHm, dxH2p, dxH3p = dOtherSpecies_dt(f, dfdt, xHm, xH2p, xH3p) .* dt
   
    x_new = x + dx
    f_new = f + df
    xH_new = xH + dxH
    xHm_new = xHm + dxHm
    xH2p_new = xH2p + dxH2p
    xH3p_new = xH3p + dxH3p
    total_new = x_new + f_new + xH_new + xHm_new + xH2p_new + xH3p_new

    if isapprox(1, total_new, atol=a_tol, rtol=r_tol) == false
        ### Add error to errorlog in future ###
        diff = 1 - total
        
        ### Assuming that the problem is xH_new ###
        xH_new = xH_new + diff
    end
    
    
    return x_new, f_new, xH_new, xHm_new, xH2p_new, xH3p_new
end
    


    
    
    
    
    



