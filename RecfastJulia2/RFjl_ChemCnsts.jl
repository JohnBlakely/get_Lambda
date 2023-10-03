const HBAR = 1.054571817e-27  # erg s
const H = 6.62607015e-27      # erg s
const C = 2.9979245800e10     # cm/s
const G_N = 6.6743e-8         # cm^3/g s^2
const ev_to_erg = 1.602176634e-12  
const K_B = 1.380649e-16      # erg/K
const K_B_EV = 8.617333262e-5
const m_SM = 9.109e-28        # g
const MH_SM = 1.673e-24        # g
const alpha_SM = 7.2973525693e-3  #
const m_SM_eV = 5.11e8
const g_ff = 1.5
const r_tol = 1e-3
const a_tol = 1e-6

atomBindE(T,r_m, r_alpha) = (((m_SM * r_m * (C^2)) * ((r_alpha * alpha_SM)^2)) / (2 * K_B * T))^(1 / 2)


## Coefficients ##
## From Table 8 in Glover and Abel 2008 and Appendix A of Glover 2015##
a_i_H_10_100 = [-16.818342, 37.383713, 58.145166, 48.656103, 20.159831, 3.847961, 0,0,0]
a_i_H_100_1000 = [-24.311209, 3.5692468, -11.332860, -27.850082, -21.328264, -4.2519023, 0,0,0]
a_i_H_1000_6000 = [-24.311209, 4.6450521, -3.7209846, 5.9369081, -5.5108047, 1.5538288, 0,0,0]
a_i_p_10_10000 = [-22.089523,  1.5714711, 0.015391166, -0.23619985, -0.51002221, 0.32168730, 0,0,0]
a_i_e_10_200 = [-21.928796,  16.815730,96.743155,  343.19180, 734.71651,  983.67576, 801.81247, 364.14446,  70.609154]
a_i_e_200_10000 = [-22.921189,  1.6802758, 0.93310622, 4.0406627, -4.7274036, -8.8077017, 8.9167183, 6.4380698, -6.3701156]
a_i_h2_100_6000 = [-23.962112, 2.09433740, -0.77151436, 0.43693353, -0.14913216, -0.033638326, 0,0,0]
a_i_HDL = [-20.584225, 5.0194035, -1.5739905, -4.7155769, 2.4714161,5.4710750, -3.9467356, -2.2148338, 1.8161874]


## Limiting Temperatures for Coefficients##
T_H_lim_1 = 100
T_H_lim_2 = 1000
T_H_lim_3 = 6000
T_p_lim = 1000000000000
T_e_lim_1 = 500
T_e_lim_2 = 100000000000
T_H2_lim = 6000
global_lower = 10 

function fitting_function(coeff, T, mode)
    T3 = T/1e3
    logt3 = log10(T3)
    logt32 = logt3^2
    logt33 = logt3^3
    logt34 = logt3^4
    logt35 = logt3^5
    logt36 = logt3^6
    logt37 = logt3^7
    logt38 = logt3^8
    a1 = coeff[1]
    a2 = coeff[2]
    a3 = coeff[3]
    a4 = coeff[4]
    a5 = coeff[5]
    a6 = coeff[6]
    a7 = coeff[7]
    a8 = coeff[8]
    a9 = coeff[9]

    if mode==1
        return ((a1) + (a2 * logt3) + (a3 * logt32) + (a4 * logt33) + (a5 * logt34) + (a6 * logt35) + (a7 * logt36) + (a8 * logt37) + (a9 * logt38))
    else
        return (a1) + (a2 * logt3) + (a3 * logt32) + (a4 * logt33) + (a5 * logt34) + (a6 * logt35)
    end
end

sigmoid(x, x_0, s) = 10/(10 + exp(-s*(x-x_0)))

function wCool(logT, logTmin, logTmax)
    x=(logT-logTmin)/(logTmax-logTmin)
    wcool = 10 ^ (200 * (sigmoid(x, -0.2, 50) * sigmoid(-x, -1.2, 50) - 1))
    if wcool < 1e-199
        return 0
    else
        return wcool
    end
end






