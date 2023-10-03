function compton_cooling(p, TM, TR, xe, nHtot)
    if TM <= 0
        return 0
    end
    4.91e-22 * 1e-13 * (TM - TR) * p.rescaling.RF_aScale^2 * p.rescaling.RF_mScale^(-3) * TR^4 * xe * nHtot
end

function atomic_cooling(p, TM, xe, xH1, nHTot)
    if TM <= 0
        return 0
    end
    alpha = RF_alpha * p.rescaling.RF_aScale
    mc_ev = RF_me_eV * p.rescaling.RF_mScale
    TM_ev = TM * RF_kb_eV
    y = sqrt(alpha^2 * mc_ev/(2*TM_ev))
    prefactor = (7.4e-18 * 1e-13 * (RF_alpha * p.rescaling.RF_aScale/.01)^2 * (p.rescaling.RF_mScale)^-0.5
        * (10. ^5/TM)^(1/2) *xe * xH1 * nHTot^2)
    uintegrand(u) = u * exp(- u^2)/(1 + (7 * y^2)/(4*u^2)) * log(4 * u/y)
    return prefactor * quadgk(uintegrand, sqrt(3)*y/2, Inf)[1] #in joules/ (m^3 s)

end

function ionization_cooling(p, TM, xe, xH1, nHTot)
    if TM <= 0
        return 0
    end
    alpha = RF_alpha * p.rescaling.RF_aScale
    mc_ev = RF_me_eV * p.rescaling.RF_mScale
    TM_ev = TM * RF_kb_eV

    y = sqrt(alpha^2 * mc_ev/(2*TM_ev))
    integrand(u) = u * exp(-u^2)/(1+ (2 * y^2/u^2)) * ( 1 - y^2/u^2 - 0.5 * (1 - y^4/u^4) * log(y^2/u^2)
        + y^2 * log(y^2/u^2)/(u^2+y^2))
    fy = quadgk(u -> integrand(u), y, Inf )[1]
    return (9e-18 * 1e-13 * (alpha/.01)^2 * (p.rescaling.RF_mScale)^(-1/2) * (1e5/TM)^(1/2) * fy
        *xe * xH1 * nHTot^2)
end

function ff_cooling(p, TM, xe, xp, nHTot)
    if TM <= 0
        return 0
    end
    alpha = RF_alpha * p.rescaling.RF_aScale
    mc_ev = RF_me_eV * p.rescaling.RF_mScale
    TM_ev = TM * RF_kb_eV

    #1e-13 from cgs to SI
    return 3.7e-27 * 1e-13 * (p.rescaling.RF_mScale)^(-3/2) * (alpha/.01)^3 * TM^(1/2) * xp * xe * nHTot^2
end

function gammaJH(J,T)
    T3 = T/1000
    return (1e-11 * sqrt(T3)/(1+60 * T3^(-4)) + 1e-12 * T3) * (0.33 + 0.9 * exp(- ((J-3.5)/0.9)^2))
end

function gammaJH2(J,T)
    T3 = T/1000
    return (3.3e-12  + 6.6e-12 * T3) * (0.276 * J^2 * exp(-(J/3.18)^(1.7)))
end
gamma10h(T) = 1e-12 * sqrt(T) * exp(-1000/T)
gamma20h(T) = 1.6e-12 * sqrt(T) * exp(-(400/T)^2)
gamma21h(T) = 4.5e-12 * sqrt(T) * exp(-(500/T)^2)

gamma10h2(T) =  1.4e-12 * sqrt(T) * exp(-12000/(T+1200))
gamma21h2(T) = 1.4e-12 * sqrt(T) * exp(-12000/(T+1200))
gamma20h2(T) = 0

function lrotLDL(gamma2, gamma3, T)
    E_rot_scale = 512/(2*3) #this is 1/(m * r) from dark chem paper, using 512 K for transition from J=2 to J=0
    delta20 = E_rot_scale * (2*(2+1)) #energy gap
    delta31 = E_rot_scale * ((3 * 4) - (2))
    delta20ergs = delta20 * 1.381e-16
    delta31ergs = delta31 * 1.381e-16
    #in erg cm^3/s
    #return 0
    return (0.25 * (5 * gamma2(T) * exp(-delta20/T) * (delta20ergs)) + 0.75 * ( 7/3 * gamma3(T) * exp(-delta31/T) *delta31ergs))

end

function lvibLDL(gamma10, gamma20, gamma21, T)
    delta10 = 5860  #energy gap between SM vibrational levels, in K
    delta21 = 5860
    delta20 = 5860*2
    delta10ergs = delta20 * 1.381e-16
    delta21ergs = delta21 * 1.381e-16
    delta20ergs = delta20 * 1.381e-16

    return gamma10(T) * exp(-delta10/T) *(delta10ergs) +
        gamma20(T) * exp(-delta20/T) * delta20ergs + gamma21(T) * exp(-delta21/T) * delta21ergs

end


function crotLTE(T, n)
    T3 = T/1000
    return n *(9.5e-22 * T3^3.76 / (1 + 0.12*T3^2.1) * exp(-(0.13/T3)^3) + 3e-24 * exp(-0.51 /T3))
end

#the per particle cooling rate in LTE from collisions with H or H2, times n(H,H2) (erg/m^3)
function cvibLTE(T, n)
    T3vib = T/1000
    return n *( 6.7e-19 * exp(-5.86/T3vib) + 1.6e-18 * exp(-11.7/T3vib))
end

function hh2cooling(p, T, nh, nh2)
     if T <= 0
        return 0
    end

    Tr = T/p.rescaling.delta_Erot
    gamma2(T) = gammaJH(2,T)
    gamma3(T) = gammaJH(3,T)

    #L from Hollenbach and McKee, units of erg cm^3 s^-1
    lrLDL = lrotLDL(gamma2, gamma3, Tr) * p.rescaling.RF_aScale * p.rescaling.RF_mScale * p.rescaling.RF_mxScale^-2
    #now the volumetric rate. We need to convert to SI
    crLDL = 1e-13 * lrLDL * nh * nh2
    crLTE = 1e-7 * crotLTE(Tr, nh2) * p.rescaling.RF_aScale^9 * p.rescaling.RF_mScale^8 * p.rescaling.RF_mxScale^-6
    if crLDL <= 0
        totalrot = 0
    else
        totalrot = crLTE/ (1 + crLTE/crLDL)
    end

    Tvib = T/p.rescaling.delta_Evib
    lvLDL =  lvibLDL(gamma10h, gamma20h, gamma21h, Tvib) * p.rescaling.RF_aScale * p.rescaling.RF_mScale^(1/4) * p.rescaling.RF_mxScale^(-5/4)
    cvLDL = 1e-13 * lvLDL * nh * nh2
    cvLTE = 1e-7 * cvibLTE(Tvib, nh2) * p.rescaling.RF_aScale^9 * p.rescaling.RF_mScale^5 * p.rescaling.RF_mxScale^-3
    if cvLDL <=0
        totalvib = 0
    else
         totalvib = cvLTE / (1 + cvLTE/cvLDL)
     end
    return totalrot + totalvib
end

function h2h2cooling(p, T, nh2)
    if T <= 0
        return 0
    end
    Tr = T/p.rescaling.delta_Erot
    gamma2(T) = gammaJH(2,T)
    gamma3(T) = gammaJH(3,T)

    #L from Hollenbach and McKee, units of erg cm^3 s^-1
    lrLDL = lrotLDL(gamma2, gamma3, Tr) * p.rescaling.RF_aScale * p.rescaling.RF_mScale * p.rescaling.RF_mxScale^-2
    #now the volumetric rate
    crLDL = 1e-13 * lrLDL * nh2^2
    crLTE = 1e-7 * crotLTE(Tr, nh2) * p.rescaling.RF_aScale^9 * p.rescaling.RF_mScale^8 * p.rescaling.RF_mxScale^-6
    if crLDL <= 0
        totalrot = 0
    else
        totalrot = crLTE/ (1 + crLTE/crLDL)
    end


    Tvib = T/p.rescaling.delta_Evib
    lvLDL = lvibLDL(gamma10h2, gamma20h2, gamma21h2, Tvib) * p.rescaling.RF_aScale * p.rescaling.RF_mScale^(1/4) * p.rescaling.RF_mxScale^(-5/4)
    cvLDL = 1e-13 * lvLDL * nh2^2

    cvLTE = 1e-7 * cvibLTE(Tvib, nh2) * p.rescaling.RF_aScale^9 * p.rescaling.RF_mScale^5 * p.rescaling.RF_mxScale^-3
    if cvLDL <=0
        totalvib = 0
    else
        totalvib = cvLTE / (1 + cvLTE/cvLDL)
     end
    #cvLDL = 0
    #1e-13 from cgs to SI
    return totalrot + totalvib
end

function hh2cooling_ldl(p, T, nh, nh2)
    if T <= 0
        return 0
    end
    Tr = T/p.rescaling.delta_Erot
    gamma2(T) = gammaJH(2,T)
    gamma3(T) = gammaJH(3,T)

    #L from Hollenbach and McKee, units of erg cm^3 s^-1
    lrLDL = lrotLDL(gamma2, gamma3, Tr) * p.rescaling.RF_aScale * p.rescaling.RF_mScale * p.rescaling.RF_mxScale^-2
    #now the volumetric rate
    crLDL = lrLDL * nh * nh2

    Tvib = T/p.rescaling.delta_Evib
    lvLDL = lvibLDL(gamma10h, gamma20h, gamma21h, Tvib) * p.rescaling.RF_aScale * p.rescaling.RF_mScale^(1/4) * p.rescaling.RF_mxScale^(-5/4)
    cvLDL = lvLDL * nh * nh2
    #cvLDL = 0
    #1e-13 from cgs to SI
    #return 1e-13 * crLDL
    return 1e-13 * (cvLDL + crLDL)
end


function h2h2cooling_ldl(p, T, nh2)
    if T <= 0
        return 0
    end
    Tr = T/p.rescaling.delta_Erot
    gamma2(T) = gammaJH2(2,T)
    gamma3(T) = gammaJH2(3,T)

    #L from Hollenbach and McKee, units of erg cm^3 s^-1
    lrLDL = lrotLDL(gamma2, gamma3, Tr) * p.rescaling.RF_aScale * p.rescaling.RF_mScale * p.rescaling.RF_mxScale^-2
    #now the volumetric rate
    crLDL = lrLDL * nh2^2

    Tvib = T/p.rescaling.delta_Evib
    lvLDL = lvibLDL(gamma10h2, gamma20h2 , gamma21h2, Tvib) *p.rescaling.RF_aScale * p.rescaling.RF_mScale^(1/4) * p.rescaling.RF_mxScale^(-5/4)

    cvLDL = lvLDL * nh2^2

    #return lvLTE + lrLTE
    #return  lvLDL + lrLDL
    return 1e-13 * (cvLDL + crLDL)

end
#function h2discooling(p, T, nh, nh2)
#    return 4.48 * 1.60218e-19 * p.rescaling.fl * h2hdis(p, T) * nh * nh2
#end









###TO Do: adiabatic and inverse compton cooling
