
"""

$(TYPEDSIGNATURES)
Evaluate the boltzmann factor given the statistical weight as
`` \\frac{g_j}{g_i} \\exp{-E/k_b T}``
"""
function Boltzmann(gj, gi, E, T)
    return  max(1e-300, gj / gi * exp(-E / (RF_kBoltz * T)))
end

"""
$(TYPEDSIGNATURES)
Evaluate the Saha equation given the statisical weights
Calculate the Boltzmann factor for detailed balance as
    ``\\frac{n_i}{n_c} = n_e\\frac{g_i}{2g_c}
    \\left(\\frac{h^2}{2 π m_e k_b T}\\right)^{3/2} \\exp{(E_{\\rm ion}k_b T)}``
    given the statistical weights, electron number density, ionization energy
    and temperature. All units are SI.
"""
function SahaBoltz(gi, gc, ne, E_ion, T)
    c1 = (RF_hPlanck/ (2 * pi *RF_mElect)) * (RF_hPlanck/RF_kBoltz)
    return ne * gi / (2 * gc) * c1^(3 / 2) * T^(-3/2) *
            exp(min(300, E_ion / (RF_kBoltz * T)))
end

@doc raw"""
`alphaH_funct(TM)`

Calculate the case B recombination rate as a function of the matter temperature
    in m^3/s as
    `` a_1 * 10^{-19} (T/1e4)^{a_2}/(1 + a_3 T^{a_4})``
    (untested so far)
"""
function alphaH_func(TM)
    a1 = 4.309
    a2 = -0.6166
    a3 = 0.6703
    a4 = 0.5300
    t = TM/1e4
    return a1 * 1e-19 * t^a2 / (1 + a3 * t^a4)
end

@doc raw"""


Calculate the case B recomination rate for netural He in m^3/s
    using the formula ``\frac{a_1}{sqrt{T_m/T_0} * (1 + \sqrt{T_m/ T_0})^{1-a_2}}
    (1+\sqrt{T_m/T_1})^{1 + a_2}``
"""
function alphaHe_func(TM)
    a1 = 10^-16.744
    a2 = 0.711
    T0 = 3.0
    T1 = 10^5.114
    sqrt_TMT0 = sqrt(TM/T0)
    sqrt_TMT1 = sqrt(TM/T1)
    return a1/ (sqrt_TMT0 * (1 + sqrt_TMT0)^(1 - a2) * (1 + sqrt_TMT1)^(1 + a2))
end

"""
evaluate_Recfast_System!(params, z, y::Array{Real, 1})
"""
function Evaluate_Recfast_System(p, logz, y)

    cosmoVals = calc_cosmoVals(p, 10^logz, y)
    calc_dxe(p, logz, cosmoVals)


end

"""
$(TYPEDSIGNATURES)
Calculate the densities, temperatures, and Hubble parameter
 """
function calc_cosmoVals(p, z, y)
    nHTot = NH(p, z)
    fHe = p.fHe

    #set variables
    xHep = y[1]
    xp = y[2]
    TM = y[3]
    TR = TCMB(p, z)
    xe = xp + xHep
    xH1 = 1 - xp
    nH1 = xH1 * nHTot
    nHe1 = nHTot * (fHe - xHep)
    Hz = H_z_loc(p, z)
    nTot = (1 + xe + fHe) * ( nHTot)

    return [nTot, nHTot, xHep, xp, TM, TR, xe, xH1, nH1, nHe1, Hz]
end

"""
$(TYPEDSIGNATURES)
Calculate the rate equations for xHep, xp, and TM.
"""
function calc_dxe(p, logz, cosmoVals)

    z = 10^logz
    fHe = p.fHe
    nTot, nHTot, xHep, xp, TM, TR, xe, xH1, nH1, nHe1, Hz = cosmoVals
    if TM < 0
        return @SVector [0,0,0]
    end



    #Compton term
    Comp = 8 * RF_sigmaT * RF_aRad * TR^4 /
    (3 * Hz * (1 + z) * RF_mElect * RF_cLight)
    #print("Compton ", Comp)

    A2s1sH = RF_Lam2s1sH * p.rescaling.fA2g


    lambda21H = RF_LyalphaH / p.rescaling.fl
    TMeval = TM / p.rescaling.fT

    TReval = TR / p.rescaling.fT
    Comp *= p.rescaling.fC

    #Boltzmann and detailed balance factors

    SH = SahaBoltz(2, 1, 1, RF_EionH2s, TMeval)
    BH = Boltzmann(1, 1, RF_E2s1sH , TMeval)


    # H recombination coefficient + fudge-factor
    alphaH = alphaH_func(TMeval)
    alphaH *= p.F * p.artificial[1]
    betaH = alphaH/SH
    #print("\n SH", SH, "\n")

    ## Photo-ionization rates at T=TR
    # This is not default, so we'll skip it for now

    # Rescale photo-ionization and recombination coefficients
    alphaH *= p.rescaling.fAp

    betaH *= p.rescaling.fBp

    # Sobolev escape for He and H
    # K = tauS/(A1g NH1)
    KH = lambda21H^3 / (Hz * 8 * pi)

    # Inhibition factors
    CH = (1 + KH * A2s1sH * nH1)/ (1 + KH * (A2s1sH + betaH) * nH1)
    dH = ((alphaH * xe * xp * nHTot - betaH * (xH1) * BH) * CH ) / ((1 + z) * Hz) *(z * log(10))

    A2s1sHe = RF_Lam2s1sHe * p.rescaling.fA2g
    lambda21He = RF_LyalphaHe / p.rescaling.fl
    BHe = Boltzmann(1, 1, RF_E2s1sHe, TMeval)
    SHe = SahaBoltz(1, 2, 1, RF_EionHe2s, TMeval)
    #The HeI energy levels 2p and 2s differ,
    #so there is a Boltzmann factor
    BHe2p2s = 1 / Boltzmann(1, 1, RF_E2p2sHe, TMeval)

    ## He recombination coefficient
    alphaHe = alphaHe_func(TMeval)
    betaHe = alphaHe/SHe
    alphaHe *= p.rescaling.fAp

    betaHe *= p.rescaling.fBp
    KHe = lambda21He^3 / (Hz * 8 * pi)
    CHe = (1 + KHe * A2s1sHe * nHe1 * BHe2p2s) /
    (1 + KHe * (A2s1sHe + betaHe) * nHe1 *BHe2p2s)
    dHe = (alphaHe * xe * xHep * nHTot - betaHe * (fHe - xHep) * BHe) * CHe/
        ((1 + z) * Hz)* (z * log(10))
    dTM =( Comp * xe * nHTot/ nTot * (TM - TR) + 2 * TM/(1+z)) * (z * log(10))

    return @SVector [dHe, dH, dTM]

end
