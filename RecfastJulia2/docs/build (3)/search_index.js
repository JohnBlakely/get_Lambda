var documenterSearchIndex = {"docs":
[{"location":"#","page":"-","title":"-","text":"CurrentModule = Recfast","category":"page"},{"location":"#","page":"-","title":"-","text":"\r\nBoltzmann\r\nSahaBoltz\r\nalphaH_func\r\nalphaHe_func\r\n\r\nH_z_loc\r\nNH\r\nTCMB\r\ncalc_Orel\r\n\r\nVariables_Fund_Consts\r\nParams\r\nSet_VFC_params!\r\nEvaluate_recombination","category":"page"},{"location":"#Main.Recfast.Boltzmann","page":"-","title":"Main.Recfast.Boltzmann","text":"Boltzmann(gj, gi, E, T) -> Any\n\n\nEvaluate the boltzmann factor given the statistical weight as fracg_jg_i exp-Ek_b T\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.SahaBoltz","page":"-","title":"Main.Recfast.SahaBoltz","text":"SahaBoltz(gi, gc, ne, E_ion, T) -> Any\n\n\nEvaluate the Saha equation given the statisical weights Calculate the Boltzmann factor for detailed balance as     fracn_in_c = n_efracg_i2g_c     left(frach^22 π m_e k_b Tright)^32 exp(E_rm ionk_b T)     given the statistical weights, electron number density, ionization energy     and temperature. All units are SI.\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.alphaH_func","page":"-","title":"Main.Recfast.alphaH_func","text":"alphaH_funct(TM)\n\nCalculate the case B recombination rate as a function of the matter temperature     in m^3/s as     a_1 * 10^-19 (T1e4)^a_2(1 + a_3 T^a_4)     (untested so far)\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.alphaHe_func","page":"-","title":"Main.Recfast.alphaHe_func","text":"Calculate the case B recomination rate for netural He in m^3/s     using the formula fraca_1sqrtT_mT_0 * (1 + sqrtT_m T_0)^1-a_2     (1+sqrtT_mT_1)^1 + a_2\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.H_z_loc","page":"-","title":"Main.Recfast.H_z_loc","text":"H_z_loc(params, z)\n\nWe calculate  the neutrino fraction as F_nu = fracrho_nurho_gamma = N_eff * frac7 8 * left(frac411right)^43\n\nthe redshfit at matter radiation equality as z_eq = frac3 * (H_0 c)^2 * Omega_m 8 pi G a (1 + F_nu) * T_0^4  - 1\n\nwhere a is the radiation constant and T_0 is the present temperature of the CMB. Then,\n\nH(z) = H_0 sqrtOmega_Lambda + (1 + z)^2 leftOmega_K + (1+z) Omega_M left(1 + frac1 +z1+ z_eqright)right\n\nNote, H(z) is given in rm s^-1\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.NH","page":"-","title":"Main.Recfast.NH","text":"NH(params, z)`\n\nWe calculate the hydrogen number density in rm m^-3 as:\n\n3 H_0^2 fracOmega_B8 pi G m_H mu_H (1+z)^3 where mu_H is 1(1-Y_p).\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.TCMB","page":"-","title":"Main.Recfast.TCMB","text":"TCMB(params, z)`\n\nReturn T_0 * (1 + z)\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.calc_Orel","page":"-","title":"Main.Recfast.calc_Orel","text":"calc_Orel(TCMB0, Nnu, h100)`\n\nCalculate the relativistic contribution as\n\nfraca T_0^4c^23H_0^2  (8 pi G) (1 + frac78 left(frac411right)^43N_eff )\n\nThis is used to calculate the dark energy when not specified.\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.Variables_Fund_Consts","page":"-","title":"Main.Recfast.Variables_Fund_Consts","text":"Struct to control the rescaling. Fields:\n\nmode\nDefault: 0\nRF_aScale\nDefault: 1.0\nRF_mScale\nDefault: 1.0\nfT\nDefault: 1.0\nfl\nDefault: 1.0\nfA2g\nDefault: 1.0\nfC\nDefault: 1.0\nfAp\nDefault: 1.0\nfBp\nDefault: 1.0\nfm\nDefault: 1.0\n\n\n\n\n\n","category":"type"},{"location":"#Main.Recfast.Params","page":"-","title":"Main.Recfast.Params","text":"Struct hodling all the cosmology parameters. Fields:\n\nYp\nDefault: 0.24\nT0\nDefault: 2.725\nOmega_M\nDefault: 0.26\nOmega_B\nDefault: 0.044\nOmega_K\nDefault: 0.0\nh100\nDefault: 0.71\nn_eff\nDefault: 3.04\nF\nDefault: 1.14\nf_ann\nDefault: 1.0e-24\ncorrection\nDefault: true\nb0\nDefault: 3.0\nnb\nDefault: -2.9\nrescaling\nDefault: VariablesFundConsts()\nfHe\nDefault: Yp / (4.0 * RFfacmHemH * (1 - Yp))\nH0\nDefault: (h100 * 100 * 100000.0) / RF_Mpc\nOmega_L\nDefault: ((1.0 - OmegaK) - OmegaM) - calcOrel(T0, neff, h100)\n\n\n\n\n\n","category":"type"},{"location":"#Main.Recfast.Set_VFC_params!","page":"-","title":"Main.Recfast.Set_VFC_params!","text":"Set_VFC_params!(r::Main.Recfast.Variables_Fund_Consts; aS, mS, mode)\n\n\nFunction to set the parameters fields. Mode: 0 - no rescaling\n\n1 - Rescale boltzmann factor exponentials (temperatures)\n\n2 - Rescale Thomson cross section\n\n3 - Rescale 2s1s 2 photon rate\n\n4 - Rescale alpha and beta coefficients\n\n5 - Rescale Lyman alpha rates\n\n6 - Rescale everything\n\n\n\n\n\nSet_VFC_params!(p::Main.Recfast.Params; aS, mS, mode)\n\n\n\n\n\n\n","category":"function"},{"location":"#Main.Recfast.Evaluate_recombination","page":"-","title":"Main.Recfast.Evaluate_recombination","text":"Evaluate_recombination(p::Main.Recfast.Params; zstart, zend) -> Any\n\n\nEvaluate the recombination history. Note that we assume everything is completely ionized at zstart. If this is not the case, we will get wrong results.\n\n\n\n\n\n","category":"function"}]
}
