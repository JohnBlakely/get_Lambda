#!/usr/bin/env julia

# ============================================================================
# Calculate the visibility function from the Recfast result
# Nov 2018
# Donghui Jeong
# ============================================================================

using PyPlot
using LaTeXStrings
using DelimitedFiles
using HCubature
using Dierckx

# ----------------------------------------------------------------------------
function read_Recfast_data(fname)
    frecfast = readdlm(fname,Float64)
    az = frecfast[:,1]
    aXe = frecfast[:,2]
    aTM = frecfast[:,5]

    return az, aXe, aTM
end
# ----------------------------------------------------------------------------
function gz(az, aXe, aTM)
    #fname = "Xe_Recfast.alp.me.all.dat"
    #az, aXe,aTM = read_Recfast_data(fname)

    fXez = Spline1D(reverse(az),reverse(aXe),k=5)

    # parameters
    Mpc2cm  = 3.0857e24
    Or      = 8.54e-5
    sigmaDT = 1.66e-23 # unit = cm^2
    # functions
    nXtotz = z -> 8.60e-8*(1+z)^3 # unit = cm^3
#    Hz     = z -> 0.00033356*0.7/Mpc2cm*Or^0.5*(1+z)^2 # unit = 1/cm
    Hz = z -> 0.00033356*0.7/Mpc2cm*(0.7+0.3*(1+z)^3+Or*(1+z)^4)^0.5 # unit = 1/cm
    dtaudz = z-> nXtotz(z)*fXez(z)/(1+z)/Hz(z)

    zarr = 10 .^range(0.0,4.7,length=4001)
    gz = z -> sigmaDT*dtaudz(z)*exp(-sigmaDT*hquadrature(dtaudz,0,z)[1])
    gzarr = [gz(z) for z in zarr]


    maxindx = argmax(gzarr)
    azfine = range(zarr[maxindx-1],zarr[maxindx+1],length=100)
    agzfine = [gz(z) for z in azfine]
    zpeak = azfine[argmax(agzfine)]
    @show zarr[maxindx-1], zarr[maxindx], zarr[maxindx+1]
    @show zpeak

    fig = figure()
    ax  = subplot(111)
    ax[:set_xlabel](L"$z$")
    ax[:set_ylabel](L"$g(z)$")
    ax[:set_xscale]("linear")
    ax[:set_yscale]("linear")
    ax[:set_xlim]([1,50000])
    ax[:set_ylim]([1e-8,3e-4])

    ax[:plot](zarr,gzarr,color="DarkRed",label="visibility function julia")
    ax[:legend](loc="upper left")
    ion()
    show()
    savefig("Recfast_gz.pdf")
    println("Type something:");readline()


end
