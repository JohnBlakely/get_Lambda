#!/usr/bin/env julia

using PyPlot
using DelimitedFiles
using LaTeXStrings



function plotRes(what::String, az, aXe, aTM)
    #fname = "Recfast++.vDJ\\outputs\\Xe_Recfast.alp.me.all.dat"
    fname = "Recfast++.vDJ\\outputs\\Xe_Recfast.alp.me.all.dat"
    azC, aXeC,aTMC = read_Recfast_data(fname)
    if(what=="Xe")
# First, plot the recombination history
        ax  = subplot(111)
        ax[:set_xlabel](L"$z$")
        ax[:set_ylabel](L"$X_e(z)$")
        ax[:set_xscale]("log")
        ax[:set_yscale]("log")
        ax[:set_xlim]([1,100000])
        ax[:set_ylim]([1e-8,1.1])

        ax[:plot](az,aXe,color="DarkRed",label="Recfast++")
        ax[:plot]([51300,51300],[1.e-10,1],color="Black",ls=":",label="Saha calculation")

        ax[:legend](loc="upper left")
        ion()
        show()
        savefig("Recfast_Xe.pdf")
        println("Type something:");readline()

        close(fig)
    elseif (what=="Tm")
# Now plot the matter temperature
        ax  = subplot(111)
        ax[:set_xlabel](L"$z$")
        ax[:set_ylabel](L"$T(z) $[K]")
        ax[:set_xscale]("log")
        ax[:set_yscale]("log")
        ax[:set_xlim]([1,100000])
        ax[:set_ylim]([1e-6,20000])

        aTgamma = similar(az)
        @. aTgamma = 2.725*0.02*(1+az)
        aTm2 = similar(az)

        @. aTm2 = 2.725*0.02*(1+35000)*((1+az)/(35000+1))^2.

        ax[:plot](azC,aTMC, "k-",label="Dark matter temperature, C")

        ax[:plot](az,aTM, "r--",label="Dark matter temperature, JL")
        ax[:plot](az,aTgamma,color="Blue",ls="--",label="Dark photon temerature")
        ax[:plot](az,aTm2,color="Green",ls=":",label=L"Decouple at $z=35000$")

        ax[:legend](loc="upper left")
        ion()
        show()
        savefig("Recfast_TM.pdf")
        println("Type something:");readline()

    else
        println("Input valid parameter (either Xe or Tm)")
    end

end
