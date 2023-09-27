# Structure of this ReadMe #

- The Atomic Dark Matter Model























This project is focused on constraining the properties of the Atomic Dark Matter Model (aDM)

In aDM, dark matter is composed of a heavy particle, a "dark proton", and a lighter particle, a "dark electron". These interact via a dark boson, a "dark photon" with a dark U(1) gauge symmetry, at a strength governed by a dark fine structure constant. Interestingly, this would allow the dark matter to radiate some of its energy (and "cool") by self-interacting and even forming bound states like dark atomic, and molecular, hydrogen. Because of this we can model the dark self-interactions using much of the same tools as the Standard Model (SM), but rescaled by the values of  
the model parameters. These are:

    r_m: the ratio of the dark electron mass to the SM electron mass
    r_M (or r_MH): the ratio of the dark proton mass to the SM proton mass
    r_alpha: the ratio of the dark fine structure constant to the SM fine structure constant
    xi: the ratio of the temperature of the dark radiation background to the CMB (since aDM is similar to the SM it also has a very complex cosmological history; in short, this parameter can be thought of as: lower = more dark hydrogen, higher = more free dark particles)
    epsilon: the fraction of DM which is aDM


Since any good model needs to agree with observations, the goal of this project is to figure out what values of those parameters make aDM agree with these observations. The observations themselves aren't relevant for this class, but are just from the matter power spectrum, galaxy cluster collisions, and that we want it to make cool things like dark black holes (using the Rees-Ostriker-Silk constraint).

So, for this project to be a success, it needs to do one very important thing:

    Calculate the dark chemistry accurately and efficiently for a collapsing halo (with a density n, temperature T, and redshift z) over the entire 5-dimensional parameter space.

This is incredibly daunting. Thankfully, we already have a chemistry solver made for aDM (DarkKROME (cite the paper)). However, it would take way too long (~ hundreds of days) to run over the parameter space, so I had to get creative.

The approach that I'm taking is to approximate the chemistry by running it until it reaches equilibrium for fixed values of n, T, and z (instead of solving the differential equations that describe them). 

The first attempt at this went fantastically (in the sense that it ran in 2 days). However, a new question arose; how long does the chemistry take to reach equilibrium? I have no idea, and I don't want to waste more resources figuring it out. So, because of this (and because it was extraordinarily inaccurate for a good amount of the parameter space) I completely scrapped this attempt.

Now I'm rewriting the chemistry solver (DarkKROME) in Julia coupled with a modified version of Recfast (cite Recfast) to calculate the background parameters (really this is just taking a chemistry solver used in this paper (cite pop3 paper) and making it for constant temperature too, then absolutely mangling it to be ran in parallel over a massive parameter space). The code has three parts:

    "generate_DarkKROME_Equilibrium_Tables.jl" (fantastic filename) This will generate the lookup tables containing the cooling rate and dark molecular fraction, which we will interpolate over in part 2. This is going to be where most of the content from this class comes in. It will need to be highly parallel, with VERY efficient memory allocation. It will also need to use HDF5 in a proficient manner to store hundreds of gigabytes of floats in a compressed, accessible manner.
    "get_Lambda.jl" This will take that data grid and access only the relevant chunks and interpolate over it. This needs to be efficient because it will be broadcasted over. 
    "constraints.jl" This will use the code in get_Lambda to determine if it satisfies the observational constraints. Not important for this class, and won't be submitted, but I'm including it for completeness.


