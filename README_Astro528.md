# Structure of this ReadMe #

- The Atomic Dark Matter Model
  * The Model
  * Review of the Parameters
- This Project
  * Goals
  * The Computational Challenges
  * Code Structure
  * How to Use
- History

# The Atomic Dark Matter Model #
## The Model ##
In the Atomic Dark Matter Model (aDM), some fraction of dark model is made up of a heavy dark particle, a "dark proton", and a lighter dark particle, a "dark electron". These two particles interact via a dark gauge boson, a "dark photon", at a strength governed by a dark fine structure constant ([0810.5126](https://arxiv.org/abs/0810.5126) & [0909.0753](https://arxiv.org/abs/0909.0753)). These interactions allow for dark matter halos to radiate away some of its energy, via the formation/evolution of dark bound states (i.e. dark "atomic and molecular hydrogen") through complex dark chemistry ([1705.10341](https://arxiv.org/abs/1705.10341), [1707.03829](https://arxiv.org/abs/1707.03829), [2106.13245](https://arxiv.org/abs/2106.13245), [2110.11964](https://arxiv.org/abs/2110.11964), [2110.11971](https://arxiv.org/abs/2110.11971)), possibly alleviating the tensions of the Cold Dark Matter Model ($\Lambda$CDM) on small scales ([1707.04256](https://arxiv.org/abs/1707.04256)). Moreover, this introduces the possibiliy for dark matter to form unique compact objects (i.e. sub-solar mass black holes ([1802.08206](https://arxiv.org/abs/1802.08206), [2009.05209](https://arxiv.org/abs/2009.05209), [2209.00064](https://arxiv.org/abs/2209.00064), Possible Candidate: [2001.01761](https://arxiv.org/abs/2001.01761) & [2309.11340](https://arxiv.org/abs/2309.11340) and dark neutron stars ([2201.05626](https://arxiv.org/abs/2201.05626))) and/or effect the formation/evolution of Standard Model (SM) small scale structure (i.e. POP III Stars ([2309.05758](https://arxiv.org/abs/2309.05758)), Milky Way Sub-Structure ([2304.09878](https://arxiv.org/abs/2304.09878)), and Black Holes ([1017.03419](https://arxiv.org/abs/1017.03419), [1812.03104](https://arxiv.org/abs/1812.03104), [2208.08557](https://arxiv.org/abs/2208.08557))). 

Because of the model's similarity to the SM, we can use much of the same tools just rescaled by the model parameters ([2106.13245](https://arxiv.org/abs/2106.13245), [2110.11964](https://arxiv.org/abs/2110.11964), [2110.11971](https://arxiv.org/abs/2110.11971)).

## Review of the Parameters Used to Study the Collapse of aDM Halos ##
- n: The density of the aDM "gas"
- T: The temperature of the aDM gas
- z: The redshift
- r_m: The ratio of the dark electron mass to the SM electron mass
- r_M (or r_MH): The ratio of the dark proton mass to the SM proton mass
- r_alpha: The ratio of the dark fine structure constant to the SM fine structure constant
- epsilon: The fraction of dark matter composed of aDM
- xi: The ratio of the temperature of the dark photon background to the CMB
   * This is important because aDM will have a complex cosmological evolution which places massive constraints on the model ([1209.5752](https://arxiv.org/abs/12909.57572), [1310.3278](https://arxiv.org/abs/1310.3278), [2110.11964](https://arxiv.org/abs/2110.11964), [2209.05209](https://arxiv.org/abs/2209.05209)), but this isn't relavent to this class. Xi also has a massive influence on the dark chemical abundances: low xi -> The Dark Recombination Epoch happened way earlier -> More dark hydrogen can form, and vice-versa.)
























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


