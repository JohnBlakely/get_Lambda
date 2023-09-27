# Structure of this ReadMe #

- The Atomic Dark Matter Model
  * The Model
  * Review of the Parameters
- This Project
  * Goals
  * The Computational Challenges
     + The Challenges
     + The Solution (Attempts)
  * Code Structure
  * How to Use
- History

# The Atomic Dark Matter Model #
## The Model ##
In the Atomic Dark Matter Model (aDM), some fraction of dark model is made up of a heavy dark particle, a "dark proton", and a lighter dark particle, a "dark electron". These two particles interact via a dark gauge boson, a "dark photon", at a strength governed by a dark fine structure constant ([0810.5126](https://arxiv.org/abs/0810.5126) & [0909.0753](https://arxiv.org/abs/0909.0753)). These interactions allow for dark matter halos to radiate away some of its energy, via the formation/evolution of dark bound states (i.e. dark "atomic and molecular hydrogen") through complex dark chemistry ([1705.10341](https://arxiv.org/abs/1705.10341), [1707.03829](https://arxiv.org/abs/1707.03829), [2106.13245](https://arxiv.org/abs/2106.13245), [2110.11964](https://arxiv.org/abs/2110.11964), [2110.11971](https://arxiv.org/abs/2110.11971)), possibly alleviating the tensions of the Cold Dark Matter Model ($\Lambda$CDM) on small scales ([1707.04256](https://arxiv.org/abs/1707.04256)). Moreover, this introduces the possibiliy for dark matter to form unique compact objects (i.e. sub-solar mass black holes ([1802.08206](https://arxiv.org/abs/1802.08206), [2009.05209](https://arxiv.org/abs/2009.05209), [2209.00064](https://arxiv.org/abs/2209.00064), Possible Candidate: [2001.01761](https://arxiv.org/abs/2001.01761) & [2309.11340](https://arxiv.org/abs/2309.11340) and dark neutron stars ([2201.05626](https://arxiv.org/abs/2201.05626))) and/or effect the formation/evolution of Standard Model (SM) small scale structure (i.e. POP III Stars ([2309.05758](https://arxiv.org/abs/2309.05758)), Milky Way Sub-Structure ([2304.09878](https://arxiv.org/abs/2304.09878)), and Black Holes ([1017.03419](https://arxiv.org/abs/1017.03419), [1812.03104](https://arxiv.org/abs/1812.03104), [2208.08557](https://arxiv.org/abs/2208.08557))). 

Because of the model's similarity to the SM, we can use much of the same tools just rescaled by the model parameters ([2106.13245](https://arxiv.org/abs/2106.13245), [2110.11964](https://arxiv.org/abs/2110.11964), [2110.11971](https://arxiv.org/abs/2110.11971)).

## Review of the Parameters Used ##
- n: The density of the aDM "gas"
- T: The temperature of the aDM gas
- z: The redshift
- r_m: The ratio of the dark electron mass to the SM electron mass
- r_M (or r_MH): The ratio of the dark proton mass to the SM proton mass
- r_alpha: The ratio of the dark fine structure constant to the SM fine structure constant
- epsilon: The fraction of dark matter composed of aDM
- xi: The ratio of the temperature of the dark photon background to the CMB
   * This is important because aDM will have a complex cosmological evolution which places massive constraints on the model ([1209.5752](https://arxiv.org/abs/12909.57572), [1310.3278](https://arxiv.org/abs/1310.3278), [2110.11964](https://arxiv.org/abs/2110.11964), [2209.05209](https://arxiv.org/abs/2209.05209)), but this isn't relavent to this class. Xi also has a massive influence on the dark chemical abundances: low xi -> The Dark Recombination Epoch happened way earlier -> More dark hydrogen can form, and vice-versa.
 - M_Halo: The aDM Halo Mass
   * Important for some of the cosmology stuff in the code, not used in the parts relevant to the class
 - x_e: The dark electron abundance
 - x_p: The dark proton abundance
 - x_H: The dark atomic hydrogen abundance
 - x_Hm: The dark $H^-$ abundance
 - x_H2: The dark molecular hydrogen abundance
 - x_H2p: The dark $H_2^+$ abundance
 - x_H3p: The dark $H_3^+$ abundance
 - x: The dark ionization fraction
 - n_bar: Background aDM density
 - T_g: Background temperature of the aDM gas
   * n_bar and T_g are used to set the lower bound on n and T at a given redshift. For simplicity, points which lie below this will be set to the background chemistry at that redshift
 

# This Project #
## Goals ##
The big-picture goal of this project is to determine what parameter values let aDM halos cool, collapse, and maybe form compact objects, in a way that satisfies a number of observational constraints.

However the goal of this project in the scope of this class is to calculate the dark chemistry at _every point_ in (n, T, z, r_m, r_M, r_alpha, epsilon, xi) space, accurately and efficiently. Then these results will be put into a data cube which the code will then interpolate over (again accurately and efficiently) to get the chemistry for any choice of those parameters.

## The Computational Challenges ##
### The Challenges ###
So the code will need to calculate the Cooling Rate, the Dark Electron Fraction, and the Dark Molecular Hydrogen Fraction at every point in the parameter space. Given that each of these will be Float64, and that the resolution required is, around, - number of elements - (35, 35, 20, 20, 20, 20, 10, 10), the total memory needed is 470 BILLION BYTES (470 GB). 

Moreover, this massive monstrosity will need to be interpolated over.

On top of that, the code will need to be stoppable and restartable, because it's expected to take past 48 hrs (because I'm an undergrad and have negative funding, I can't buy more resources). This means that there will need to be a method to dump the data, know where to start again, and be able to construct the data into a single 8d data grid

And then the final big issue, the dark chemistry solver we have (DarkKROME: [2110.11971](https://arxiv.org/abs/2110.11971)) is slow and innacurate for large parts of the parameter space boundaries. If it was used as is for this project it would:

 1. Write out a ton of files each of the ~50,000 times it'll have to be compiled
 2. Take over a year, maybe. Considering running it over 1/400 of the parameter space took ~1.5 days
 3. Not be thread safe at large since four of the eight dimensions will need to be recompiled at every point; and everytime it is recompiled, it writes in a number of files

Plus it would make this project like 85% Fortran 90.

### The Solution (Attempt) ###
So the central part of this code is going to be the same used in [2309.05758](https://arxiv.org/abs/2309.05758), where they modified, an already modified [2110.11964](https://arxiv.org/abs/2110.11964) version of Recfast [astro-ph/9909275](https://arxiv.org/abs/astro-ph/9909275) (The first modification was done by James Gurian, who rewrote it in Julia and for the aDM model) to not only calculate the evolution of the dark recombinatination epoch, but also evolve th chemistry at constant density. 

Here are the major modifications that will be done to the actual chemistry solver from James Gurian:

1. It will be modified to run at constant T and n for the age of the universe at that z
2. For a given n and T, it will return the required values at each every z. So it won't need to restarted at each z
3. (_Maybe_) The chemistry will be even further simplified by following [astro-ph/9603007](https://arxiv.org/abs/astro-ph/9603007), and getting rid of every ODE for the chemical abundances EXCEPT the dark x_e, x_p, x_H, x_H2. Then by requiring charge neutrality (x_e = x_p = x) and that x_H = (1-x-2x_H2), we can go from the 9 ODEs that DarkKROME uses (techinically it's 21, but only 9 are used for aDM) to only 2

Of course this doesn't include the actual code which executes this, and manages the results. The general overview of the code which generates the tables is the following:

- There will be code that runs the chemistry evolution at a single (n, T) (from the ICs calculated higher up) and return the (n, T, z) result
  * Then after this, the result will be written into an individual HDF5 dataset, titles by its point in the 5d parameter space. This will then be compressed using GZIP. Finally it will update an outside (still in the HDF5 file) 5d array of binary types from 0 to 1, at the corresponding coordinate
     + This will make sure that at any point we don't have much more than 10-20 MB of allocated memory
     + Also it will keep the HDF5 from exploding in size
     + Finally, it makes it so there doesn't need to be a designated dump functionality. And restarting will require only reading in that outside dataset, determine what cartesian indexes are 0's, and then going through the parameters and only run over the ones which are at the corresponding points.
- There will be code which takes that and assigns it to threads (maybe those of a single core, that way 38 can run on points in the 5d space) over the (n, T) space. Then instead of passing it in as a variable, it will write it out and return nothing
- There will be code which executes James's Recfast code to calculate the chemical ICs. Then it will pass in (or defined globally with respect to the parallel code, somehow) the initial abundances, n_bar, and T_g, to the parallel code. Then, finally, it executes the parallel code.
- There will be (if I run parallel over the cpus) code which runs the previous in parallel
- In the output file, there will be the HDF5 file containing this, a text file which keeps track of the last point in 5d space succesfully ran, as well as the time, percenct done, and when the code has been stopped and restarted, and an error file containing any error messages generated and their corresponding point in 8d space
- Also some code will read in the config file and write some of its contents of attributes of the outside binary dataset

For the code which interpolates over the table, the overview is:

- There will be code which takes an input and finds the first indexes in the range of the evaluated dimensions (from the HDF5 file) then finds the first pair of points that contain the input variable (unless it equals one of them)
  * Done for all 8 parameters
- Then there will be code which opens the two + chunks which contain the inputs, then read these in as a (2,2,2,2,2,2,2,2) array.
  * I could attempt to do something more non-linear to interpolate and call in 4 total chunks
- Then there will be something that finds fractional "index" the point has wrt the new array
- Then there will finally be some code which interpolates over it



## Code Structure ##
Parent Directory is get_Lambda
- Output Directory
  * darkKROME_Equilibrium_Table.hdf5
  * log.txt
  * error_log.txt
  * stop
- config.jl
- utilities.jl
- recfast2.jl
- create_DarkKROME_Equilibrium_Tables.jl
- get_Lambda.jl
- constraints.jl
- observational_Constraints.jl

