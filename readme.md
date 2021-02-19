  **************************************
# 3+1D relativistic viscous hydro code *
  **************************************

 This code is inspired from the reference below.
                                                                             
[vHLLE] : a 3D viscous hydrodynamic code
by Iurii Karpenko  et. al.
Comput. Phys. Commun. 185 (2014), 3016 [arXiv:1312.4160]	
                                                                                                                                                       
The original version of vHLLE can be found in
[ https://github.com/yukarpenko/vhlle ]


# In this code, What's from vHLLE ?

 1. The classes cnvrt.cpp, cell.cpp, ideal calculation part 
    of hydro.cpp is taken from vHLLE .

 2. Algorithm to solve the evolution 
    is same as vHLLE. 


# 2 major differences between this code and vHLLE :

 1. In the hydro viscous part calculation, the gmunu metric
    tensor = { 1.0, -1.0, -1.0, -1/pow(tau,2) }
    where as in vHLLE gmunu = {1.0, -1.0, -1.0, -1.0}.
    Due to this, terms are different than vHLLE. 
 
 2. The hydro regulation scheme to make ideal part of T^{\mu \nu}
    to be always greater than viscous part of the same is 
    also different. The modified regulation scheme is as of
    MUSIC. One can find this in hrr.cpp file. 

#  some extras here :

 1. temperature dependent bulk viscosity is 
    hard coded. The form of zeta/s(T) is same
    as taken in MUSIC.

 2. Tempertaure dependent eta/s is also introduced.
    turn on t_etas_flag in input file to take it.

 3. Latttice EoS is taken (s95pv1). set eos = 1 in input
    file to take it during hydro evolution.

 4. The external input file that one take in MUSIC
    for initial condition setup can be used here as well
    because the reading format is same.

 5. In hypersurface finding one can skip steps.

#  IC MODELS (not from vHLLE) :
   optical and monte carlo Glauber model
   or one can read an external file.


# How to give run ?

  1. Go to example inputs folder, see one example input file.
     Then make a input.md file.
  2. (In terminal) make clean
  3. (In terminal) make -j2
  3. (In terminal) ./vischydro 1 example_inputs/input.md


#  Current status ...











