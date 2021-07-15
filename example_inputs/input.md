#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::                                                                                           :::
#:::  Set your input parameters here for hydro simulation.                                     :::
#:::  The code will not read any line of this file starting with symbol "#" or any empty line. :::
#:::                                                                                           :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Eq. of state
eos 1


#initial condition mode
#[Info] ic_mode == -2/-1/0/1/2 (read gubser form file/gubser/ optical Glauber/ MC Glauber/ read from file) 
ic_mode 0


# if ic_mode 2 choosen
init_file_name example_inputs/optical_glauber_ic_dist.dat 


#collision energy
#[Info] in GeV
SNN 2760.0

#impact parameter
bmin  4.0899
bmax  4.0901


#If [Optical Glauber] is choosen fill below :
#[Info] currently in optical glauber model only same species collision is possible
species Pb
opt_eps0 25.8


#If [MC Glauber] is choosen fill below :
projectile U
target U
mc_eps0 5.0
DELTA 0.4



#Hydro grid setup :
#[Info] set etamax = +1E-15, etamin =-1E-15 and neta = 1 for 2+1D evolution
xmax  +10.0 
xmin  -10.0 
ymax  +10.0
ymin  -10.0
etamax  +10
etamin  -10
nx  201
ny  201
neta 101


#if 3+1D choosen (if neta != 1)
eta_platue  3.4
eta_fall 2.0



#dtau, hydro starting time, freezeout condition, max run time of hydro :
dtau 0.02
tau0 0.4
Tfreeze  0.150
tauMax  30.0


# eta/s and zeta/s #
#[Info] zetaS is not a value it's a 'flag'
#[Info] set zetas 0 to turn off bulk else set 1.
etas_flag 1
zetas_flag 0
t_etas_flag 0
etas  0.095

#save in hydro file
save_every_N_steps 5000

# avg. hypersurface over :
skip_fo_tau 4
skip_fo_x 4
skip_fo_y 4
skip_fo_eta 2


# :: END :: #


bmin  2.24
bmax  2.26
tau0  0.4
Tfreeze  0.150
