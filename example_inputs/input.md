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
SNN 200.0

#impact parameter
bmin  6.2999
bmax  6.3001


#If [Optical Glauber] is choosen fill below :
#[Info] currently in optical glauber model only same species collision is possible
species Au
opt_eps0 5.0


#If [MC Glauber] is choosen fill below :
projectile U
target U
mc_eps0 5.0
DELTA 0.4



#Hydro grid setup :
#[Info] set etamax = +1E-15, etamin =-1E-15 and neta = 1 for 2+1D evolution
xmax  +12.0 
xmin  -12.0 
ymax  +12.0
ymin  -12.0
etamax  +1e-15
etamin  -1e-15
nx  241
ny  241
neta 1


#if 3+1D choosen (if neta != 1)
eta_platue  2.59
eta_fall 0.4



#dtau, hydro starting time, freezeout condition, max run time of hydro :
dtau 0.02
tau0 0.6
Tfreeze  0.140
tauMax  30.0


# eta/s and zeta/s #
#[Info] zetaS is not a value it's a 'flag'
#[Info] set zetas 0 to turn off bulk else set 1.
etas_flag 0
zetas_flag 0
t_etas_flag 1
etas  0.08

#save in hydro file
save_every_N_steps 5000

# avg. hypersurface over :
skip_fo_tau 4
skip_fo_x 4
skip_fo_y 4
skip_fo_eta 2


# :: END :: #



