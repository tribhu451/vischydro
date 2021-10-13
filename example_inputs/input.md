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
init_file_name hydro_output/optical_glauber_ic_dist.dat 


#collision energy
#[Info] in GeV
SNN 200.0

#impact parameter
bmin  2.3499
bmax  2.3501


#initial entropy scaling factor
entropy_scale_factor 16.8


#If [Optical Glauber] is choosen fill below :
#[Info] currently in optical glauber model only same species collision is possible
species Au



#If [MC Glauber] is choosen fill below :
projectile U
target U
DELTA 0.4



#Hydro grid setup :
#[Info] set etamax = +1E-15, etamin =-1E-15 and neta = 1 for 2+1D evolution
xmax  +10.0 
xmin  -10.0 
ymax  +10.0
ymin  -10.0
etamax  +1e-10
etamin  -1e-10
nx  201
ny  201
neta 1


#if 3+1D choosen (if neta != 1)
eta_platue  3.4
eta_fall 2.0



#dtau, hydro starting time, freezeout condition, max run time of hydro :
dtau 0.02
tau0 0.4
eps_freeze_flag 1
eps_freeze 0.24
Tfreeze  0.150
tauMax  30.0


# eta/s and zeta/s #
#[Info] zetaS is not a value it's a 'flag'
#[Info] set zetas 0 to turn off bulk else set 1.
etas_flag 1
zetas_flag 0
t_etas_flag 0
etas  0.08

#save in hydro file
save_every_N_steps 5000


# Below inputs are not required for vischydro_for_therminator
# avg. hypersurface over :
skip_fo_tau 2
skip_fo_x 2
skip_fo_y 2
skip_fo_eta 1


# :: END :: #


