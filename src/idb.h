#pragma once
#include<iostream>
#include<string>
using std::string;

// IDB : 'I'nput 'D'ata 'B'ase
 
typedef struct IDB
{
  int nx;
  int ny;
  int neta;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double etamin;
  double etamax;
  double dx;
  double dy;
  double deta;
  double dtau;
  double tau0;
  int eps_freeze_flag;
  double eps_freeze ; 
  double Tfreeze;
  double tauMax;
  double ic_mode;
  double SNN;
  
  string species;
  string projectile;
  string target;  

  double DELTA;
  double entropy_scale_factor;
  double bmin,bmax;

  int eos ; 

  double eta_platue;
  double eta_fall ;
  
  double impact_parameter;

  int save_every_N_steps;
  string init_file_name;     // filename of the initial density profile

  int etas_flag; 
  int zetas_flag; 
  int t_etas_flag; 
  double etas;

  int skip_fo_tau; 
  int skip_fo_x; 
  int skip_fo_y; 
  int skip_fo_eta; 

  
}idb;


