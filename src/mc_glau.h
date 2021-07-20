#pragma once

#include<fstream>
#include<string>
#include<sstream>
#include "global.h"
#include<iostream>
#include<fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "grid.h"
#include "eos.h"
#include "idb.h"
	
using namespace std;

class mc_glau
{
  
 public:
  mc_glau(idb *IDB);
  ~mc_glau();
  void set_ic(grid* f, EoS* eos);
  
  
 private:
  idb *IDB;
  
  int A; // mass no. of projrctile nucleus
  int B; // mass no. of target nucleus
  double sigma; // energy-> cross-section
  
  double p_radius; //wood-Saxon parameters       
  double p_dlt;        
  double p_beta2;
  double p_beta4;

  double t_radius; //wood-Saxon parameters       
  double t_dlt;        
  double t_beta2;
  double t_beta4;
  
  double npp=2.25;// two-component energy deposition
  double X_hard=0.14;

  
  double DELTA; //gussian smearing sigma
  double eps0 ; //energy density scaling factor
  double bmin,bmax; //impact parameter range

  double theta(double _a);

void generate_nucleus(double* X1, double* Y1,double* Z1,int A,
                double R, double dlt, double BETA2, double BETA4, double etaA, double psiA);

void calculate_npart_ncoll(double* vxA,double* vyA,double* vxB,double* vyB, int &Npart, 
       int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y);

void shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
                     double zhi,double* X2, double* Y2, double* Z2 );

void set_mc_glau_params()
{
    // projectile nucleus
    if(IDB->projectile == "Au") {A= 197; p_radius = 6.37; p_dlt=0.53; p_beta2 = 0.0; p_beta4 =0.0;}
    else if(IDB->projectile == "U") {A= 238; p_radius = 6.81; p_dlt=0.54; p_beta2 = 0.28; p_beta4 =0.093;}
    else {cout<<"projectile not recognized, it's : "<<IDB->projectile<<endl; exit(1);}
    
    // target nucleus
    if(IDB->target == "Au") {B= 197; t_radius = 6.37; t_dlt=0.53; t_beta2 = 0.0; t_beta4 =0.0;}
    else if(IDB->target == "U"){B =238;  t_radius = 6.81; t_dlt=0.54; t_beta2 = 0.28; t_beta4 =0.093;}
    else {cout<<"target not recognized, it's : "<<IDB->target<<endl; exit(1);}

    // collision energy
    if(IDB->SNN == 200.0){sigma =4.2;}else{cout<<"SNN(energy) not recognised"<<endl; exit(1);}

    // impact parameter range
    bmin = IDB->bmin; bmax = IDB->bmax;

    DELTA = IDB->DELTA; eps0 = IDB->mc_eps0; 

}



};
















