#pragma once

#include <iostream>
#include<TMath.h>
#include<TRandom3.h>
#include<TF1.h>
#include<fstream>
#include<string>
#include<sstream>
#include "grid.h"
#include "eos.h"
#include "idb.h"

using namespace std;

class opt_glau
{
  
 public:
  opt_glau(idb *IDB);
  ~opt_glau();
  double npart();
  double ncoll();
  void set_ic(grid* , EoS* eos);
  void opt_glau_calc();
  
  
 private:
  idb *IDB;
  double beta2;
  double beta4;
  double a ;
  double R ;
  double sigma;
  double A;
  double B;
  double n_pp = 2.25;
  double X_hard = 0.14;
  
  double mThetaA =0;
  double mPhiA = 0;
  double mThetaB =0;
  double mPhiB = 0;
  double s0 ;
  double pi = 3.1415927;
  double b;
  double norm_const;
  double bmin;
  double bmax;
  
  double Modf_WoodSaxon(const double *x);
  double Norm_Function(double* x,double* p);

  double Get_Norm_Constant();
  double npartxy(double x,double y);
  double ncollxy(double x,double y);
  double theta(double _a);
  
  double npart_tab(double* x,double* p);
  double ncoll_tab(double* x,double* p);
  
  void set_opt_glau_params()
  { 
    if(IDB->species == "Au") {A= 197; B=197; R = 6.38; a=0.535; beta2 = 0.0; beta4 =0.0;} // see ideal hdro MUSIC paper or Bozek's ideal hydro paper
    else if(IDB->species == "U") {A =238; B = 238; R = 6.81; a=0.53; beta2 = 0.28; beta4 =0.093;}
    else if(IDB->species == "Pb") {A =208; B = 208; R = 6.66; a=0.45; beta2 = 0.0; beta4 =0.0;} //see arXiv: 1409.8375 [Table. 1]
    else {cout<<"species not recognized, it's : "<<IDB->species<<endl; exit(1);}

    s0 = IDB->entropy_scale_factor;
    bmin = IDB->bmin; bmax = IDB->bmax;
    if(IDB->SNN == 62.4){sigma = 3.155;} 
    else if(IDB->SNN == 200.0){sigma = 4.2;} 
    else if(IDB->SNN == 2760.0){sigma = 6.4;} 
    else if(IDB->SNN == 5500.0){sigma = 7.2 ;} 
    else{cout<<"SNN(energy) not recognised"<<endl; exit(1);}
  }
  
  
};
