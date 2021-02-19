#pragma once
#pragma once
#include <iostream>
#include "fluid.h"
#include "eos.h"

using namespace std;

class gubser {

public :
gubser(){}

void set_ic(fluid* f, EoS* eos, double tau, bool is_first_time) // setting initial condition
{

  double eps_0 =1.0;
  double n_0 = 0.5;
  double q = 1.0 ;
  cell* c;
  
  cout<<"After time "<<tau<<" gubser analytical data is recorded"<<endl;
  ofstream File1;
  char name[200];
  sprintf(name,"analytical_dist_%f.dat",tau);	
  File1.open(name);
  
  
  for(int i=0; i<f->get_nx(); i++)
      for(int j=0; j<f->get_ny(); j++)
         for(int k=0; k<f->get_nz(); k++){
          c = f->get_cell(i,j,k);            
	  double x_ = f->get_x(i);
	  double y_ = f->get_y(j);
          double rt = TMath::Sqrt(x_*x_+y_*y_);
          double eps = ((eps_0*TMath::Power(2.0*q,8.0/3.0))/(TMath::Power(tau,4.0/3.0)))*(TMath::Power((1+2*q*q*(tau*tau+rt*rt)+q*q*q*q*TMath::Power(tau*tau-rt*rt,2)),-4.0/3.0));
          double nb = ((n_0*4*q*q)/tau)*(TMath::Power((1+2*q*q*(tau*tau+rt*rt)+q*q*q*q*TMath::Power(tau*tau-rt*rt,2)),-1.0));
	  
          double k_ = TMath::ATanH((2.0*q*q*tau*rt)/(1+q*q*tau*tau+q*q*rt*rt));
          double vx =  (x_/rt)*(TMath::TanH(k_)); 
          double vy =  (y_/rt)*(TMath::TanH(k_)); 
          double vz = 0.0;
          double nq = 0.; double ns=0.0;
	  
	  if(j == int(f->get_ny()/2)){File1<<f->get_x(i)<<"\t"<<f->get_y(j)<<"\t"<<rt<<"\t"<<eps<<"\t"<<nb<<"\t"<<vx<<"\t"<<vy<<endl;}
	  if(is_first_time)c->set_prim_var(eos,tau,eps, nb, nq,  ns,  vx,  vy,  vz);
        }

  
  File1.close();
}


};
