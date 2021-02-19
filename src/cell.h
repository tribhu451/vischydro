#pragma once
#include<iostream>
#include "TMath.h"
#include "global.h"
#include "cnvrt.h"

using namespace std;

class cell
{
  
private:
  int ix, iy, iz; //cell position
  double Q[7];    // conserved quantities at tau 
  double Qh[7];   // conserved quantities at (tau+0.5*dtau)
  double Qprev[7]; // conserved quantities at (tau-dt)
  double flux[7];  // flux of 7 conserved quantities
  double minmod(double dl, double dr);





         ////////////////////
        ///// Viscous //////
       ////////////////////



  double Pi;double PiH; double PiH0; double Pi0; 
  double piH0[10];double piH[10]; double pi[10];double pi0[10];
  int visc_correct_flag;
  
  
  cnvrt *CN;
  
  public:
  cell();
  ~cell(){};

  cell* next_cell[3];
  cell* prev_cell[3];

  inline void set_cnvrt(cnvrt* _CN){ CN = _CN; }

  inline void set_pos(int _ix, int _iy, int _iz) {ix = _ix;iy = _iy;iz = _iz;} //sets the cell position
  inline int get_ix(){return ix;} //returns cell position number along X-axis
  inline int get_iy(){return iy;} //returns cell position number along Y-axis
  inline int get_iz(){return iz;} //returns cell position number along Z-axis

  inline void set_prev_cell(int i, cell* c){prev_cell[i-1]=c;}  //sets previous cell adress. prev_cell[0] means
                                                                // previous cell along X-axis and prev_cell[1] means previous cell along Y-axis 
  inline void set_next_cell(int i, cell* c){next_cell[i-1]=c;}  //sets next cell adress

  cell* get_prev_cell(int i){return prev_cell[i-1];}     //returns previous cell adress
  cell* get_next_cell(int i){return next_cell[i-1];}     //returns next cell adress

  inline void save_Q_prev(void){ for(int i = 0; i<7; i++) Qprev[i] = Q[i]; }
  inline void clear_flux(void){ for(int i = 0; i<7; i++) flux[i] = 0.;}

  inline void get_Q(double* _Q)  {for(int i=0; i<7; i++) _Q[i]=Q[i];}
  inline void get_Qh(double* _Qh){ for(int i=0; i<7; i++) _Qh[i]=Qh[i];}


  void set_prim_var(EoS* eos,double tau, double _eps,double _nb, double _nq,
                   double _ns, double _vx, double _vy, double _vz); //set the calculational frame variables like E, Mx, My and R.
  void get_physical_var(EoS *eos, double tau, double &_e, double &_p, double &_nb, 
                        double &_nq, double &_ns, double &_vx, double &_vy, double &_vz);

  void get_left_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void get_right_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void get_left_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void get_right_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir);
  void get_center_var_prev(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz);
  void get_center_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz);
  void add_flux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns);

  void update_Q_to_Qh_by_flux();
  void update_by_flux();

 inline void set_Q(double *_Q) 
  {
   for (int i = 0; i < 7; i++) Q[i] = _Q[i];
   if (Q[T_] < 0.) 
     {
       for (int i = 0; i < 7; i++) Q[i] = 0.;
     }
   }





         ////////////////////
        ///// Viscous //////
       ////////////////////




 inline void swap(int i, int j, int &k, int &l){l=i ;  k=j; } 
 // avoid confusion !!! it's equivalent to write in code \pi^{xy}- 
 // instead of \pi^{yx}

 int indexpi(int i,int j){
    if(i>3 || j>3 || i<0 || j<0)
       {std::cout<<"error for pi index , exiting ..."<<std::endl; exit(1); return 0 ;}
    if (j>i) swap(i,j,i,j);
    if(j==0) return i ;
    if(j==1) return 4*j+(i-1) ;
    if(j==2) return 3*j + (i-1) ;
    if(j==3) {return 9;} 
    else{cout<<"PI index error"<<endl; exit(1); return 0;}
 }  // return pi^{i,j} index 


  inline void set_piH0(int i,int j, double val)
       {piH0[indexpi(i,j)] = val;}

  inline void set_pi0(int i,int j, double val)
       {pi0[indexpi(i,j)] = val;}

  inline void set_piH(int i,int j, double val)
       { piH[indexpi(i,j)] = val;}
  
  inline void set_pi(int i,int j, double val)
       {pi[indexpi(i,j)] = val;}
  

  inline void set_PiH0(double val){PiH0 = val;}
  inline void set_Pi0(double val){Pi0 = val;}
  inline void set_PiH(double _val){PiH = _val; }
  inline void set_Pi(double _val){Pi = _val; }

  inline void add_piH0(int i,int j, double val)
    { piH0[indexpi(i,j)] += val;}

  inline void add_pi0(int i,int j, double val)
    {pi0[indexpi(i,j)] += val;}
 
  inline void add_Pi0(double val){Pi0 += val;}
  inline void add_PiH0(double val){PiH0 += val;}


  inline double get_Pi(){return Pi; }
  inline double get_pi(int i,int j)
    {return pi[indexpi(i,j)];}

  inline double get_PiH0(){return PiH0; }
  inline double get_piH0(int i,int j)
    { return piH0[indexpi(i,j)];}

  inline double get_PiH(){return PiH; }
  inline double get_piH(int i,int j)
    { return piH[indexpi(i,j)];}


  inline double get_Pi0(){return Pi0; }
  inline double get_pi0(int i,int j)
    {return pi0[indexpi(i,j)];}

void update_by_visc_flux()
{
  if(fabs(flux[0]) <= 0.5*Q[0])
    {
      for (int i = 0; i < 7; i++) Q[i] += flux[i];
    } 
  else if (flux[0]!=0.)
   {
     double fac;
     fac = fabs(0.5*Q[0]/flux[0]);
     for (int i = 0; i < 7; i++) Q[i] += fac*flux[i];
   }

}


 inline void set_visc_correct_flag(int value) { visc_correct_flag = value; }
 inline int get_visc_correct_flag() { return visc_correct_flag; }



  // specially for hypersurface finding //


  double Q_fo_prev[7];
  double pi_fo_prev[10];
  double Pi_fo_prev;


  inline void get_Q_fo_prev(double* _Q) 
     {for(int i=0; i<7; i++) _Q[i] = Q_fo_prev[i];}

  inline double get_pi_fo_prev(int i,int j)
     {return pi_fo_prev[indexpi(i,j)];}

  inline double get_Pi_fo_prev()
     {return Pi_fo_prev;}

inline void get_center_var_fo_prev(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz)
{
 double _Q[7];
 for (int i = 0; i < 7; i++) _Q[i] = Q_fo_prev[i] / tau;
 CN->CALC_2_LRF(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}


 inline void save_for_fo()
 {
  for (int i = 0; i < 7; i++){ Q_fo_prev[i] = Q[i]; }
  for (int i = 0; i < 10; i++) { pi_fo_prev[i] = pi[i]; }
  Pi_fo_prev = Pi ;
 }


};








