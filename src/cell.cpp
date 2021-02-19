#include "cell.h"



using std::cout;
using std::endl;


// cell quantities
cell::cell()
{
  for(int i = 0; i<7; i++)   //tau*T^00, tau*T^0x, tau*T^0y,
     {                       //tau*tau*T^0z, tau*nb, tau*nq, tau*ns (7 conserved quantities)
	Q[i] = 0.;
	Qh[i] = 0.;
	Qprev[i] = 0.;
	flux[i] = 0.;
      }


         ////////////////////
        ///// Viscous //////
       ////////////////////


for(int i=0; i<10; i++){pi[i] =0.0;piH0[i]=0.0; pi0[i]=0.0;}
Pi =0.0; PiH0 = 0.0; Pi0 = 0.0;

visc_correct_flag = 0;


}





//[Info] the minmod function
double cell::minmod(double a, double b) {
 if (a * b <= 0.) return 0.;
 if (fabs(a) > fabs(b))
  return b;
 else
  return a;
}


//[Info] set the conserved quantities from initial condition
void cell::set_prim_var(EoS* eos,double tau, double _eps,double _nb,
                                 double _nq, double _ns, double _vx, double _vy, double _vz)
{
  double gamma2 = 1.0/(1.0- (_vx*_vx+_vy*_vy+_vz*_vz));
  double gamma = TMath::Sqrt(gamma2);
  double p = eos->pressure(_eps,_nb, _nq, _ns);
  Q[T_] = tau*(((_eps + p)*gamma2) - p);
  Q[X_] = tau*((_eps+p)*gamma2*_vx);
  Q[Y_] = tau*((_eps+p)*gamma2*_vy);
  Q[Z_] = tau*((_eps+p)*gamma2*_vz);
  Q[NB_] = tau*_nb*gamma;
  Q[NQ_] = tau*_nq*gamma;
  Q[NS_] = tau*_ns*gamma;
}


//[Info] get physical variables like eps,pressure from conserved quantities
void cell::get_physical_var(EoS *eos, double tau, double &_e, double &_p,
                                      double &_nb, double &_nq, double &_ns, double &_vx,
                                      double &_vy, double &_vz)
{
  double _Q[7];
  for(int i =0; i<7 ; i++){_Q[i] = Q[i]/tau; }
  CN->CALC_2_LRF( eos,  _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}




//[info] get physical variables of the cell present to the left in direction-
//       "dir(x,y or z)" of the current cell 

//[format] the format of writting the name of next 4/5 functions is as follows : 
//[format] get_(position)_var_(time)
void cell::get_left_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)

{

 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->get_Q(Qr);
 prev_cell[dir - 1]->get_Q(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Q[i] - dQ[i]) / tau;
 CN->CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);

}



void cell::get_right_var(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{
 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->get_Q(Qr);
 prev_cell[dir - 1]->get_Q(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Q[i] + dQ[i]) / tau;
 CN->CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}


void cell::get_left_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{
 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->get_Qh(Qr);
 prev_cell[dir - 1]->get_Qh(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Qh[i] - dQ[i]) / tau;
 CN->CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}



void cell::get_right_varH(EoS *eos, double tau, double &_e, double &_p,
                                  double &_nb, double &_nq, double &_ns, double &_vx,
                                  double &_vy, double &_vz, int dir)
{
 double Qr[7], Ql[7], dQ[7];

 next_cell[dir - 1]->get_Qh(Qr);
 prev_cell[dir - 1]->get_Qh(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Qh[i] + dQ[i]) / tau;
 CN->CALC_2_LRF(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);

}


//[info] get the physical variables of the current cell but at previous time step
void cell::get_center_var_prev(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz)
{
 double _Q[7];
 for (int i = 0; i < 7; i++) _Q[i] = Qprev[i] / tau;
 CN->CALC_2_LRF(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}


//[info] get the physical variables of the current cell but at half time step
void cell::get_center_varH(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz)
{
 double _Q[7];
 for (int i = 0; i < 7; i++) _Q[i] = Qh[i] / tau;
 CN->CALC_2_LRF(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}




void cell::add_flux(double Ft, double Fx, double Fy, double Fz, double Fnb,
                     double Fnq, double Fns)
{

  if(std::isinf(Ft) or std::isnan(Ft)) 
    {
      std::cout<<"at pos -> "<<get_ix()<<"\t"<<get_iy()<<"\t"<<get_iz()<<" : "<<Ft<<std::endl;
      std::cout << "Cell::addFlux inf/nan\n"; exit(1);
    }
flux[T_] += Ft;
flux[X_] += Fx;
flux[Y_] += Fy;
flux[Z_] += Fz;
flux[NB_] += Fnb;
flux[NQ_] += Fnq;
flux[NS_] += Fns;
}

void cell :: update_Q_to_Qh_by_flux()
{
   for(int i = 0; i<7; i++) Qh[i] = Q[i] + flux[i];
}

void cell :: update_by_flux()
{
if(Q[0] + flux[0] < 0.) return;
for(int i = 0; i<7; i++) Q[i] += flux[i];
}



         ////////////////////
        ///// Viscous //////
       ////////////////////







