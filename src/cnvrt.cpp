#include <iostream>
#include "TMath.h"
#include "cnvrt.h"
#include "global.h"


using std::cout;
using std::endl;

cnvrt::cnvrt()
{
}

cnvrt::~cnvrt()
{
}

void cnvrt::CALC_2_LRF(EoS *eos, double* Q, double &e, double &p, double &nb,
                 double &nq, double &ns, double &vx, double &vy, double &vz) {
 // conserved -> primitive transtormation requires
 // a numerical solution to 1D nonlinear algebraic equation:
 // v = M / ( Q_t + p(Q_t-M*v, n) )       (A.2)
 // M being the modulo of the vector {Q_x, Q_y, Q_z}.
 // Bisection/Newton methods are used to solve the equation.
 const int MAXIT = 100;        // maximum number of iterations
 const double dpe = 1. / 3.;   // dp/de estimate for Newton method
 const double corrf = 0.9999;  // corrected value of M
 // when it brakes the speed of light limit, M>Q_t
 double v, vl = 0., vh = 1., dvold, dv, f, df;

 double M = sqrt(Q[X_] * Q[X_] + Q[Y_] * Q[Y_] + Q[Z_] * Q[Z_]);
 if (Q[T_] <= 0.) {
  e = 0.;
  p = 0.;
  vx = vy = vz = 0.;
  nb = nq = ns = 0.;
  return;
 }
 if (M == 0.) {
  e = Q[T_];
  vx = 0.;
  vy = 0.;
  vz = 0.;
  nb = Q[NB_];
  nq = Q[NQ_];
  ns = Q[NS_];
  p = eos->pressure(e, nb, nq, ns);
  return;
 }
 if (M > Q[T_]) {
  Q[X_] *= corrf * Q[T_] / M;
  Q[Y_] *= corrf * Q[T_] / M;
  Q[Z_] *= corrf * Q[T_] / M;
  M = Q[T_] * corrf;
 }

 v = 0.5 * (vl + vh);
 e = Q[T_] - M * v;
 if (e < 0.) e = 0.;
 nb = Q[NB_] * sqrt(1 - v * v);
 nq = Q[NQ_] * sqrt(1 - v * v);
 ns = Q[NS_] * sqrt(1 - v * v);
 p = eos->pressure(e, nb, nq, ns);
 f = (Q[T_] + p) * v - M;
 df = (Q[T_] + p) - M * v * dpe;
 dvold = vh - vl;
 dv = dvold;
 for (int i = 0; i < MAXIT; i++) {
  if ((f + df * (vh - v)) * (f + df * (vl - v)) >= 0. ||
      fabs(2. * f) > fabs(dvold * df)) {  // bisection
   dvold = dv;
   dv = 0.5 * (vh - vl);
   v = vl + dv;
   //			cout << "BISECTION v = " << setw(12) << v << endl ;
  } else {  // Newton
   dvold = dv;
   dv = f / df;
   v -= dv;
   //			cout << "NEWTON v = " << setw(12) << v << endl ;
  }
  if (fabs(dv) < 0.00001) break;

  e = Q[T_] - M * v;
  if (e < 0.) e = 0.;
  nb = Q[NB_] * sqrt(1 - v * v);
  nq = Q[NQ_] * sqrt(1 - v * v);
  ns = Q[NS_] * sqrt(1 - v * v);
  p = eos->pressure(e, nb, nq, ns);
  f = (Q[T_] + p) * v - M;
  df = (Q[T_] + p) - M * v * dpe;

  if (f > 0.)
   vh = v;
  else
   vl = v;
  if (nb != nb)
   cout << "transformCV:nbInf " << i << "  " << e << "  " << nb << "  " << nq
        << "  " << ns << "  " << p << "  " << v << endl;
  
 }  // for loop
    //----------after
    // v = 0.5*(vh+vl) ;
 vx = v * Q[X_] / M;
 vy = v * Q[Y_] / M;
 vz = v * Q[Z_] / M;
 e = Q[T_] - M * v;
 p = eos->pressure(e, nb, nq, ns);
 nb = Q[NB_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
 nq = Q[NQ_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
 ns = Q[NS_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
 
 if (e < 0. || sqrt(vx * vx + vy * vy + vz * vz) > 1.) {
  cout << Q[T_] << "  " << Q[X_] << "  " << Q[Y_] << "  " << Q[Z_] << "  "
       << Q[NB_] << endl;
  cout << "transformRF::Error\n";

  cout<<"\n\n velocity is "<<sqrt(vx * vx + vy * vy + vz * vz) <<"\n\n"<<endl;
 }
 if (!(nb < 0. || nb >= 0.)) {
  cout << "transformRF::Error nb=#ind\n";
   return ;
 }
}


void cnvrt::LRF_2_CALC(double eps, double p, double nb, double nq, double ns, double vx,
                 double vy, double vz, double* Q)  //converts energy density, velocity and number density to E,Mx,My
{
  double gamma2 = 1.0/(1.0- (vx*vx+vy*vy+vz*vz));
  double gamma = TMath::Sqrt(gamma2);
  Q[T_] = ((eps + p)*gamma2) - p;
  Q[X_] = (eps+p)*gamma2*vx;
  Q[Y_] = (eps+p)*gamma2*vy;
  Q[Z_] = (eps+p)*gamma2*vz;
  Q[NB_] = nb*gamma;
  Q[NQ_] = nq*gamma;
  Q[NS_] = ns*gamma;
}


/*
void CALC_2_LRF(EoS *eos, double* Q, double &eps, double &prs, double &nb,
                 double &nq, double &ns, double &vx, double &vy, double &vz) 
{
   double E = Q[T_];
   double Mx = Q[X_];
   double My = Q[Y_];
   double Mz = Q[Z_];
   double R_nb = Q[NB_];
   double R_nq = Q[NQ_];
   double R_ns = Q[NS_];
   double corrf = 0.9999;

   double M=TMath::Sqrt(Mx*Mx+My*My+Mz*Mz);

   if (E <= 0.)
    {
     eps = 0.;    // energy density
     vx  = 0.;    // velocity
     vy =0;
     vz =0;
     nb = 0.;    
     nq = 0.;    
     ns = 0.;    
     prs = 0.;
     return;
    }

   if(M == 0.)
    {
     eps = E;
     vx = 0.;
     vy =0.;
     vz =0.;
     nb =R_nb;
     nq =R_nq;
     ns =R_ns;
     prs = eos->pressure(eps,nb,nq,ns);
     return;
     }

     if(M>E){Mx *= corrf*E/M;  My *= corrf*E/M; Mz *= corrf*E/M; M = E*corrf;}

     double x0; 
     double x1;
     double epsl;
     double vl;
     
     x0=0.5;
     int _count = 0;
     do
       {
         eps = E-M*x0;
         nb = R_nb/(1.0/(1.-(x0*x0)));
         nq = R_nq/(1.0/(1.-(x0*x0)));
         ns = R_ns/(1.0/(1.-(x0*x0)));
         prs = eos->pressure(eps,nb,nq,ns);
         double f = (E+prs)*x0-M;
         double df = (E+prs)-x0*M*(1./3.);  //  -->dp/de = 1/3 ;
	 x1 = x0-(f/df);
	 epsl=abs(x1-x0);
	 x0=x1;
         _count ++;
         if(_count > 1000){cout<<"maximum iteration step(1000) to find out vl exceeds, exiting ..."<<endl; exit(1); }
       }
     while(epsl>0.01);

     vl=x0;
     eps= E - (M*vl);   // energy density
     vx=vl*(Mx/M);   // velocity
     vy=vl*(My/M);   // velocity
     vz=vl*(Mz/M);   // velocity
     nb= R_nb* sqrt(1 - vx * vx - vy * vy - vz * vz);
     nq= R_nq* sqrt(1 - vx * vx - vy * vy - vz * vz);
     ns= R_ns* sqrt(1 - vx * vx - vy * vy - vz * vz);
     prs = eos->pressure(eps,nb,nq,ns);


    if( vx >= 1.0 || vy >= 1.0 || vz >= 1.0  ){cout<<"[Info] Error in calculating velocity  : "<<"vx = "<<vx<<", vy = "<<vy<<endl;  exit(1);}
    if( eps<0.  ){ cout<<"[Info] Error in calculating energy density : Energy density = "<<eps<<endl;  eps = 0; cout<<"--> Setting, energy density = 0"<<endl; }
}

*/


