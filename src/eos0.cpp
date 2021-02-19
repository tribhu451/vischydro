#include "eos0.h"


EoS0::EoS0(){    cout<<"\n[Info] allocating ideal EoS...\n"<<endl;}
EoS0::~EoS0(){}
double EoS0::pressure( double eps, double _nb, double _nq, double _ns)
                 {return eps/3.0;} // in GeV/fm^{3}


double EoS0::entropy( double eps,double _nb, double _nq, double _ns) // entropy in fm^{-3}
                {  
  if (eps == 0.0) {return 0;}
  else{
return  ( 4.0 * eps ) / ( 3.0 * temperature(eps ,_nb ,_nq ,_ns) ); 
}
}


double EoS0::temperature( double eps,double _nb, double _nq, double _ns) // T in GeV
{ 	

  if (eps == 0.0) {return 0;}
  else{
  double gg = 16.0; // # gluon degeneracy
  double nc = 3.0 ; double nf = 2.5 ; // # colours and flavours 
  double aux = 3 * ( gg + 7./2. * nc * nf ) * 3.1415927 * 3.1415927 /90. ;
  return pow( ( eps * 0.007645 ) / aux , 0.25);  // page 163 of book by C. Y. Wong Eq. (9.2) - (9.7) .
}
}

double EoS0::cs(double eps ,double _nb, double _nq, double _ns){return 1./TMath::Sqrt(3);}

double EoS0::cs2(double eps ,double _nb, double _nq, double _ns){return 1./3.;}

double EoS0::cs_(){return 1./TMath::Sqrt(3);} // HLLE Ideal flux will call me .
 
double EoS0::cs2_(){return 1./3.;}  // HLLE Ideal flux will call me .

double EoS0::temp_2_eps( double T,double _nb, double _nq, double _ns)  
{
  double gg = 16.0; // # gluon degeneracy
  double nc = 3.0 ; double nf = 2.5 ; // # colours and flavours 
  double aux = 3 * ( gg + 7./2. * nc * nf ) * 3.1415927 * 3.1415927 /90. ;
  return aux/0.007645 * pow(T,4.0);
 
}

