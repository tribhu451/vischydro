#pragma once

#include <iostream>
#include <iomanip>
#include "eos.h"
#include "global.h"
#include "idb.h"

// this class contains the information about the transport coefficients
// of the fluid: eta/s, zeta/s and the corresponding relaxation times,
// taupi (\tau_\pi) and tauPi (\tau_\Pi)
class trancoeff {
  double etaS,  taupi, tauPi;

  idb* InData;  

  // flags
  int zetas_flag;
  int etas_flag;
  int t_etas_flag;
  
  EoS *eos;  // EoS instance is needed optionally for zeta/s parametrization,
  // which depends on the speed of sound
  
  double zetaS(double e, double T);
  double EtaS(double e, double T);
  
 public:
  
  trancoeff( EoS *_eos, idb*);
  ~trancoeff();
  
  // returns (optionally temperature dependent) eta/s and zeta/s
  void getEta(double e, double T, double &_etaS, double &_zetaS);
  
  // returns shear and bulk relaxation times
  void getTau(double e, double T, double &_taupi, double &_tauPi);
  
  // deltapipi, taupipi, lambdapiPi * divided by tau_pi * !
  void getOther(double e, double nb, double nq, double ns, 
		double &deltapipi, double &taupipi, double &lambdapiPi, double &phi7)
  {
    deltapipi = 4./3.;  taupipi = 10./7.;  lambdapiPi = 6./5.;
    phi7 = 9./70./eos->pressure(e, nb, nq, ns);
   
    //deltapipi = 4./3.;  taupipi = 0.;  lambdapiPi = 0.;
    //phi7 = 0.;
    
    if(std::isinf(phi7)) phi7=0.0;
    
  }
  
  
  void getOtherBulk(double e, double nb, double nq, double ns, 
		    double &delPiPi, double &lamPipi) {
    
    if(zetas_flag > 0)
      {
	delPiPi = 2./3.;  lamPipi = 8./5.*(1. / 3. - eos->cs2(e,nb,nq,ns));
      }
    else
      {
	delPiPi = 0. ;  lamPipi = 0.;
      }
    
  }

  
  // isViscous tells whether the fluid is viscous or inviscid
  inline bool isViscous()
  {
    if (etas_flag > 0 || zetas_flag > 0)
      return true;
    else
      return false;
  }

  
};
