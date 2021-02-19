///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                       //
//   conversion from E,M,R (calculational frame variables) to energy density, vx (Local rest frame variables) ... etc    //
//                                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "eos.h"

class cnvrt {

public:
cnvrt();
~cnvrt();

void CALC_2_LRF(EoS *eos, double* Q, double &eps, double &prs, double &nb,
                 double &nq, double &ns, double &vx, double &vy, double &vz) ;// converts E,Mx,My,R to energy density,velocity and number density

void LRF_2_CALC(double e, double p, double nb, double nq, double ns, double vx,
                 double vy, double vz, double* Q) ; //converts energy density, velocity and number density to E,Mx,My

} ;
