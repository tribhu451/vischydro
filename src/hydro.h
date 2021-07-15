/////////////////////////////////////////////////////////
//     The viscous part is not                         //
//     following the tilde notation                    //
//     as written in vHLLE.                            //
//     means pi^{\eta \eta} != tau*tau*pi^{\eta \eta}  //
/////////////////////////////////////////////////////////


#pragma once

#include <cmath>
#include "TMath.h"
#include "cnvrt.h"
#include "grid.h"
#include "cell.h"
#include "TMath.h"
#include "global.h"
#include "trancoeff.h"
#include <fstream>
#include "hrr.h"


using std::cout;
using std::endl;

class hydro
{

public:
hydro(EoS* , grid* ,idb* ,double , double ,cnvrt* , trancoeff* );
~hydro();
void evolve();
void hlle_flux(cell* left, cell* right, int direction, int mode);
void sourcestep(int mode, int ix, int iy,int iz, double _tau);
inline double get_tau() { return tau; }
void set_dtau(double deltaTau);  // change the timestep
inline double get_dtau() { return dt; }  // return current value of timestep
void visc_flux(cell *left, cell *right, int direction);
void visc_source_step(int ix, int iy, int iz) ;
void setNSvalues();

private:
grid *f;
EoS *eos;
trancoeff* trcoef;
double dt;
double tau;
idb* IDB;
cnvrt* CN;

       //// viscous part ////

void NSquant(int ix,int iy,int iz, double dmunu[4][4], double &PiNS, double piNS[4][4] );
void ISformal();
hrr* hr; // hydro regulation

inline double Gamma(int a, int b, int c, double _tau) // \Gamma^{a}_{bc}
{
 if(a==3 && b == 0 && c == 3){return 1.0/_tau;}
 else if(a==3 && b == 3 && c == 0){return 1.0/_tau;}
 else if(a==0 && b == 3 && c == 3){return _tau;}
 else {return 0;}
}


double sign(double x) {
 if (x > 0)
  return 1.;
 else if (x < 0.)
  return -1.;
 else
  return 0.;
}


};

