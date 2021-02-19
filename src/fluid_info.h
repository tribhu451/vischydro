#pragma once
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <cmath>
#include "grid.h"
#include "eos.h"
#include "idb.h"
#include <iomanip>
#include <string>

using namespace std;

class fluid_info{

  public :
  fluid_info(grid* , EoS* , idb* );
  ~fluid_info();
  void anisotropy_in_xy_plane(double time);
  void midslice_xy_dist(double time);
  void conservation_check(double time);

 private :
 grid *f;
 EoS* eos;
 idb *IDB;
 ofstream File1;

};
