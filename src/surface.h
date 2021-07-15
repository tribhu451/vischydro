#pragma once
#include <iostream>
#include <iomanip>      // std::setprecision
#include <cmath>
#include "cnvrt.h"
#include "grid.h"
#include "idb.h"
#include "cornelius.h"
#include "eos.h"

class surf
{

public:

surf(EoS* , idb* , grid* , cnvrt*  );
~surf();

void get_surface2p1d(double );
void get_surface3p1d(double );

idb* IDB;
grid* f;
Cornelius *cornelius;  // instance of Cornelius to calculate the hypersurface
cnvrt* CN;
EoS* eos;
std::fstream freeze_file;

private:

};
