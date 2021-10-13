#pragma once
#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "idb.h"
#include "grid.h"
#include "eos.h"
#include "cell.h"
#include "opt_glau.h"
#include "mc_glau.h"
#include "fluid_info.h"
#include "read_ic.h"
#include "hydro.h"
#include "eos0.h"
#include "eos2.h"
#include "eos1.h"
#include "surface.h"
#include "freeze.h"
#include "trancoeff.h"


using std::cout;
using std::endl;
using namespace std::chrono;

class master {

public:

  master(idb* );
  ~master();
  void init();
  void run_hydro();
  EoS* eos;
  idb* IDB;
  grid* g;
  fluid_info* out;
  hydro* h;
  cnvrt* CN;
  surf* sf;
trancoeff* trcoef;


private:

int check_to_stop(grid* , EoS* , double  );
void save_for_fo(int );
void find_hyper_surface(int , double);



} ;
