// viscous hydro [ Tribhuban ] [ Feb. 01 2021 ]

// This class will create the entire grid
// will store adress of each cell in a vector
// and returns the adress of a particular cell.

#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "idb.h"
#include "cell.h"
#include <string>
#include <sstream>
#include <fstream>
#include "cnvrt.h"
#include "trancoeff.h"
#include "eos.h"

class grid {
  
 
 private :

  int nx,ny,neta;
  double xmin, xmax, ymin, ymax,etamin, etamax;
  double dx, dy,deta, dtau; // dt <=> d_tau 
  
  istringstream* iss;
  char   buff[200];


  


 public:
  grid(idb*, cnvrt*, trancoeff*, EoS* );
  ~grid();
  
  idb* IDB;
  cell* Cell;
  trancoeff* trcoef;
  EoS* eos;
  
  void make_grid();
  cell* new_cell(int, int, int );
  std::vector < vector < vector<cell> > > tube; // a cell type vector

  void correct_imaginary_cells() ;
  void correct_imaginary_cells_full() ;

  cnvrt* CN;

  inline double get_nx(){return nx;}  
  inline double get_ny(){return ny;}
  inline double get_neta(){return neta;}
  inline double get_xmax() { return xmax; }
  inline double get_ymax() { return ymax; }
  inline double get_etamax() { return etamax; }
  inline double get_xmin() { return xmin; }
  inline double get_ymin() { return ymin; }
  inline double get_etamin() { return etamin; }
  inline double get_dx(){return dx; } 
  inline double get_dy(){return dy;}
  inline double get_deta(){return deta;}
  inline double get_x(int ix) {return  xmin + ix * dx; }      //returns x-coordinate of a cell
  inline double get_y(int iy) {return  ymin + iy * dy; }      //returns y-coordinate of a cell
  inline double get_eta(int iz) {return  etamin + iz * deta; }      //returns eta-coordinate of a cell

  inline void set_dtau(int _dtau) {dtau = _dtau; }      
  inline double get_dtau() {return dtau; }      


inline cell* get_cell(int ix, int iy, int iz)
{
  // border cells next or prev are they themselves :)
  ix = ix > 0 ? ix : 0;
  ix = ix < IDB->nx ? ix : IDB->nx - 1;
  iy = iy > 0 ? iy : 0;
  iy = iy < IDB->ny ? iy : IDB->ny - 1;
  iz = iz > 0 ? iz : 0;
  iz = iz < IDB->neta ? iz : IDB->neta - 1;
  
  //return &tube[ix][iy][iz];
  return &Cell[ix + nx * iy + nx * ny * iz];
}


void getCMFvariablesUmunu(cell *c, double tau, double &e, double &nb, 
				 double &nq, double &ns, double &vx, double &vy, double &Y){
  double p; double vz;
  c->get_physical_var(eos, tau, e, p, nb, nq, ns,vx, vy, vz);
  double eta = get_eta(c->get_iz());
  Y = eta + 0.5*log((1.+vz)/(1.-vz));
  double gf=1./sqrt(1.-vx*vx-vy*vy-vz*vz);
  vx = vx*gf;
  vy = vy*gf;
}

 
  
} ;





