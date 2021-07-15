#pragma once
#include<iostream>
#include<fstream>
#include "cell.h"
#include "global.h"
#include "eos.h"
#include "trancoeff.h"
#include "cornelius.h"
#include "input_data.h"




using namespace std;

class fluid
{
  
 private:
  EoS* eos ;
  cell *Cell;            
  int nx, ny,nz;
  double xmin, xmax, ymin, ymax,zmin, zmax;
  double dx, dy,dz, dt; // dt <=> d_tau 
  Cornelius *cornelius;  // instance of Cornelius to calculate the hypersurface
  std::ofstream freeze_file;
  InputData* DATA;

 public:

  fluid(EoS* _eps,  trancoeff* _tr, InputData* );
  ~fluid();
  trancoeff* trcoef;
  double get_nx(){return nx;}  
  double get_ny(){return ny;}
  double get_nz(){return nz;}
  inline double get_xmax() { return xmax; }
  inline double get_ymax() { return ymax; }
  inline double get_zmax() { return zmax; }
  inline double get_xmin() { return xmin; }
  inline double get_ymin() { return ymin; }
  inline double get_zmin() { return zmin; }
  double get_dx(){return dx; } 
  double get_dy(){return dy;}
  double get_dz(){return dz;}
  double get_dt(){return dt;}
  double get_x(int ix) {return  xmin + ix * dx; }      //returns x-coordinate of a cell
  double get_y(int iy) {return  ymin + iy * dy; }      //returns y-coordinate of a cell
  double get_z(int iz) {return  zmin + iz * dz; }      //returns eta-coordinate of a cell
  void getCMFvariablesUmunu(cell *c, double tau, double &e, double &nb, 
			    double &nq, double &ns, double &vx, double &vy, double &Y);

  void correct_imaginary_cells();
  void correct_imaginary_cells_full();
  

 inline cell *get_cell(int ix, int iy, int iz) {
  ix = ix > 0 ? ix : 0; //reult = (condition) ? expressionTrue : expressionFalse;
  ix = ix < nx ? ix : nx - 1; //https://www.w3schools.com/cpp/cpp_conditions_shorthand.asp
  iy = iy > 0 ? iy : 0;
  iy = iy < ny ? iy : ny - 1;
  iz = iz > 0 ? iz : 0;
  iz = iz < nz ? iz : nz - 1;
  return &Cell[ix + nx * iy + nx * ny * iz];
 }


void get_surface3d(double tau); // freezeout hypersurface;
void get_surface2d(double tau); // freezeout hypersurface;
  
};
