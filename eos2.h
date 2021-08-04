#pragma once
#include <iostream>
#include "TMath.h"
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include "eos.h"

using std::istringstream; 
using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::ios;


using namespace std;

class EoS2 : public EoS 
{
  
 public:
  EoS2();
  ~EoS2();
  void check_eos();
  double pressure( double e,double _nb, double _nq, double _ns);
  double entropy( double ,double , double , double );
  double temperature( double ,double , double , double );
  double mu_b( double ,double , double , double );
  double cs2(double ,double , double , double );
  double cs(double e,double _nb, double _nq, double _ns)
            {return sqrt(cs2(e, _nb,  _nq,  _ns));};


  double temp_2_eps( double T,double _nb, double _nq, double _ns)  ;
  double entr_2_eps(double e,double _nb, double _nq, double _ns)  ;
  
 private:
  void is_file_exist(fstream& file);
    
  double ***pressure_tb;
  double ***temperature_tb;
  double ***mu_B_tb;
  
  double nb_bounds[7];
  double e_bounds[7];
  double nb_spacing[7];
  double e_spacing[7];
  double nb_length[7];
  double e_length[7];
  
  
  double **mtx_malloc(const int n1, const int n2)
  {
    double **d1_ptr; 
    d1_ptr = new double *[n1];
    
    for(int i=0; i<n1; i++) 
      d1_ptr[i] = new double [n2];
    
    for(int i=0; i<n1; i++) 
      for(int j=0; j<n2; j++) 
        d1_ptr[i][j] = 0.0;
    
    return d1_ptr;
  }
  
  
  
  void mtx_free(double **m, const int n1, const int n2)
  {
    for (int j = 0; j < n1; j++) 
      delete [] m[j];
    delete [] m;
  }
  
  
  
  int get_table_idx(double e) const
  {
    //double local_ed = e*hbarc;  // [GeV/fm^3]
    double local_ed = e;  // [GeV/fm^3]
    for (int itable = 1; itable < 7; itable++)
      {
	if (local_ed < e_bounds[itable])
	  {
	    return(itable - 1);
	  }
      }
    return(std::max(0, 6));
  }
  
  
  
  double interpolate2D(double e, double rhob, int table_idx, double ***table) const
  {
    // This is a generic bilinear interpolation routine for EOS at finite mu_B
    // it assumes the class has already read in
    //        P(e, rho_b), T(e, rho_b), s(e, rho_b), mu_b(e, rho_b)
    // as two-dimensional arrays on an equally spacing lattice grid
    // units: e is in GeV/fm^3, rhob is in 1/fm^3
    //double local_ed = e*hbarc;  // [GeV/fm^3]
    double local_ed = e;
    double local_nb = rhob;     // [1/fm^3]
    
    double e0       = e_bounds[table_idx];
    double nb0      = nb_bounds[table_idx];
    double delta_e  = e_spacing[table_idx];
    double delta_nb = nb_spacing[table_idx];
    
    int N_e  = e_length[table_idx];
    int N_nb = nb_length[table_idx];
    
    // compute the indices
    int idx_e  = static_cast<int>((local_ed - e0)/delta_e);
    int idx_nb = static_cast<int>((local_nb - nb0)/delta_nb);
    
    // treatment for overflow, use the last two points to do extrapolation
    idx_e  = std::min(N_e - 2, idx_e);
    idx_nb = std::min(N_nb - 2, idx_nb);
    
    // check underflow
    idx_e  = std::max(0, idx_e);
    idx_nb = std::max(0, idx_nb);
    
    double frac_e    = (local_ed - (idx_e*delta_e + e0))/delta_e;
    double frac_rhob = (local_nb - (idx_nb*delta_nb + nb0))/delta_nb;
    
    double result;
    double temp1 = table[table_idx][idx_nb][idx_e];
    double temp2 = table[table_idx][idx_nb][idx_e + 1];
    double temp3 = table[table_idx][idx_nb + 1][idx_e + 1];
    double temp4 = table[table_idx][idx_nb + 1][idx_e];
    result = ((temp1*(1. - frac_e) + temp2*frac_e)*(1. - frac_rhob)
              + (temp3*frac_e + temp4*(1. - frac_e))*frac_rhob);
    return(result);
  }
  
  


  
};








