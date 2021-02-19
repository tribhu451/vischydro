#pragma once
#include <iostream>
#include "TMath.h"
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "eos.h"

using std::istringstream; 
using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::ios;


using namespace std;

class EoS1 : public EoS 
{



public:
EoS1();
~EoS1();
void check_eos();
double pressure( double eg,double _nb, double _nq, double _ns);
double entropy( double eg,double _nb, double _nq, double _ns);
double temperature( double eg,double _nb, double _nq, double _ns);
double cs(double eg,double _nb, double _nq, double _ns){return TMath::Sqrt(cs2(eg, _nb,  _nq,  _ns));};
double cs2(double eg,double _nb, double _nq, double _ns);
double cs2_(){return 1./3.;};
double cs_(){return sqrt(1./3.) ;};
double temp_2_eps( double T,double _nb, double _nq, double _ns)  ;
double temp_2_prs( double T,double _nb, double _nq, double _ns)  ;
double temp_2_entr( double T,double _nb, double _nq, double _ns)  ;

private:
  void is_file_exist(fstream& file);

  fstream infile1_d;
  fstream infile2_d;
  fstream infile3_d;
  fstream infile4_d;
  fstream infile5_d;
  fstream infile6_d;
  fstream infile7_d;

  fstream infile1_t;
  fstream infile2_t;
  fstream infile3_t;
  fstream infile4_t;
  fstream infile5_t;
  fstream infile6_t;
  fstream infile7_t;

  fstream file_teosdens; // for Teosdens

 // energy density starting value (e0_n)  ||
 // intervals in enrg. dens. value (de_n) ||
 // number of energy dens. points (ne_n)  ||
  double e0_1,de_1; int ne_1;
  double e0_2,de_2; int ne_2;
  double e0_3,de_3; int ne_3;
  double e0_4,de_4; int ne_4;
  double e0_5,de_5; int ne_5;
  double e0_6,de_6; int ne_6;
  double e0_7,de_7; int ne_7;

  // arrays to store data from file-1
  double eps_1[501]={0.0};
  double prs_1[501]={0.0};
  double ent_1[501]={0.0}; 
  double tmp_1[501]={0.0};

  // arrays to store data from file-2
  double eps_2[501]={0.0};
  double prs_2[501]={0.0};
  double ent_2[501]={0.0}; 
  double tmp_2[501]={0.0};

  // arrays to store data from file-3
  double eps_3[501]={0.0};
  double prs_3[501]={0.0};
  double ent_3[501]={0.0}; 
  double tmp_3[501]={0.0};

  // arrays to store data from file-4
  double eps_4[501]={0.0};
  double prs_4[501]={0.0};
  double ent_4[501]={0.0}; 
  double tmp_4[501]={0.0};

  // arrays to store data from file-5
  double eps_5[501]={0.0};
  double prs_5[501]={0.0};
  double ent_5[501]={0.0}; 
  double tmp_5[501]={0.0};

  // arrays to store data from file-6
  double eps_6[501]={0.0};
  double prs_6[501]={0.0};
  double ent_6[501]={0.0}; 
  double tmp_6[501]={0.0};

  // arrays to store data from file-7
  double eps_7[501]={0.0};
  double prs_7[501]={0.0};
  double ent_7[501]={0.0}; 
  double tmp_7[501]={0.0};


  double teos_e[750]={0.0}; // array to store data of
  double teos_p[750]={0.0}; // Teosdens file
  double teos_s[750]={0.0}; 
  double teos_t[750]={0.0};

  double nb,mub,mus,muq;

  istringstream* iss;
  char   buff[200];


  gsl_interp_accel *acc;

// temperature array pointers
gsl_spline *spline_temp_1; 
gsl_spline *spline_temp_2; 
gsl_spline *spline_temp_3; 
gsl_spline *spline_temp_4;
gsl_spline *spline_temp_5; 
gsl_spline *spline_temp_6;
gsl_spline *spline_temp_7; 

// pressure array pointers
gsl_spline *spline_prs_1;
gsl_spline *spline_prs_2;
gsl_spline *spline_prs_3;
gsl_spline *spline_prs_4;
gsl_spline *spline_prs_5;
gsl_spline *spline_prs_6;
gsl_spline *spline_prs_7;

// entropy array pointers
gsl_spline *spline_ntrpy_1; 
gsl_spline *spline_ntrpy_2; 
gsl_spline *spline_ntrpy_3; 
gsl_spline *spline_ntrpy_4;
gsl_spline *spline_ntrpy_5; 
gsl_spline *spline_ntrpy_6;
gsl_spline *spline_ntrpy_7; 


// spline for Teosdens
gsl_spline *spline_teos_e; 
gsl_spline *spline_teos_p; 
gsl_spline *spline_teos_s; 
gsl_spline *spline_teos_t; 

};








