#pragma once
#include<iostream>
#include "TMath.h"
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
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

class EoS0 : public EoS 
{

public:
EoS0();
~EoS0();

double pressure( double eg,double _nb, double _nq, double _ns);
double entropy( double eg,double _nb, double _nq, double _ns);
double temperature( double eg,double _nb, double _nq, double _ns);
double cs_();
double cs2_();
double cs2(double eg,double _nb, double _nq, double _ns);
double cs(double eg,double _nb, double _nq, double _ns);
double temp_2_eps( double T,double _nb, double _nq, double _ns)  ;

};
