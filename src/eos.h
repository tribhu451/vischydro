#pragma once    //The #pragma once directive has a very simple concept.
                //The header file containing this directive is included- 
                //only once even if the programmer includes it multiple times during a compilation.
                // https://www.geeksforgeeks.org/pragma-directive-in-c-c/
#include<iostream>
#include "TMath.h"
#include<cmath>

using std::cin;
using std::cout;
using std::endl;
using std::fstream;

class EoS
{

public:
virtual ~EoS(){};   //A virtual function is a member function which is declared-
                    // within a base class and is re-defined(Overriden) by a derived class.

                    //https://www.geeksforgeeks.org/virtual-function-cpp/?ref=lbp
                    //https://www.geeksforgeeks.org/virtual-functions-and-runtime-polymorphism-in-c-set-1-introduction/?ref=lbp
virtual double pressure( double eg,double _nb, double _nq, double _ns){return 0;};
virtual double entropy( double eg,double _nb, double _nq, double _ns){return 0;};
virtual double temperature( double eg,double _nb, double _nq, double _ns){return 0;};
virtual double temp_2_eps( double T,double _nb, double _nq, double _ns) {return 0;} ;
virtual double temp_2_prs( double T,double _nb, double _nq, double _ns) {return 0;} ;
virtual double temp_2_entr( double T,double _nb, double _nq, double _ns) {return 0;} ;
virtual double cs(double eg,double _nb, double _nq, double _ns){return 0;};
virtual double cs2( double eg,double _nb, double _nq, double _ns){return 0;};
virtual double cs_(){return 0;};
virtual double cs2_(){return 0;};
virtual void check_eos(){};
};
