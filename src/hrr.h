
// hydro regulation routine

#pragma once
#include <iostream>
#include <cmath>

class hrr
{

private :

double rho_max = 5.0;
double zeta0 = 0.1;
inline double max(double a, double b)
 {
   if(a>b) {return a;}
   else    {return b;}
 }


public :
hrr();
~hrr();

void hrr_vhlle(double pi[4][4], double Pi, 
      double quant[6],double pi_update[4][4], double &Pi_update, bool &rescaled);

void hrr_music(double pi[4][4], double Pi, 
      double quant[6],double pi_update[4][4], double &Pi_update, bool &rescaled);

void hrr_music2(double pi[4][4], double Pi, 
      double quant[6],double pi_update[4][4], double &Pi_update, bool &rescaled);


//void check_it(double pi[4][4],double u[5]);

};


