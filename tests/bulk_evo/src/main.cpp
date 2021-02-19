#include <iostream>
#include <fstream>
#include <cmath>
#include<gsl/gsl_sf.h>

using std::cout;
using std::endl;


int main()
{

std::ofstream File0;
std::ofstream File1;

File0.open("analytical.dat");

double Pi0 = 0.0;
double zeta = 1000.0; //GeV/fm^2
double tau0 = 1.0; //fm
double tauPi = 1.0; //fm

double Pi;

for(double tau=tau0; tau<10; tau=tau+0.01){
    Pi = Pi0*exp((tau0-tau)/tauPi)  +  
         (zeta/tauPi) * exp(-tau/tauPi) *
           ( gsl_sf_expint_Ei(tau0/tauPi) - gsl_sf_expint_Ei(tau/tauPi) ) ;
  File0<<tau<<"\t"<<Pi<<endl;
   }

File0.close();



File1.open("numerical.dat");

Pi = Pi0; //Pi^{n}
double tau =tau0;
double dt = 0.00100;
int nstep = (10-tau0)/dt ;

for(int istep = 0; istep<nstep; istep++)
{
Pi = Pi0 - (dt/tauPi)*(Pi0 - zeta/tau) ;
// Pi -= (4.0*Pi0)/(3.0*tau)*dt ;   //this will be added for additional term \frac{4 \Pi} {3 \tau}
tau += dt;
File1<<tau<<"\t"<<-Pi<<endl;
Pi0 = Pi ;
}

File1.close();

return 0;
}


