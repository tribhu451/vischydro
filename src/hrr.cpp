#include "hrr.h"
#include "TMath.h"


hrr::hrr()
{
};

hrr::~hrr()
{
};

void hrr::hrr_vhlle(double pi[4][4], double Pi,
     double quant[6],double pi_update[4][4], double &Pi_update, bool &rescaled)
{

//double tau = quant[0];
double eps = quant[1];
double prs = quant[2];
double vx  = quant[3];
double vy  = quant[4];
double vz  = quant[5];

	  double maxT0 = max((eps + prs) / (1. - vx * vx - vy * vy - vz * vz) - prs,
			     (eps + prs) * (vx * vx + vy * vy + vz * vz) /
			     (1. - vx * vx - vy * vy - vz * vz) +
			     prs);
	  // double maxpi = max(fabs(pi[1][1]),fabs(pi[2][2])) ;
	  double maxpi = 0.;
	  for (int i = 0; i < 4; i++)
	    for (int j = 0; j < 4; j++)
	      if (fabs(pi[i][j]) > maxpi) maxpi = fabs(pi[i][j]);
	   rescaled = false;
	  if (maxT0 / maxpi < 1.0)
	    {
	      for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		  {
		    pi_update[i][j] = 0.1 * pi[i][j] * maxT0 / maxpi;
		  }
	      rescaled = true;
	    }
      else
    {

     for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        pi_update[i][j] = pi[i][j] ;

    }


	  if (fabs(Pi) > prs)
	    {
	      if (Pi != 0.) Pi_update = 0.1 * Pi / fabs(Pi) * prs;
	      rescaled = true;
	    }
          else 
          {
              Pi_update = Pi;
          }

}



void hrr::hrr_music(double pi[4][4], double Pi,
     double quant[6],double pi_update[4][4], double &Pi_update, bool &rescaled)
{
double eps_warning = 0.015; //~100MeV

double tau = quant[0];
double eps = quant[1];
double prs = quant[2];
double vx  = quant[3];
double vy  = quant[4];
double vz  = quant[5];

double u[4] = {0.,0.,0.,0.};
 u[0] = 1.0 / sqrt ( 1.0 - vx*vx - vy*vy - vz*vz );
 u[1] = u[0]*vx;
 u[2] = u[0]*vy;
 u[3] = u[0]*vz;

double gmunu[4] = {1.,-1.,-1.,-1./(tau*tau)};


double pimunu2 = 0.0; // pi^{\mu \nu} * \pi_{\mu \nu}
for(int i=0; i<4; i++) 
 for(int j=0; j<4; j++)
  pimunu2 += pi[i][j] * ( (1.0/gmunu[i]) * (1.0/gmunu[j]) * pi[i][j] ) ;
 
double ideal = eps*eps + 3.0*prs*prs  ;

 double rho_pi = 1.0 ;
 double rho_Pi = 1.0 ;
 double factor = 300.0*tanh(eps);

 pimunu2 = TMath::Abs(pimunu2)+1e-15;
 ideal = TMath::Abs(ideal)+1e-15;

 rho_pi = 1/factor * sqrt(pimunu2/ideal);
 rho_Pi = 1/factor * sqrt( (3*Pi*Pi) / ideal );
   if(std::isinf(rho_pi) or std::isnan(rho_pi) )
    {
     std::cout<<"rho pi  = inf/nan [wrong hrr]"<<std::endl;
     std::cout<<"pimunu2 = "<<pimunu2<<std::endl;
     std::cout<<"ideal = "<<ideal<<std::endl;
     exit(1);
    }

 if(rho_pi > 0.15){

   if(eps > eps_warning)std::cout<<"[Warning] hrr pi required at eps = "<<eps<<std::endl;

   for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++){
       pi_update[i][j] = 0.1/rho_pi * pi[i][j];
       }

       rescaled = true;
  } 
  else
    {
   for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++)
       pi_update[i][j] = pi[i][j];

    }


 if(rho_Pi > 0.15)
   {
    if(eps > eps_warning)std::cout<<"[Warning] hrr PI required at eps = "<<eps<<std::endl;
    Pi_update = 0.1/rho_Pi * Pi;

     rescaled = true;
  } 
 else
   {
    Pi_update = Pi ;

   }

 for(int i = 0; i<4 ; i++)
  for(int j=0; j<4 ; j++){
   if(std::isinf(pi[i][j]) or std::isnan(pi[i][j]) )
    {
     std::cout<<"1.0/rho_pi = " << 1.0/rho_pi << std::endl;
     std::cout<<"pi("<<i<<","<<j<<") = inf/nan [wrong hrr]"<<std::endl;
    }
 } 

}



void hrr::hrr_music2(double pi[4][4], double Pi,
     double quant[6],double pi_update[4][4], double &Pi_update, bool &rescaled)
{
double eps_warning = 0.015; //~100MeV
double xi = 0.001;

double tau = quant[0];
double eps = quant[1];
double prs = quant[2];
//double vx  = quant[3];
//double vy  = quant[4];
//double vz  = quant[5];

double gmunu[4] = {1.,-1.,-1.,-1./(tau*tau)};


double pimunu2 = 0.0; // pi^{\mu \nu} * \pi_{\mu \nu}
for(int i=0; i<4; i++) 
 for(int j=0; j<4; j++)
  pimunu2 += pi[i][j] * ( (1.0/gmunu[i]) * (1.0/gmunu[j]) * pi[i][j] ) ;
 
double ideal = eps*eps + 3.0*prs*prs  ;

 double rho_pi = 1.0 ;
 double rho_Pi = 1.0 ;
 double factor = 100.*(1./(exp(-(eps - eps_warning)/xi) + 1.)
                          - 1./(exp(eps_warning/xi) + 1.));

 pimunu2 = TMath::Abs(pimunu2)+1e-15;
 ideal = TMath::Abs(ideal)+1e-15;

 rho_pi = 1/factor * sqrt(pimunu2/ideal);
 rho_Pi = 1/factor * sqrt( (3*Pi*Pi) / ideal );

   if(std::isinf(rho_pi) or std::isnan(rho_pi) )
    {
     std::cout<<"rho pi  = inf/nan [wrong hrr]"<<std::endl;
     std::cout<<"pimunu2 = "<<pimunu2<<std::endl;
     std::cout<<"ideal = "<<ideal<<std::endl;
     exit(1);
    }

 if(rho_pi > 0.1){

   if(eps > eps_warning)std::cout<<"[Warning] hrr pi required at eps = "<<eps<<std::endl;

   for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++){
       pi_update[i][j] = 0.1/rho_pi * pi[i][j];
       }

       rescaled = true;
  } 
  else
    {
   for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++)
       pi_update[i][j] = pi[i][j];

    }


 if(rho_Pi > 0.1)
   {
    if(eps > eps_warning)std::cout<<"[Warning] hrr PI required at eps = "<<eps<<std::endl;
    Pi_update = 0.1/rho_Pi * Pi;

     rescaled = true;
  } 
 else
   {
    Pi_update = Pi ;

   }

 for(int i = 0; i<4 ; i++)
  for(int j=0; j<4 ; j++){
   if(std::isinf(pi[i][j]) or std::isnan(pi[i][j]) )
    {
     std::cout<<"1.0/rho_pi = " << 1.0/rho_pi << std::endl;
     std::cout<<"pi("<<i<<","<<j<<") = inf/nan [wrong hrr]"<<std::endl;
    }
 } 

}







