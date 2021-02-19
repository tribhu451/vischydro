
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "Vector3D.h"
#include "sVector3D.h"


class par{
 public:

 

  static constexpr double etaadd= .5; // range in eta = etabeam+etaadd
  static constexpr double dtfixed=0.025;
  static constexpr double dtsave=0.2;
  static constexpr double initime=0.25;



  //  static const double thetaflux=1.1; // flux control parameter (obsolete)
  //  static const double damp=0.005;
  static const int ensolve_iterations=15;
  static const int distiter=18;
  static constexpr double Tmintau=0.08;
  static constexpr double Tsc=0.008;
 
 
  static constexpr double rmin=0;
  static const int ntheta=91;
  static const int nzeta=91;
  static const int nphi=257;
  //  static const int ntheta=21;
  //  static const int nzeta=21;
  //   static const int nphi=17;
 

  static constexpr double nucleonmass=0.939;
  static constexpr double ren=0.19732;   // Gev*fm
  static constexpr double Pi=3.1415926535897932384626433;

 
};

class evolve{
 public:
  double time,dtau,dtsave;
  double etavis(double en, EoS* EOS,grid* grid);
  void ini(grid* grid,double tau0,double maxtime,double dtauin);
  double zetavis(double en,EoS* EOS, grid* grid);
  double endens(double d,double phi, double xi, double theta);
  double value(double t,double x, double y, double eta, 
	       sVector3D* tab[]);
  double dist(double enF,double phi, double xi, double theta);
  void hypersurface(const char* filename,int eventcount);
  void gethypersurface(grid* grid, EoS* EOS,double Tf);
  void put(int step,grid* grid,double tau,EoS* EOS);
  int stepsave,maxstep,timesaved;
  double tmin,lambda,Tfreeze,cs2freeze,etafreeze;
  double  maxtimesaved,distmax;
  double maxx,maxy,maxeta;
  int nnx,nny,nneta;
  int Nw,eventcount;
  double rds;
  double etaovs,zetaovs;
void writeVector3D(std::ofstream& oh,const char* name,const char* unit,
			   Vector3D* hyp);
  void setAuAu200ev(int cent,int &nmin, int &nmax);
  double ebegamma(double x,double y,double eta);
  double readeventquark(std::ifstream& eventfile,int ievent);
  void setPbPbquark(int cent,double &nmin, double &nmax,EoS *EOS,double t0);
  void setebe(grid *f, EoS *eos, double tau0) ;
  void   writefinal(grid *f,EoS *eos,double eta,double t0,char* dirname);
  void   writespacetimehistory(grid *f,EoS *eos,double t0,char* dirname);
 char system[10];
  char device[10];
  double sqrts,b,cmin,cmax;
  double muB,muS,muI,muC;
  double inieps2_2,inieps3_3,psi3;
  double T0center,s0center;
  Vector3D *tux,*tuy,*tue,*ten;
  Vector3D *tuxo,*tuyo,*tueo;
  sVector3D**  energymap;
  sVector3D**  uxmap;
  sVector3D**  uymap;
  sVector3D**  uemap;
  sVector3D**  pixxmap;
  sVector3D**  pixymap;
  sVector3D**  pixemap;
  sVector3D**  piyymap;
  sVector3D**  piyemap;
  sVector3D**  pieemap;
  sVector3D**  bulkmap;
  Vector3D *disthyper,*uxhyper,*uyhyper,*Yhyper;
  Vector3D *pixxhyper,*pixyhyper,*pixehyper,*piyyhyper,*piyehyper,*pieehyper
    ,*bulkhyper;

 double cnw[2500];
  double gw[2500];
  double cnwsum[50000];
  double xnwsum[50000],ynwsum[50000];//,cnw[2500]; // wounded nucleons
  double xnw[2500],ynw[2500];//,cnw[2500]; // wounded nucleons
 private:
   double RDSwaga;
 double sig0,etas,x0,sigma,alpha,rho0;

 
};
