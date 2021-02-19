
#define FSIDEONE 0
#define RDSFLAGA 0

// #Ddefine FSIDEONE 0


#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include "eos.h"
#include "trancoeff.h"
#include "grid.h"
#include "freeze.h"
#include "cell.h"

using std::cout;
using std::endl;


double HeavisideTheta(double x){
  if (x>0){return 1;}
  else {return 0;}
}

inline double minimum(double a, double b){
  if (a>b) 
    return b;
  else
    return a;
}
inline double maximum(double a, double b){
  if (a>b) 
    return a;
  else
    return b;
}

void evolve::ini(grid* grid,double t0,double maxtime,double dtauin){
 
  maxx=grid->get_xmax(); maxy=grid->get_ymax(); maxeta=grid->get_etamax();
  nnx=(grid->get_nx()+1)/2; 
  nny=(grid->get_ny()+1)/2; 
  nneta=(grid->get_neta()+1)/2;
  tmin=t0;lambda=maxx/maxeta;
  dtau=dtauin;
  maxstep=(maxtime-t0)/dtau+1;
  stepsave=(int) (par::dtsave/dtau);
  dtsave=stepsave*dtau;
  timesaved=(maxstep/stepsave)+1;
  maxtimesaved=(timesaved-1)*dtsave;
  distmax=maximum(maxtimesaved,maxx);
  cout<<" time step "<<dtau<<endl;
  cout<<" save every "<< stepsave<<" step"<<endl;
  cout<<" calculating "<<maxstep<<" timesteps"<<endl;
  cout<<" saving "<<timesaved<<" timesteps"<<endl;
  cout<<" maxtimesaved "<<maxtimesaved<<endl;
  energymap= new sVector3D*[timesaved];
  uxmap= new sVector3D*[timesaved];
  uymap= new sVector3D*[timesaved];
  uemap= new sVector3D*[timesaved];
  pixxmap= new sVector3D*[timesaved];
  pixymap= new sVector3D*[timesaved];
  pixemap= new sVector3D*[timesaved];
  piyymap= new sVector3D*[timesaved];
  piyemap= new sVector3D*[timesaved];
  pieemap= new sVector3D*[timesaved];
  bulkmap= new sVector3D*[timesaved];
 
  for (int i=0;i<timesaved;i++){
    energymap[i]=new  sVector3D(" energymap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    uxmap[i]=new  sVector3D(" uxmap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    uymap[i]=new  sVector3D(" uymap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    uemap[i]=new  sVector3D(" uemap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    pixxmap[i]=new  sVector3D(" pixxmap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    pixymap[i]=new  sVector3D(" pixymap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    pixemap[i]=new  sVector3D(" pixemap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    piyymap[i]=new  sVector3D(" piyymap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    piyemap[i]=new  sVector3D(" piyemap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    pieemap[i]=new  sVector3D(" pieemap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    bulkmap[i]=new  sVector3D(" bulkmap ",-maxx,maxx,nnx
		   ,-maxy,maxy,nny
		   ,-maxeta,maxeta,nneta);
    
  }
  disthyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  uxhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  uyhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  Yhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  pixxhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  pixyhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  pixehyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  piyyhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  piyehyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  pieehyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
  bulkhyper=new Vector3D(" hypersurface radius ",0.0,par::Pi/2.,par::nzeta
			 ,0,2.*par::Pi,par::nphi
			 ,0,par::Pi,par::ntheta);
}


void evolve::gethypersurface(grid* grid, EoS* EOS,double Tf){
  Tfreeze=Tf;
  cout<<"calculating freeze-out hypersurface at Tf="<<Tfreeze<<endl;

  double Efreeze=EOS->temp_2_eps(Tfreeze,0,0,0);
  cout<<"calculating freeze-out hypersurface at Ef="<<Efreeze<<endl;
  double dzeta=par::Pi/2./(par::nzeta-1);
  double dphi=par::Pi*2./(par::nphi-1);
  double dtheta=par::Pi/(par::ntheta-1);
  etafreeze=etavis(EOS->temp_2_eps(Tf,0,0,0),EOS,grid)/EOS->entropy(Efreeze,0,0,0);
  grid->trcoef->getEta(EOS->temp_2_eps(Tf,0,0,0),Tf,etaovs,zetaovs);
  cout<<"calculating freeze-out hypersurface at etaf="<<etafreeze<<endl;
  cs2freeze=EOS->cs2(Efreeze,0,0,0);
  //cout<<EOS->DPE(Efreeze)<<endl;
  cout<<"calculating freeze-out hypersurface at cs2="<<cs2freeze<<endl;

  double entalpyfreeze=Efreeze+EOS->pressure(Efreeze,0,0,0);
  double pressfreeze=EOS->pressure(Efreeze,0,0,0);
  cout<<pressfreeze<<endl;
  for (int i=0;i<par::nzeta;i++){
    double zeta=i*dzeta;
    //    cout<<zeta<<endl;
    for (int j=0;j<par::nphi;j++){
      double phi=j*dphi;
      for (int k=0;k<par::ntheta;k++){
	double theta=k*dtheta;
	//	cout<<Efreeze<<" "<<phi<<" "<<zeta<<" "<<theta<<endl;
	double radius=dist(Efreeze,phi,zeta,theta);
	//	cout<<i<<" "<<j<<" "<<k<<" "<<radius<<endl;
	disthyper->put(radius/par::ren,i,j,k);
	double rho=par::rmin+radius*cos(zeta)*sin(theta);
	double x=rho*cos(phi);double y=rho*sin(phi);
	double t=tmin+radius*sin(zeta)*sin(theta);
	double eta=radius*cos(theta)/lambda;
	double ux=value(t,x,y,eta,uxmap);
	double uy=value(t,x,y,eta,uymap); //
	double ue=value(t,x,y,eta,uemap);
	uxhyper->put(ux,i,j,k);
	uyhyper->put(uy,i,j,k);
	double Y=ue;
	Yhyper->put(Y,i,j,k);	double o;
	o=value(t,x,y,eta,pixxmap)/entalpyfreeze;
	pixxhyper->put(o,i,j,k);
	o=value(t,x,y,eta,pixymap)/entalpyfreeze;
	pixyhyper->put(o,i,j,k);
	o=value(t,x,y,eta,pixemap)/entalpyfreeze;
	pixehyper->put(o,i,j,k);
	o=value(t,x,y,eta,piyymap)/entalpyfreeze;
	piyyhyper->put(o,i,j,k);
	o=value(t,x,y,eta,piyemap)/entalpyfreeze;
	piyehyper->put(o,i,j,k);
	o=value(t,x,y,eta,pieemap)/entalpyfreeze;
	pieehyper->put(o,i,j,k);
	o=value(t,x,y,eta,bulkmap)/pressfreeze;
	bulkhyper->put(o,i,j,k);
	  
	
      }
    }
  }
    
}

void evolve::put(int dumpstep,grid* h,double tau,EoS* EOS){
  //  for(int k=0;k<nneta;k++){int i=h->getNX()/2; int j=h->getNY()/2;
  //  Cell *c = h->getCell(i,j,k);
    //    double eta=h->getZ(k);
  //  double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;
  //   h->getCMFvariables(c,tau, e, nb, nq, ns, vx, vy, vz);
    //    cout<<k<<" "<<eta<<" "<<vz<<endl;
  //  }
  int ix,iy,iz;
  double maxe=0;
  for (int i=0;i<nnx;i++){int ii=2*i;
    for(int j=0;j<nny;j++){int jj=2*j;
      for(int k=0;k<nneta;k++){int kk=2*k;
	cell *c = h->get_cell(ii,jj,kk);
	double o;
	double eta=h->get_eta(kk);	
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;
	h->getCMFvariablesUmunu(c,tau, e, nb, nq, ns, vx, vy, vz);
	if (maxe<e){maxe=e;ix=ii;iy=jj;iz=kk;}
  	energymap[dumpstep]->put(e,i,j,k);
	uxmap[dumpstep]->put(vx,i,j,k);
	uymap[dumpstep]->put(vy,i,j,k);
	//		uemap[dumpstep]->put(vz-eta,i,j,k);
	uemap[dumpstep]->put(vz,i,j,k);
	pixxmap[dumpstep]->put(c->get_pi(1,1),i,j,k);
	pixymap[dumpstep]->put(c->get_pi(1,2),i,j,k);
	o=c->get_pi(1,3)*tau*cosh(eta)+c->get_pi(0,1)*sinh(eta);
	pixemap[dumpstep]->put(o,i,j,k);
	piyymap[dumpstep]->put(c->get_pi(2,2),i,j,k);
	o=c->get_pi(2,3)*tau*cosh(eta)+c->get_pi(0,2)*sinh(eta);
	piyemap[dumpstep]->put(o,i,j,k);
	o=c->get_pi(3,3)*pow(tau,2)*cosh(eta)*cosh(eta)
	  +c->get_pi(0,0)*sinh(eta)*sinh(eta)+
	  c->get_pi(0,3)*tau*sinh(2.*eta);
	pieemap[dumpstep]->put(c->get_pi(3,3)*pow(tau,2),i,j,k);
	bulkmap[dumpstep]->put(c->get_Pi(),i,j,k);
      }
    }
 }
/*  int ic=(h->get_nx()-1)/2,jc=(h->get_ny()-1)/2,kc=(h->get_nz()-1)/2;
  // cout<<i<<" "<<j<<" "<<k<<endl;
  cell *c = h->get_cell(ic,jc,kc);
  double   eta=h->get_z(kc);//eta=0;	
  double e, p, nb, nq, ns,  vx, vy, vz;
  h->getCMFvariablesUmunu(c,tau, e, nb, nq, ns, vx, vy, vz);
  cout<<" central temperature "<<EOS->TE(e)<<" maximal temp. "<<EOS->TE(maxe)<<
    " at "<<h->getX(ix)<<" "<<h->getY(iy)<<" "<<h->getZ(iz)<<endl; 
 int i=ic-5;int  j=ic-5;int k=kc-5;
  c = h->get_cell(i,j,k);
  eta=h->get_z(k);//	eta=0;
  //double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;
  h->getCMFvariablesUmunu(c,tau, e, nb, nq, ns, vx, vy, vz);
  double u0,ux,uy,uz;
  ux=vx;
  uy=vy;
  uz=sqrt(1+ux*ux+uy*uy)*sinh(vz);
  u0=sqrt(1+ux*ux+uy*uy)*cosh(vz);
  //  cout<<i<<" "<<j<<" "<<k<<" "<<h->getX(i)<<" "<<h->getY(j)
  //    <<" "<<eta<<endl;
  double o00,o01,o02,o03,o11,o12,o13,o22,o23,o33;
  o00=c->getpi(0,0)*cosh(eta)*cosh(eta)
    +c->getpi(3,3)*sinh(eta)*sinh(eta)+
    c->getpi(0,3)*sinh(2.*eta);
  o01=c->getpi(0,1)*cosh(eta)+c->getpi(1,3)*sinh(eta);
  o02=c->getpi(0,2)*cosh(eta)+c->getpi(2,3)*sinh(eta);
  o03=c->getpi(0,3)*cosh(2.*eta)
    +c->getpi(0,0)*sinh(eta)*cosh(eta)+
    c->getpi(3,3)*sinh(eta)*cosh(eta);
  o11=c->getpi(1,1);
  o12=c->getpi(1,2);
  o22=c->getpi(2,2);
  o13=c->getpi(1,3)*cosh(eta)+c->getpi(0,1)*sinh(eta);
  o23=c->getpi(2,3)*cosh(eta)+c->getpi(0,2)*sinh(eta);
  o33=c->getpi(3,3)*cosh(eta)*cosh(eta)
    +c->getpi(0,0)*sinh(eta)*sinh(eta)+
	  c->getpi(0,3)*sinh(2.*eta);
  p=EOS->PE(e);
  cout<<ux<<" "<<uy<<" "<<uz<<" "<<e<<endl;
  cout<<" trace "<<o00-o11-o22-o33<<endl;
  cout<<u0*o00-ux*o01-uy*o02-uz*o03<<" ";
  cout<<u0*o01-ux*o11-uy*o12-uz*o13<<" ";
  cout<<u0*o02-ux*o12-uy*o22-uz*o23<<" ";
  cout<<u0*o03-ux*o13-uy*o23-uz*o33<<" pipi^1/2 ";
  cout<<sqrt(o00*o00+o11*o11+o22*o22+o33*o33
	     +2.*o12*o12+2.*o13*o13+2.*o23*o23
	     -2.*o01*o01-2.*o02*o02-2.*o03*o03)<<endl;
  i=ic-5;j=ic-5;k=kc+5;
  c = h->getCell(i,j,k);
  eta=h->getZ(k);//	eta=0;
  //double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;
  h->getCMFvariablesUmunu(c,tau, e, nb, nq, ns, vx, vy, vz);
  // double u0,ux,uy,uz;
  ux=vx;
  uy=vy;
  uz=sqrt(1+ux*ux+uy*uy)*sinh(vz);
  u0=sqrt(1+ux*ux+uy*uy)*cosh(vz);
  //  cout<<i<<" "<<j<<" "<<k<<" "<<h->getX(i)<<" "<<h->getY(j)
  //    <<" "<<eta<<endl;
  //  double o00,o01,o02,o03,o11,o12,o13,o22,o23,o33;
  o00=c->getpi(0,0)*cosh(eta)*cosh(eta)
    +c->getpi(3,3)*sinh(eta)*sinh(eta)+
    c->getpi(0,3)*sinh(2.*eta);
  o01=c->getpi(0,1)*cosh(eta)+c->getpi(1,3)*sinh(eta);
  o02=c->getpi(0,2)*cosh(eta)+c->getpi(2,3)*sinh(eta);
  o03=c->getpi(0,3)*cosh(2.*eta)
    +c->getpi(0,0)*sinh(eta)*cosh(eta)+
    c->getpi(3,3)*sinh(eta)*cosh(eta);
  o11=c->getpi(1,1);
  o12=c->getpi(1,2);
  o22=c->getpi(2,2);
  o13=c->getpi(1,3)*cosh(eta)+c->getpi(0,1)*sinh(eta);
  o23=c->getpi(2,3)*cosh(eta)+c->getpi(0,2)*sinh(eta);
  o33=c->getpi(3,3)*cosh(eta)*cosh(eta)
    +c->getpi(0,0)*sinh(eta)*sinh(eta)+
	  c->getpi(0,3)*sinh(2.*eta);
  p=EOS->PE(e);
  cout<<ux<<" "<<uy<<" "<<uz<<" "<<e<<endl;
  cout<<" trace "<<o00-o11-o22-o33<<endl;
  cout<<u0*o00-ux*o01-uy*o02-uz*o03<<" ";
  cout<<u0*o01-ux*o11-uy*o12-uz*o13<<" ";
  cout<<u0*o02-ux*o12-uy*o22-uz*o23<<" ";
  cout<<u0*o03-ux*o13-uy*o23-uz*o33<<" pipi^1/2 ";
  cout<<sqrt(o00*o00+o11*o11+o22*o22+o33*o33
	     +2.*o12*o12+2.*o13*o13+2.*o23*o23
	     -2.*o01*o01-2.*o02*o02-2.*o03*o03)<<endl;*/
};


void evolve::setAuAu200ev(int cent,int &nmin, int &nmax){
  T0center=0.209;
  sqrts=200.;
  strcpy(system,"AuAut"); 
  strcpy(device,"RHIC");
  muB=28.5/.165;muI=-0.9/.165;muS=6.9/.165;muC=0.;
 switch (cent){
 case 1: cmin=0.; cmax=5.; b=2.1; nmin=325; nmax=1000; break;
  case 2: cmin=5.; cmax=10.; b=3.1; nmin=276; nmax=324; break;
  case 3: cmin=10.; cmax=20.; b=5.0; nmin=199; nmax=275;  break;
  case 4: cmin=20.; cmax=30.; b=7.1; nmin=143; nmax=198;  break;
  case 5: cmin=30.; cmax=40.; b=8.7; nmin=96; nmax=142; break;
  case 6: cmin=40.; cmax=50.; b=9.9; nmin=65; nmax=95; break;
  case 7: cmin=50.; cmax=60.; b=10.9; nmin=40; nmax=64; break;
  case 8: cmin=60.; cmax=70.; b=11.9; nmin=22; nmax=39; break;
  case 10: cmin=35.; cmax=35.; b=100.; nmin=100; nmax=100; break;
  }  
}

double evolve::value(double t,double x,double y,double eta,
	      sVector3D* tab[]){
  if (fabs(x)>maxx||fabs(y)>maxy||fabs(eta)>maxeta
      ||(t-tmin)>maxtimesaved||t<tmin){
    return 0.0;
  }
  else{
    int it=(int) ((t-tmin)/dtsave);
    double ddt=t-it*dtsave-tmin;
    it=minimum(it,timesaved-1);
    return (1.-ddt)*tab[it]->Interpolate(x,y,eta)+
      ddt*tab[it+1]->Interpolate(x,y,eta);
  }
}



double evolve::dist(double enF,double phi, double xi, double theta){
  double distit0=0,distit1=distmax;
  for (int i=1;i<par::distiter;i++){
    if (endens((distit0+distit1)/2.,phi,xi,theta)<enF){
      distit1=(distit0+distit1)/2.;
      //  cout<<i<<" "<<distit0<<" "<<distit1<<" "<<endens(distit1,phi,xi,theta,grid)-enF<<endl;
    }
    else{
      distit0=(distit0+distit1)/2.;
    
      //     cout<<i<<" "<<distit0<<" "<<distit1<<" "<<endens(distit0,phi,xi,theta,grid)-enF<<endl;
}
  }
  return (distit0+distit1)/2.;
}



double evolve::endens(double d,double phi, double xi, double theta){
  double rho=par::rmin+d*cos(xi)*sin(theta);
  double x=rho*cos(phi);double y=rho*sin(phi);
  double t=tmin+d*sin(xi)*sin(theta);
  double eta=d*cos(theta)/lambda;
  return value(t,x,y,eta,energymap);
}



double evolve::etavis(double en,EoS* EOS,grid* grid){
  // return 1e-10;
  //   return minimum(3.* EOS.pe(en) *t0/4./par::ren,par::alpha *EOS.se(en))/
  //     (1.+exp(-(EOS.te(en)-par::Tmintau)/par::Tsc));
  //   return par::alpha *EOS.se(en)/
  //    (1.+exp(-(EOS.te(en)-par::Tmintau)/par::Tsc));
  // if (par::dissipativeHG){ return
  //     (1./(1.+exp(-(EOS.te(en)-par::Tcut)/par::deltaTcut)))*
  //   (par::alpha+(par::alphaHG-par::alpha)/(1+exp((EOS.te(en)-par::THG)
  //						  /par::deltaTHG)))*EOS.se(en);
  //  }
  // else{
  double etaS, zetaS;
  grid->trcoef->getEta(en,EOS->temperature(en,0,0,0),etaS,zetaS);
   return etaS*EOS->entropy(en,0,0,0);
   //  }
  // return par::alpha *EOS.se(en);
}

double evolve::zetavis(double en,EoS* EOS,grid* grid){
  // return 1e-10;
  //   return minimum(3.* EOS.pe(en) *t0/4./par::ren,par::alpha *EOS.se(en))/
  //     (1.+exp(-(EOS.te(en)-par::Tmintau)/par::Tsc));
  //   return par::alpha *EOS.se(en)/
  //    (1.+exp(-(EOS.te(en)-par::Tmintau)/par::Tsc));
  // if (par::dissipativeHG){ return
  //     (1./(1.+exp(-(EOS.te(en)-par::Tcut)/par::deltaTcut)))*
  //   (par::alpha+(par::alphaHG-par::alpha)/(1+exp((EOS.te(en)-par::THG)
  //						  /par::deltaTHG)))*EOS.se(en);
  //  }
  // else{
  double etaS, zetaS;
  grid->trcoef->getEta(en,EOS->temperature(en,0,0,0),etaS,zetaS);
   return zetaS*EOS->entropy(en,0,0,0);
   //  }
  // return par::alpha *EOS.se(en);
}


void evolve::hypersurface(const char* dirname,int eventcount1){
  std::ofstream oh;
  char fullfilenameloc[128]; 
  char filename[58]; 
  strcpy(fullfilenameloc,dirname);
//  strcat(fullfilenameloc,".");
//  strcat(fullfilenameloc,"./"); 
  strcat(fullfilenameloc,device);    
  strcat(fullfilenameloc,system);   
    if (RDSFLAGA==1){
      
      sprintf(filename,"%.0fgrzad%.1fe%.2fTi%.0ft%.2fTf%.0fe%03i.xml",
	      sqrts,b,etafreeze,T0center*1000,tmin,Tfreeze*1000
	      ,eventcount1);   
      
 }
    else{
      sprintf(filename,"%.0fsrzad%.1fe%.2fTi%.0ft%.2fTf%.0fe%03i.xml",
	      sqrts,b,etafreeze,T0center*1000,tmin,Tfreeze*1000
	      ,eventcount1);    
    }

strcat(fullfilenameloc,filename);
cout<<fullfilenameloc<<endl;


  oh.open(fullfilenameloc);
  if((oh) && (oh.is_open())) {
    oh<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"<<endl;
    oh<<"<HYPERSURFACE filename=\""<<fullfilenameloc<<"\" version=\"1.0\
\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
xsi:noNamespaceSchemaLocation=\"http://therminator2.ifj.edu.pl/\
hypersurface.xsd\">"<<endl;
    oh<<"<PARAMETERS>"<<endl;

    oh<<"  <PARAMETER name=\"Tau_i\" unit=\"fm\">"<<tmin<<"</PARAMETER>"<<endl;
    oh<<"  <PARAMETER name=\"Lambda\" unit=\"fm\">"<<lambda<<"</PARAMETER>"
      <<endl;
    oh<<"  <PARAMETER name=\"Temperature\" unit=\"MeV\">"
      <<Tfreeze*1000.<<"</PARAMETER>"<<endl;
    oh<<"  <PARAMETER name=\"Mu_B\" unit=\"MeV\">"<<muB*Tfreeze<<"</PARAMETER>"
      <<endl;
    oh<<"  <PARAMETER name=\"Mu_I\" unit=\"MeV\">"<<muI*Tfreeze<<"</PARAMETER>"
      <<endl;
    oh<<"  <PARAMETER name=\"Mu_S\" unit=\"MeV\">"<<muS*Tfreeze<<"</PARAMETER>"
      <<endl;
    oh<<"  <PARAMETER name=\"Mu_C\" unit=\"MeV\">"<<muC*Tfreeze<<"</PARAMETER>"
      <<endl;
    oh<<"  <PARAMETER name=\"device\">"<<device<<"</PARAMETER>"
      <<endl;
    oh<<"  <PARAMETER name=\"colliding_system\">"<<system<<"</PARAMETER>"
      <<endl;
    oh<<" <PARAMETER name=\"colliding_energy\" unit =\" GeV \" >"
      <<sqrts<<"</PARAMETER>"<<endl;
    oh<<"  <PARAMETER name=\"centrality_min\" unit=\"%\">"<<cmin
      <<"</PARAMETER>"<<endl;
    oh<<"  <PARAMETER name=\"centrality_max\" unit=\"%\">"
      <<cmax<<"</PARAMETER>"<<endl;
    oh<<"  <PARAMETER name=\"impact_parameter\" unit=\"fm\">"
      <<b<<"</PARAMETER>"<<endl;
    oh<<"  <PARAMETER name=\"temperature_at_center\" unit=\"MeV\">"
      <<T0center<<"</PARAMETER>"<<endl;
   oh<<"  <PARAMETER name=\"eps2\" unit=\"\">"<<inieps2_2<<"</PARAMETER>"
      <<endl;
   oh<<"  <PARAMETER name=\"eps3\" unit=\"\">"<<inieps3_3<<"</PARAMETER>"
      <<endl;
   oh<<"  <PARAMETER name=\"psi3\" unit=\"\">"<<psi3<<"</PARAMETER>"
      <<endl;
   oh<<"  <PARAMETER name=\"cs^2\" unit=\"\">"<<cs2freeze<<"</PARAMETER>"
      <<endl;
   oh<<"  <PARAMETER name=\"Nw\" unit=\"\">"<<Nw<<"</PARAMETER>"
     <<endl;
     oh<<"  <PARAMETER name=\"eta/s\" unit=\"\">"<<etaovs<<"</PARAMETER>"
      <<endl;
     oh<<"  <PARAMETER name=\"zeta/s\" unit=\"\">"<<zetaovs<<"</PARAMETER>"<<endl;
       oh<<"  <PARAMETER name=\"eta/s(freeze)\" unit=\"\">"<<etafreeze<<"</PARAMETER>" <<endl;


    oh<<"</PARAMETERS>"<<endl;


    oh<<"<DESCRIPTION>"<<endl;


    oh<<"</DESCRIPTION>"<<endl;

    writeVector3D(oh,"Distance","GeV^-1",disthyper);
    writeVector3D(oh,"FluidUx"," ",uxhyper);
    writeVector3D(oh,"FluidUy"," ",uyhyper);
    writeVector3D(oh,"FluidRap"," ",Yhyper);
      writeVector3D(oh,"Pixx"," ",pixxhyper);
      writeVector3D(oh,"Pixy"," ",pixyhyper);
      writeVector3D(oh,"Pixe"," ",pixehyper);
      writeVector3D(oh,"Piyy"," ",piyyhyper);
      writeVector3D(oh,"Piye"," ",piyehyper);
      writeVector3D(oh,"Piee"," ",pieehyper);
      writeVector3D(oh,"PI"," ",bulkhyper);
    
  }
}


void evolve::writeVector3D(std::ofstream& oh,const char* name,const char* unit,
			   Vector3D* hyp){
  oh<<"<VECTOR3D name=\""<<name<<"\" unit=\""<<unit<<"\">"<<endl;
  oh<<"  <AXIS name=\"Zeta\">"<<endl;
  oh<<"    <DETAIL name=\"min\" unit=\"rad\">"<<hyp->GetXMin()<<"</DETAIL>"
    <<endl;
  oh<<"    <DETAIL name=\"max\" unit=\"rad\">"<<hyp->GetXMax()<<"</DETAIL>"
    <<endl;
  oh<<"    <DETAIL name=\"pts\">"<<hyp->GetXPts()<<"</DETAIL>"<<endl;
  oh<<"  </AXIS>"<<endl;
  oh<<"  <AXIS name=\"Phi\">"<<endl;
  oh<<"    <DETAIL name=\"min\" unit=\"rad\">"<<hyp->GetYMin()<<"</DETAIL>"
    <<endl;
  oh<<"    <DETAIL name=\"max\" unit=\"rad\">"<<hyp->GetYMax()<<"</DETAIL>"
    <<endl;
  oh<<"    <DETAIL name=\"pts\">"<<hyp->GetYPts()<<"</DETAIL>"<<endl;
  oh<<"  </AXIS>"<<endl;
  oh<<"  <AXIS name=\"Theta\">"<<endl;
  oh<<"    <DETAIL name=\"min\" unit=\"rad\">"<<hyp->GetZMin()<<"</DETAIL>"
    <<endl;
  oh<<"    <DETAIL name=\"max\" unit=\"rad\">"<<hyp->GetZMax()<<"</DETAIL>"
    <<endl;
  oh<<"    <DETAIL name=\"pts\">"<<hyp->GetZPts()<<"</DETAIL>"<<endl;
  oh<<"  </AXIS>"<<endl;
  oh<<"  <DATA points=\""<<hyp->GetXPts()*hyp->GetYPts()*hyp->GetZPts()
<<"\">"<<endl;
  for (int ix=0;ix<hyp->GetXPts();ix++){
    for (int iy=0;iy<hyp->GetYPts();iy++){
      for(int ieta=0;ieta<hyp->GetZPts();ieta++){
	oh<<(*hyp)(ix,iy,ieta)<<endl;
      }
    }
  }
  oh<<"  </DATA>"<<endl;
  oh<<"</VECTOR3D>"<<endl;



}

/*
double evolve::readeventquark(std::ifstream& eventfile,int ievent){
  //  std::cout<<"eventfile"<<std::endl;
  std::cout<<" event n. "<<ievent<<std::endl;
  eventcount=ievent;
  float rds,p1,p2,p3,p4;float rdssum=0;
  eventfile>>rds>>Nw>>p1>>p2>>p3>>p4;rds=rds*.5;
  std::cout<<"event with  "<<Nw<<" wounded quarks   RDS="<<rds;
  for(int i=0;i<Nw;i++){
    eventfile>>xnw[i]>>ynw[i]>>cnw[i]>>p2>>gw[i];
    gw[i]=fabs(gw[i]);
    //gw[i]=1.;
    //rdssum+=fabs(gw[i]);
    rdssum+=gw[i];
    //   std::cout<<xnw[i]<<" "<<ynw[i]<<" "<<cnw[i]<<std::endl;
  }
  std::cout<<" rds sum " <<rdssum <<std::endl;
  eventfile>>p2;
  return rdssum;
}

void evolve::setPbPbquark(int cent,double &nmin, double &nmax,EoS *EOS,double t0){  T0center=0.188; etas=0.;
  float efac;
  if (t0<.255){
    efac=0.88;
  }
  else{
    efac=1.0;
  }
  s0center=EOS->SE(EOS->ET(T0center));
  // T0center=EOS.TE(EOS.ES((0.6*efac*s0center/t0)));
  cout<<"central temperature parameter "<<T0center<<endl;
  sqrts=2760;
  //  x0=2.4; 
  x0=2.3; 
  sig0=1.4;
  strcpy(system,"PbPbq");
  strcpy(device,"LHC");
  muB=0;muI=0.;muS=0.;muC=0.;
 switch (cent){
 case 1: cmin=0.; cmax=5.; b=2.1; nmin=983.5; nmax=20000; break;
  case 2: cmin=5.; cmax=10.; b=3.1; nmin=818.9; nmax=983.5; break;
  case 3: cmin=10.; cmax=20.; b=5.0; nmin=568.8; nmax=818.9;  break;
  case 4: cmin=20.; cmax=30.; b=7.1; nmin=384.0; nmax=568.8;  break;
  case 5: cmin=30.; cmax=40.; b=8.7; nmin=248.0; nmax=384.0; break;
  case 6: cmin=40.; cmax=50.; b=9.9; nmin=149.3; nmax=248.0; break;
  case 7: cmin=50.; cmax=60.; b=10.9; nmin=81.7; nmax=149.3; break;
  case 8: cmin=60.; cmax=70.; b=11.9; nmin=40.3; nmax=81.7; break;
  case 9: cmin=70.; cmax=80.; b=12.9; nmin=17.5; nmax=40.3; break;
  case 10: cmin=0.; cmax=1.; b=101.; nmin=1151.2; nmax=20000; break;
  case 11: cmin=5.; cmax=6.; b=106; nmin=950.7; nmax=983.5; break;
  case 12: cmin=10.; cmax=11.; b=111.; nmin=818.9; nmax=818.9; break;
  case 13: cmin=20.; cmax=21.; b=121.; nmin=546.5; nmax=568.8; break;
  case 14: cmin=30.; cmax=31.; b=131.; nmin=368.4; nmax=384.0; break;
  case 15: cmin=40.; cmax=41.; b=121.; nmin=546.5; nmax=248.0; break;
  case 16: cmin=50.; cmax=51.; b=131.; nmin=140.8; nmax=149.3; break;
  case 17: cmin=60.; cmax=61.; b=131.; nmin=76.6; nmax=81.7; break;
  case 18: cmin=70.; cmax=71.; b=131.; nmin=37.3; nmax=40.3; break;
  case 19: cmin=80.; cmax=81.; b=131.; nmin=16.13; nmax=17.5; break;

  }  
}


double evolve::ebegamma(double x,double y,double eta){

 double ng;
 
  const double   sz=0.3;
 ng=s0center*(2.*0.4*0.4)/2./sz/sz;
double   etabeam=log(sqrts/par::nucleonmass);

  double etarad=etabeam-etas;
  double fsidepl;
  double fsidemi;
  //  std::cout<<x0<<" "<<sig0<<" "<<std::endl;
  if (eta<-etarad){
    fsidepl=0.;fsidemi=2.;}
  else 
    if (eta>etarad){
          fsidepl=2.;fsidemi=0.;}
    else{
      fsidepl=1.+eta/etarad;
      fsidemi=1.-eta/etarad;
    }
    double longprof=exp(-((fabs(eta) - x0)*(fabs(eta) - x0)
			  * HeavisideTheta(fabs(eta) - x0)/(2.* sig0*sig0)));
    double sum=0.;
   for (int i=0;i<Nw;i++){
	sum=sum+ng
	  *exp(-((x-xnw[i])*(x-xnw[i])+(y-ynw[i])*(y-ynw[i]))/2./sz/sz)
	  *(cnw[i]>0 ? fsidepl : fsidemi)*longprof*2.*gw[i];
      } 
     return sum;
}


void evolve::setebe(grid *f, EoS *eos, double tau0) {
  double e, nb, nq, vx = 0., vy = 0., vz = 0.;
  Cell *c;

 
  //--------------
  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);	
	e=eos->ES(ebegamma(x,y,eta));
        if (e < 0.5) e = 0.0;
        vx = vy = 0.0;
        nb = nq = 0.0;
        vz = 0.0;

      avv_num += sqrt(vx * vx + vy * vy) * e;
      avv_den += e;

        c->setPrimVar(eos, tau0, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0.);
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
        c->saveQprev();

        if (e > 0.) c->setAllM(1.);
      }
  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "total energy = " << Etotal *f->getDx() * f->getDy() * f->getDz() *
                                   tau0 << endl;
}


void evolve::writefinal(grid *f,EoS *eos,double eta,double t0,char* dirname){
  cout<<"writting hypersurface at eta="<<eta<<endl;
  char filenamex[30];char fullfilenamelocx[120];
  char filenamey[30];char fullfilenamelocy[120];
  char filenameux[30];char fullfilenamelocux[120];
  char filenameuy[30];char fullfilenamelocuy[120];
  sprintf(filenamex,"hypx%2.2f.dat",eta);
  sprintf(filenamey,"hypy%2.2f.dat",eta);
  strcpy(fullfilenamelocx,dirname);
  strcat(fullfilenamelocx,filenamex);
  strcpy(fullfilenamelocy,dirname);
  strcat(fullfilenamelocy,filenamey);
  sprintf(filenameux,"hypux%2.2f.dat",eta);
  sprintf(filenameuy,"hypuy%2.2f.dat",eta);
  strcpy(fullfilenamelocux,dirname);
  strcat(fullfilenamelocux,filenameux);
  strcpy(fullfilenamelocuy,dirname);
  strcat(fullfilenamelocuy,filenameuy);
  std::cout<<fullfilenamelocx<<std::endl;
  std::ofstream filedumpx(fullfilenamelocx);
  std::ofstream filedumpy(fullfilenamelocy);
  std::ofstream filedumpux(fullfilenamelocux);
  std::ofstream filedumpuy(fullfilenamelocuy);
  for (double xd=- f->getmaxX() ;xd<=f->getmaxX()+0.00001*f->getmaxX();xd=
	 xd+(f->getmaxX()*2.)/(f->getNX()-1)){
    for (double t=t0;t<=maxtimesaved+0.01*f->getDt()+t0;t=t+par::dtsave){
 
      double  o=eos->TE(value(t,xd,0.,eta,energymap));
      filedumpx<<xd<<" "<<t<<" "<<" "<<o<<std::endl;
      o=eos->TE(value(t,0.,xd,eta,energymap));
      filedumpy<<xd<<" "<<t<<" "<<" "<<o<<std::endl;
      o=(value(t,xd,0.0,eta,uxmap));
      filedumpux<<xd<<" "<<t<<" "<<" "<<o<<std::endl;
      o=(value(t,0.,xd,eta,uymap));
      filedumpuy<<xd<<" "<<t<<" "<<" "<<o<<std::endl;
    }
  }
}


void evolve::writespacetimehistory(grid *f,EoS *eos,double t0,char* dirname){
  cout<<"writing spacetime history"<<endl;
  char filename[30];char fullfilenameloc[120];
  sprintf(filename,"spacetimehistory.dat");
  strcpy(fullfilenameloc,dirname);
  strcat(fullfilenameloc,filename);
  std::cout<<fullfilenameloc<<std::endl;
  std::ofstream filedump(fullfilenameloc);

  for(double t=t0; t<=maxtimesaved+0.01*f->getDt()+t0; t=t+2.5*par::dtsave){
   for(double xd=-f->getmaxX(); xd<=f->getmaxX()+0.00001*f->getmaxX();
              xd=xd+(f->getmaxX()*2.)/(f->getNX()-1)){
    for(double yd=-f->getmaxY(); yd<=f->getmaxY()+0.00001*f->getmaxY();
               yd=yd+(f->getmaxY()*2.)/(f->getNY()-1)){
     for(double ed= -f->getmaxZ(); ed<=f->getmaxZ()+0.00001*f->getmaxZ();
                ed=ed+(f->getmaxZ()*2.)/(f->getNZ()-1)){
       double  temp=eos->TE(value(t,xd,yd,ed,energymap));
       double  fldvx=value(t,xd,yd,ed,uxmap);
       double  fldvy=value(t,xd,yd,ed,uymap);
       double  fldve=value(t,xd,yd,ed,uemap);
       filedump<<t<<" "<<xd<<" "<<yd<<" "<<ed<<" "<<temp<<" "<<fldvx<<" "<<fldvy<<" "<<fldve<<std::endl;
       }
      }
     }
  }
}

*/



