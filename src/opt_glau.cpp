///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Optical Glauber Model.                                           //
//  calculation in cartesian coordinate.                             //
//  same species i.e A=B.                                            //
//  No orientation i.e theta = phi =0 ;                              //
//  R.H.S of Normalisation const. is 1 instead of A/B. That will-    //
//  be managed during calculation of npart and ncoll.                //
//                                                                   //
///////////////////////////////////////////////////////////////////////


#include "opt_glau.h"
#include<TF1.h>
#include<cmath>
#include<TF2.h>
#include<TMath.h>
#include "Math/Functor.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "global.h"
#include <fstream>

using namespace std;

opt_glau::opt_glau(idb *_IDB)
{
 IDB = _IDB; 
 set_opt_glau_params();   // setting of initial data like species,energy etc...
 norm_const = Get_Norm_Constant();
}


opt_glau::~opt_glau(){}


double opt_glau::Modf_WoodSaxon(const double* x){
  
  double theta=0;
  double phi=0;
  
  double A1 = (beta2/4.0)*(TMath::Sqrt(5.0/pi));
  double A2 = ((3.0*beta4)/16.0)*(TMath::Sqrt(1.0/pi));
  double ZP = (- ( (TMath::Sin(theta))*x[0] ) -  ((TMath::Cos(theta))*(TMath::Sin(phi))*x[1])  +  
	       ((TMath::Cos(theta))*(TMath::Cos(phi))*x[2]));            
  double rsq = ((x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2])) ;
  double RP = R* (  1+
		    (A1*(((3*ZP*ZP)/rsq)-1))+
		    A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3    ) 
		    );
  
  
  double r = TMath::Sqrt(rsq);
  return 1/(TMath::Exp((r - RP)/a)+1);
}



double opt_glau::Get_Norm_Constant()
{  
  ROOT::Math::Functor wf(this, &opt_glau::Modf_WoodSaxon,3);
  double min[3] = {-3*R,-3*R,-3*R};
  double max[3] = {3.1*R,3.1*R,3.2*R};
  ROOT::Math::IntegratorMultiDim pd;
  pd.SetFunction(wf);
  double val = pd.Integral(min,max);
  return 1/val;
}




double opt_glau::Norm_Function(double* x,double* p)
{
  double A1 = (beta2/4.0)*(TMath::Sqrt(5.0/pi));
  double A2 = ((3.0*beta4)/16.0)*(TMath::Sqrt(1.0/pi));
  double ZP = (- ( (TMath::Sin(p[2]))*(TMath::Cos(p[3]))*p[0] ) -  ((TMath::Sin(p[3]))*p[1])  +  
	       ((TMath::Cos(p[2]))*(TMath::Cos(p[3]))*x[0]));            
  double rsq = ((x[0]*x[0])+p[0]*p[0]+p[1]*p[1]);
  double RP = R* ( 1+
		   (A1*(((3*ZP*ZP)/rsq)-1))+
		   A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3 ) 
		   
		   );
  double r = TMath::Sqrt(rsq);
  return norm_const/(TMath::Exp((r - RP)/a)+1);
  
}


double opt_glau::npartxy(double x,double y)
{
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x+(b/2.0)); 
  f1->SetParameter(1,y); 
  f1->SetParameter(2,mThetaA); 
  f1->SetParameter(3,mPhiA); 
  double TA = f1->Integral(-3*R,3*R,1.0E-5);
  
  f1->SetParameter(0,x-(b/2.0)); 
  f1->SetParameter(2,mThetaB); 
  f1->SetParameter(3,mPhiB); 
  double TB = f1->Integral(-3*R,3*R,1.0E-5);
  double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
  delete f1;
  return  va;
}


double opt_glau::ncollxy(double x,double y)
{
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x+(b/2.0)); 
  f1->SetParameter(1,y); 
  f1->SetParameter(2,mThetaA); 
  f1->SetParameter(3,mPhiA); 
  double TA = f1->Integral(-3*R,3*R,1.0E-05);

  f1->SetParameter(0,x-(b/2.0)); 
  f1->SetParameter(2,mThetaB); 
  f1->SetParameter(3,mPhiB); 
  double TB = f1->Integral(-3*R,3*R,1.0E-05);
  delete f1;
  return  A*B*sigma*TA*TB;
}



void opt_glau::set_ic(grid* f, EoS* eos)
{

  cout<<"      Optical Glauber Model    "<<endl;
  cout<<"      *********************    "<<endl; 
  cout<<IDB->species<<"+"<<IDB->species<<" at "<<IDB->SNN<<"GeV("<<sigma<<"fm)"<<endl;
  cout<<"deformation parameter (Beta-2) : "<<beta2<<" and (Beta-4) : "<<beta4<<endl;

  TRandom* t1=new TRandom();
  t1->SetSeed(0);
  long kss=t1->GetSeed();
  gRandom->SetSeed(kss);
  TF1* f1= new TF1("f1","x",0.0,25.0);

  b = f1->GetRandom(bmin,bmax);
  cout<<"Impact parameter : "<<b<<" fm"<<endl;
  double Npart = npart();
  double Ncoll = ncoll();
  cout<<"No. of participants : "<<Npart<<endl;
  cout<<"No. of binary collisions : "<<Ncoll<<endl;
  cout<<"Multiplicity scaling factor eps0 : "<<eps0<<endl;

  std::ofstream File0;
  File0.open("hydro_output/optical_glauber_ic_dist.dat");

  // output will be a input to music
  File0<<"#"<<"\t"<<"optical_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<IDB->neta<<"\t"<<"nx="<<"\t"<<IDB->nx<<"\t"<<"ny="<<"\t"<<IDB->ny
       <<"\t"<<"deta="<<"\t"<<IDB->deta<<"\t"<<"dx="<<"\t"<<IDB->dx<<"\t"<<"dy="<<"\t"<<IDB->dy<<endl;

  double total_deposited = 0.0;  // total deposited energy
  
  cell* c;
  for(int i=0; i<IDB->nx; i++)
    for(int j=0; j<IDB->ny; j++)
      for(int k=0; k<IDB->neta; k++)
	{
          c = f->get_cell(i,j,k);         
	  
	  double x_ = IDB->xmin + i*IDB->dx;
	  double y_ = IDB->ymin + j*IDB->dy;
          double eta = IDB->etamin + k*IDB->deta;
          double eps;

          double nchxy =(1-X_hard) * npartxy(x_, y_) + X_hard * ncollxy(x_, y_); 
                                                              // two component energy deposition

          double H_eta = exp(  - pow( fabs(eta) - IDB->eta_platue / 2.0, 2 )  /  
                                  ( 2 * pow(IDB->eta_fall,2) ) *  theta(fabs(eta)-IDB->eta_platue/2) );


          eps = nchxy*H_eta;
                                                                   // rapidity distribution 
                                                                   //https://arxiv.org/pdf/0902.4121.pdf  (eqn_2.12)
          
          total_deposited = total_deposited + eps*eps0 ;
          //if(abs(x_)<f->get_dx() && abs(y_)<f->get_dy()){cout<<"energy at (0,0,0) : "<<eps0*eps<<endl;}

	  
          double nb= 0; double nq = 0; double ns =0; double vx=0; double vy=0; double vz= 0;
          double utau = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz);
          double ux = utau*vx;
          double uy = utau*vy;
          double uz = utau*vz; 

          if(abs(eta)<0.0001){  // output will be a input to music
          //cout<<eta<<endl;
          File0<<eta<<"\t"<<x_<<"\t"<<y_<<"\t"<<eps0*eps<<"\t"<<utau<<"\t"<<ux<<"\t"<<uy<<"\t"<<uz<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<endl; // writing the dist in a file
          }


          c->set_prim_var(eos,IDB->tau0,eps0*eps, nb, nq,  ns,  vx,  vy,  vz);
	}
  
  
  cout<<"total amount of deposited energy is = "<<total_deposited<<endl;
  cout<<"\n";
  
}

double opt_glau::theta(double _a)
{
  if(_a > 0){return 1;}else{return 0;}
}




double opt_glau::npart_tab(double* x,double* p)
{
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x[0]+(b/2.0)); 
  f1->SetParameter(1,x[1]); 
  f1->SetParameter(2,mThetaA); 
  f1->SetParameter(3,mPhiA); 
  double TA = f1->Integral(-3*R,3*R,1.0E-5);
  
  f1->SetParameter(0,x[0]-(b/2.0)); 
  f1->SetParameter(2,mThetaB); 
  f1->SetParameter(3,mPhiB); 
  double TB = f1->Integral(-3*R,3*R,1.0E-5);
  double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
  delete f1;
  return  va;
}

double opt_glau::npart()
{
  TF2 *f3;
  f3=new TF2("T_AB",this,&opt_glau::npart_tab,-4*R,4*R,-4*R,4*R,0);
  double TAAB=f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-02);
  double N_Part=TAAB;
  delete f3;
  return N_Part;
}


double opt_glau::ncoll_tab(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-09);
        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-09);
        delete f1;
        return  TA*TB;

}



double opt_glau::ncoll(){

        TF2* f2;
        f2= new TF2("TAB",this,&opt_glau::ncoll_tab,-4*R,4*R,-4*R,4*R,0);
        double T_AB = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
        double N_coll= A*B*sigma*(T_AB);
        delete f2;
        return N_coll;
}






