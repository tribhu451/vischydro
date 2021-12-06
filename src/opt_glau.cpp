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

opt_glau::opt_glau(idb *_IDB){
  IDB = _IDB; 
  set_opt_glau_params();   // setting of initial data like species,energy etc...
  norm_const = Get_Norm_Constant();
}


opt_glau::~opt_glau(){}

double opt_glau::theta(double _a)
{
  if(_a < 0){
   return 0;
  }
  else{
   return 1;
  }
}

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



double opt_glau::Get_Norm_Constant(){ 
  ROOT::Math::Functor wf(this, &opt_glau::Modf_WoodSaxon,3);
  double min[3] = {-3*R,-3*R,-3*R};
  double max[3] = {3.1*R,3.1*R,3.2*R};
  ROOT::Math::IntegratorMultiDim pd;
  pd.SetFunction(wf);
  double val = pd.Integral(min,max);
  return 1/val;
}




double opt_glau::Norm_Function(double* x,double* p){
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




double opt_glau::npartxy(double x,double y){
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x+(b/2.0)); 
  f1->SetParameter(1,y); 
  f1->SetParameter(2,mThetaA); 
  f1->SetParameter(3,mPhiA); 
  double TA = f1->Integral(-3*R,3*R,1.0E-9);
  
  f1->SetParameter(0,x-(b/2.0)); 
  f1->SetParameter(2,mThetaB); 
  f1->SetParameter(3,mPhiB); 
  double TB = f1->Integral(-3*R,3*R,1.0E-9);
  double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
  delete f1;
  return  va;
}


double opt_glau::npartxy_min(double x,double y)
{
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x+(b/2.0));
  f1->SetParameter(1,y);
  f1->SetParameter(2,mThetaA);
  f1->SetParameter(3,mPhiA);
  double TA = f1->Integral(-3*R,3*R,1.0E-9);

  f1->SetParameter(0,x-(b/2.0));
  f1->SetParameter(2,mThetaB);
  f1->SetParameter(3,mPhiB);
  double TB = f1->Integral(-3*R,3*R,1.0E-9);
  double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) - (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
  delete f1;
  return  va;
}



double opt_glau::ncollxy(double x,double y){
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x+(b/2.0)); 
  f1->SetParameter(1,y); 
  f1->SetParameter(2,mThetaA); 
  f1->SetParameter(3,mPhiA); 
  double TA = f1->Integral(-3*R,3*R,1.0E-09);
  
  f1->SetParameter(0,x-(b/2.0)); 
  f1->SetParameter(2,mThetaB); 
  f1->SetParameter(3,mPhiB); 
  double TB = f1->Integral(-3*R,3*R,1.0E-09);
  delete f1;
  return  A*B*sigma*TA*TB;
}

double opt_glau::npart_tab(double* x,double* p)
{
  TF1* f1;
  f1=new TF1 ("TA",this,&opt_glau::Norm_Function,-3*R,3*R,4);
  f1->SetParameter(0,x[0]+(b/2.0)); 
  f1->SetParameter(1,x[1]); 
  f1->SetParameter(2,mThetaA); 
  f1->SetParameter(3,mPhiA); 
  double TA = f1->Integral(-3*R,3*R,1.0E-9);
  
  f1->SetParameter(0,x[0]-(b/2.0)); 
  f1->SetParameter(2,mThetaB); 
  f1->SetParameter(3,mPhiB); 
  double TB = f1->Integral(-3*R,3*R,1.0E-9);
  double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
  delete f1;
  return  va;
}


double opt_glau::npart()
{
  TF2 *f3;
  f3=new TF2("T_AB",this,&opt_glau::npart_tab,-3*R,3*R,-3*R,3*R,0);
  double TAAB=f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-05);
  delete f3;
  return TAAB;
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
        return  A*B*sigma*TA*TB;

}



double opt_glau::ncoll(){

        TF2* f2;
        f2= new TF2("TAB",this,&opt_glau::ncoll_tab,-3*R,3*R,-3*R,3*R,0);
        double T_AB = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-05);
        delete f2;
        return T_AB;
}



void opt_glau::set_ic(grid* f, EoS* eos){
  
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
  cout<<"Impact parameter :\t"<<b<<" (fm)"<<endl;
  double Npart = npart();
  double Ncoll = ncoll();
  cout << "No. of participants : " << Npart << endl;
  cout << "No. of binary collisions : " << Ncoll << endl;
  cout << "Entropy scaling factor s0 : " << s0 << endl;
  cout << "Hardness factor in two component glauber : " << X_hard << endl;
  cout << "n_pp for two component glauber : " << n_pp << endl;

  if (IDB->neta != 1){
    cout << "Eta plateau : " << IDB->eta_platue << endl ; 
    cout << "Eta fall : " << IDB->eta_fall << endl ; 
  }

  // boost_invariant_ic(f,eos);
  // rapidity_shifted_ic(f,eos);
  // rapidity_tilted_ic(f,eos,2.2,1);
  // energy_momentum_conserving_ic_by_chun_shen( f, eos, 0.2, 0 ) ;

  
}



void opt_glau::boost_invariant_ic(grid* f, EoS* eos){
  
  std::ofstream File0;
  File0.open("hydro_output/optical_glauber_ic_dist.dat");
  
  // output will be a input to music
  File0<<"#"<<"\t"<<"optical_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<IDB->neta<<"\t"<<"nx="<<"\t"<<IDB->nx<<"\t"<<"ny="<<"\t"<<IDB->ny
       <<"\t"<<"deta="<<"\t"<<IDB->deta<<"\t"<<"dx="<<"\t"<<IDB->dx<<"\t"<<"dy="<<"\t"<<IDB->dy<<endl;
  
  double total_deposited_entropy = 0.0;  // total deposited energy
  double total_deposited_energy  = 0.0;  // total deposited energy
  
  cell* c;
  for(int i=0; i<IDB->nx; i++)
    {
      for(int j=0; j<IDB->ny; j++)
	{
	  double x_ = IDB->xmin + i*IDB->dx;
	  double y_ = IDB->ymin + j*IDB->dy;
	  double nchxy = n_pp * ( (1-X_hard) *  ( npartxy(x_, y_) / 2.0 ) + X_hard * ncollxy(x_, y_) ) ; 
	  
	  for(int k=0; k<IDB->neta; k++)
	    {
              c = f->get_cell(i,j,k);         
	    	   
              double eta;
              if (IDB->neta == 1) { eta = 0.0 ;                      } 
              else                { eta = IDB->etamin + k*IDB->deta; }
	      
	      double H_eta = exp(  - pow( fabs(eta) - IDB->eta_platue / 2.0, 2 )  /  
				   ( 2 * pow(IDB->eta_fall,2) ) *  theta(fabs(eta)-IDB->eta_platue/2) );
	      // rapidity distribution 
	      //https://arxiv.org/pdf/0902.4121.pdf  (eqn_2.12)
	      
	      
	      double nb= 0; double nq = 0; double ns =0; 
	      
	      double eps = eos->entr_2_eps(s0*nchxy*H_eta,nb,nq,ns);  // entropy converted to energy density
	      
	      double vx=0; double vy=0; double vz= 0;
	      double utau = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz);
	      double ux = utau*vx;
	      double uy = utau*vy;
	      double uz = utau*vz; 
	      
	      if(abs(eta)<0.0001){  // output will be a input to music
		//cout<<eta<<endl;
		File0 << eta << "\t" << x_ << "\t" << y_ << "\t" << eos->entropy(eps,nb,nq,ns) 
                      << "\t" << utau << "\t" << ux << "\t" << uy << "\t" << uz
                      << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl; // writing the dist in a file
	      }
	      
	      
              
	      total_deposited_entropy += s0*nchxy*H_eta ;
	      total_deposited_energy  += eps  ;
	      //if(abs(x_)<f->get_dx() && abs(y_)<f->get_dy()){cout<<"energy at (0,0,0) : "<<eps0*eps<<endl;}
	      
	      
	      c->set_prim_var(eos,IDB->tau0,eps, nb, nq,  ns,  vx,  vy,  vz);
	    }
	}
    }
  
  cout<<"total deposited entropy : "<< total_deposited_entropy << endl;
  cout<<"total deposited energy  : "<< total_deposited_energy << endl;
  cout<<"\n";
    
}







void opt_glau::rapidity_shifted_ic(grid* f, EoS* eos){

  std::cout <<"[Info] A rapidity shifted IC condition..." << std::endl ; 
  std::ofstream File0;
  File0.open("hydro_output/optical_glauber_ic_dist.dat");
  
  // output will be a input to music
  File0<<"#"<<"\t"<<"optical_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<IDB->neta<<"\t"<<"nx="<<"\t"<<IDB->nx<<"\t"<<"ny="<<"\t"<<IDB->ny
       <<"\t"<<"deta="<<"\t"<<IDB->deta<<"\t"<<"dx="<<"\t"<<IDB->dx<<"\t"<<"dy="<<"\t"<<IDB->dy<<endl;
  
  double total_deposited_entropy = 0.0;  // total deposited energy
  double eta, x_, y_, vN ;
  double eta_sh = 0 ; 
  double nchxy = 0 ; 
  cell* c;

  for(int k=0; k<IDB->neta; k++){
    eta = IDB->etamin + k*IDB->deta;
    for(int i=0; i<IDB->nx; i++){
      for(int j=0; j<IDB->ny; j++){

	x_ = IDB->xmin + i*IDB->dx;
	y_ = IDB->ymin + j*IDB->dy;
	
	
        if( k==0 ){ // calculate once at mid rapidity 
	   vN = - sqrt( 1 - (4*0.938*0.938) / (200*200) );  // velocity (vz) of each neucleon
	   eta_sh = 0.5 * log ( ( npartxy(x_, y_) + vN*npartxy_min(x_, y_) ) / ( npartxy(x_, y_) - vN*npartxy_min(x_, y_) ) ) ;
	   nchxy = n_pp * ( (1-X_hard) *  ( npartxy(x_, y_) / 2.0 ) + X_hard * ncollxy(x_, y_) ) ;
        }

	double H_eta = exp(  - pow( fabs(eta-eta_sh) - IDB->eta_platue / 2.0, 2 )  /
			     ( 2 * pow(IDB->eta_fall,2) ) *  theta(fabs(eta-eta_sh)-IDB->eta_platue / 2 ) );

	
	
	c = f->get_cell(i,j,k);         
		
	double nb= 0; double nq = 0; double ns =0; 

        // entropy converted to energy density
	double eps = eos->entr_2_eps(s0*nchxy*H_eta,nb,nq,ns);  

	
	double vx=0; double vy=0; double vz= 0;
	double utau = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz);
	double ux = utau*vx;
	double uy = utau*vy;
	double uz = utau*vz; 
	
	 // output will be a input to music
	 File0 << eta << "\t" << x_ << "\t" << y_ << "\t" << s0*nchxy*H_eta 
	       << "\t" << utau << "\t" << ux << "\t" << uy << "\t" << uz
	       << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl; 
	
	total_deposited_entropy += s0*nchxy*H_eta ;
	
	c->set_prim_var(eos,IDB->tau0,eps, nb, nq,  ns,  vx,  vy,  vz);
      }
    }
  }
  
  cout<<"total deposited entropy : " << total_deposited_entropy << endl;

  
}


void opt_glau::rapidity_tilted_ic(grid* f, EoS* eos, double etam, int baryon_density_flag){

  std::cout <<"[Info] A rapidity tilted IC condition..." << std::endl ; 

  double baryon_depo_peak, baryon_depo_right_fall, baryon_depo_left_fall ;   
  baryon_depo_peak        =  3.5 ; 
  baryon_depo_right_fall  =  0.1 ; 
  baryon_depo_left_fall   =  2.0 ;  

  // output will be a input to music
  std::ofstream File0;
  File0.open("hydro_output/optical_glauber_ic_dist.dat");
  File0<<"#"<<"\t"<<"optical_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<IDB->neta<<"\t"<<"nx="<<"\t"<<IDB->nx<<"\t"<<"ny="<<"\t"<<IDB->ny
       <<"\t"<<"deta="<<"\t"<<IDB->deta<<"\t"<<"dx="<<"\t"<<IDB->dx<<"\t"<<"dy="<<"\t"<<IDB->dy<<endl;

  const int ix = ( IDB->nx * 1.0 );
  const int iy = ( IDB->ny * 1.0 );
  double total_deposited_entropy = 0.0;  // total deposited energy
  double total_nb = 0.0 ; 
  double eta, x_, y_, fA, fB ;
  double NA[ix][iy] ; double NB[ix][iy]; double nn_coll[ix][iy];   
  double nchxy = 0 ; 
  cell* c;

  if ( baryon_density_flag != 0 ){
     cout << "[Info] baryon deposition peak position : " << baryon_depo_peak       << endl ; 
     cout << "[Info] baryon deposition right fall    : " << baryon_depo_right_fall << endl ; 
     cout << "[Info] baryon deposition left fall     : " << baryon_depo_left_fall  << endl ; 
  }
 else
   {
      cout <<  "[Info] Baryon flag = 0 " << endl ;  
   }


  for(int k=0; k<IDB->neta; k++){
    eta = IDB->etamin + k*IDB->deta;
    cout << "iz : " << k << "  eta : " <<  eta  << endl ; 
    for(int i=0; i<IDB->nx; i++){
      for(int j=0; j<IDB->ny; j++){
	x_ = IDB->xmin + i*IDB->dx;
	y_ = IDB->ymin + j*IDB->dy;
        if( k==0 ){ // calculate once 
             NA[i][j] =  0.5 * ( npartxy(x_, y_) + npartxy_min(x_, y_) ) ;
             NB[i][j] =  0.5 * ( npartxy(x_, y_) - npartxy_min(x_, y_) ) ;
             nn_coll[i][j] =  ncollxy(x_, y_) ;
        }
             fA = f_F(-eta,etam);
             fB = f_F(eta,etam);
             nchxy = n_pp * ( ( NB[i][j]*fB + NA[i][j]*fA ) * (1-X_hard)  +  X_hard * nn_coll[i][j] );
             double H_eta = exp(  - pow( fabs(eta) - IDB->eta_platue / 2.0, 2 )  /  
               ( 2 * pow(IDB->eta_fall,2) ) *  theta(fabs(eta)-IDB->eta_platue/2) );

	c = f->get_cell(i,j,k);         
		
	double nb = 0 ;
        double nq = 0 ;
        double ns = 0 ; 
        if(baryon_density_flag != 0 ){
            set_eta_0_nb(baryon_depo_peak);
            set_sigma_eta_nb_plus(baryon_depo_right_fall);
            set_sigma_eta_nb_minus(baryon_depo_left_fall);
            nb += baryon_density_eta_envelop_profile_0(eta) * NB[i][j] ; 
            set_eta_0_nb(-baryon_depo_peak);
            set_sigma_eta_nb_plus(baryon_depo_left_fall);
            set_sigma_eta_nb_minus(baryon_depo_right_fall);
            nb += baryon_density_eta_envelop_profile_0(eta) * NA[i][j] ; 
        }


        // entropy converted to energy density
	double eps = eos->entr_2_eps(s0*nchxy*H_eta,nb,nq,ns);  
	double vx = 0 ; 
        double vy = 0 ; 
        double veta = 0;
        // double veta  =  ( eps * sinh(y_L[i][j]) ) / ( eps * cosh(y_L[i][j]) + eos->pressure( eps, nb, nq, ns ) ) ;
	double utau = 1.0 / sqrt( 1.0 - vx*vx - vy*vy - veta*veta );
	double ux   = utau*vx;
	double uy   = utau*vy;
	double ueta = utau * veta / IDB->tau0;
        
	 // output will be a input to music
	 File0 << eta << "\t" << x_ << "\t" << y_ << "\t" << s0*nchxy*H_eta 
	       << "\t" << utau << "\t" << ux << "\t" << uy << "\t" << ueta
	       << "\t" << nb << "\t" << "0" << "\t" << "0" << endl; 
	
        total_nb += nb ; 
	total_deposited_entropy += s0*nchxy*H_eta ;
	
	c->set_prim_var(eos,IDB->tau0,eps, nb, nq,  ns,  vx,  vy,  veta );
      }
    }
  }
  
  cout<<"[Info] total deposited entropy : " << total_deposited_entropy * IDB->dx * IDB->dy * IDB->deta << endl;
  cout<<"[Info] total deposited nb : " << total_nb * IDB->dx * IDB->dy * IDB->deta << endl;
 
}





// chun shen proposed IC by conserving energy and momentum
// Ref : PHYSICAL REVIEW C 102, 014909 (2020)
// ref : arXiv : 2106.08125 
void opt_glau::energy_momentum_conserving_ic_by_chun_shen(grid* f, EoS* eos, double frac, int baryon_density_flag ){

  std::cout << "[Info] An IC with Energy Momentum conservation condition..." << std::endl ;

  double baryon_depo_peak, baryon_depo_right_fall, baryon_depo_left_fall ;   
  baryon_depo_peak        =  3.5 ; 
  baryon_depo_right_fall  =  0.1 ; 
  baryon_depo_left_fall   =  2.0 ;  

  // output will be a input to music 
  std::ofstream out_file;
  out_file.open("hydro_output/optical_glauber_ic_dist.dat");
  out_file<<"#"<<"\t"<<"optical_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<
               "\t"<<IDB->neta<<"\t"<<"nx="<<"\t"<<IDB->nx<<"\t"<<"ny="<<"\t"<<IDB->ny
                     <<"\t"<<"deta="<<"\t"<<IDB->deta<<"\t"<<"dx="<<"\t"<<IDB->dx<<"\t"<<"dy="<<
                         "\t"<<IDB->dy<<endl;

  // Print informations
  cout << "[Info] f parameter : " << frac << endl ; 
  if ( baryon_density_flag != 0 ){
     cout << "[Info] baryon deposition peak position : " << baryon_depo_peak       << endl ; 
     cout << "[Info] baryon deposition right fall    : " << baryon_depo_right_fall << endl ; 
     cout << "[Info] baryon deposition left fall     : " << baryon_depo_left_fall  << endl ; 
  }
 else
   {
      cout <<  "[Info] Baryon flag = 0 " << endl ;  
   }

  const int ix  = ( IDB->nx * 1.0 );
  const int iy  = ( IDB->ny * 1.0 );
  double mN     = 0.938 ;                              // mass of nucleon in GeV
  double y_beam = acosh( IDB->SNN / ( 2 * mN ) ) ; // beam rapidity. s_{NN} in GeV .  
  double eta, x_, y_, C_eta ;
  double TA[ix][iy] ; 
  double TB[ix][iy];           // Thickness functions 
  double M[ix][iy]  ;
  double y_CM[ix][iy] ;
  double y_L[ix][iy] ;
  cell* c;
  double total_deposited_entropy = 0.0 ;            // total deposited energy
  double total_nb                = 0.0 ;

  for(int k=0; k<IDB->neta; k++){
    eta = IDB->etamin + k*IDB->deta;
    std::cout << "iz : " << k 
          << "  eta : " <<  eta  << std::endl ; 
    for(int i=0; i<IDB->nx; i++){
      for(int j=0; j<IDB->ny; j++){
	x_ = IDB->xmin + i*IDB->dx ;
	y_ = IDB->ymin + j*IDB->dy ;
        if( k==0 ){ // calculate once 
             TA[i][j]    =  0.5 * ( npartxy(x_, y_) + npartxy_min(x_, y_) ) ; // TA
             TB[i][j]    =  0.5 * ( npartxy(x_, y_) - npartxy_min(x_, y_) ) ; // TB

             M[i][j]     =  mN * sqrt( TA[i][j] * TA[i][j] +  TB[i][j] * TB[i][j] + 2 * TA[i][j] * TB[i][j] * cosh( 2 * y_beam ) ) ; // M(x,y) -> invariant mass
             y_CM[i][j]  =  atanh( ( TB[i][j] - TA[i][j] ) / ( TB[i][j] + TA[i][j] ) * tanh( y_beam ) ); // center of mass rapidity 
             y_L[i][j]   =  frac * y_CM[i][j] ; 

             C_eta       =  exp( IDB->eta_platue ) * TMath::Erfc( -IDB->eta_fall / sqrt(2.) ) 
                               + exp( -IDB->eta_platue ) * TMath::Erfc( IDB->eta_fall / sqrt(2.) ) ; 
         
             M[i][j]    /=  IDB->tau0 * 
                              ( 2 * sinh(IDB->eta_platue) + sqrt( TMath::Pi() / 2) * IDB->eta_fall * exp( IDB->eta_fall * IDB->eta_fall / 2 ) * C_eta ) ; 
        }


 	double H_eta     =  exp(  - pow( fabs(eta - ( y_CM[i][j] - y_L[i][j] ) ) - IDB->eta_platue  , 2 )  /
			     ( 2 * pow(IDB->eta_fall,2) ) *  theta(fabs(eta - ( y_CM[i][j] - y_L[i][j] ) )-IDB->eta_platue ) );

	
	c  =  f->get_cell(i,j,k); 
        		
	double nb    =  0 ;
        double nq    =  0 ;
        double ns    =  0 ; 

        if(baryon_density_flag != 0 ){
            set_eta_0_nb(baryon_depo_peak);
            set_sigma_eta_nb_plus(baryon_depo_right_fall);
            set_sigma_eta_nb_minus(baryon_depo_left_fall);
            nb += baryon_density_eta_envelop_profile_0(eta) * TB[i][j] ; 
            set_eta_0_nb(-baryon_depo_peak);
            set_sigma_eta_nb_plus(baryon_depo_left_fall);
            set_sigma_eta_nb_minus(baryon_depo_right_fall);
            nb += baryon_density_eta_envelop_profile_0(eta) * TA[i][j] ; 
        }
   
	double eps   =  M[i][j] * H_eta ;  
	double vx    =  0 ; 
        double vy    =  0 ;
        double veta  =  ( eps * sinh(y_L[i][j]) ) / ( eps * cosh(y_L[i][j]) + eos->pressure( eps, nb, nq, ns ) ) ;
	double utau  =  1.0 / sqrt( 1.0 - vx * vx - vy * vy - veta * veta ) ;
	double ux    =  utau * vx ;
	double uy    =  utau * vy ;
	double ueta  =  utau * veta / IDB->tau0 ;

	 // output will be a input to music
	 out_file << eta << "\t" << x_ << "\t" << y_ << "\t" << eps 
	            << "\t" << utau << "\t" << ux << "\t" << uy << "\t" << ueta
	                 << "\t" << nb << "\t" << "0" << "\t" << "0" << endl; 
	
        total_nb += nb ; 
	total_deposited_entropy += M[i][j] * H_eta ;
	
	c->set_prim_var( eos, IDB->tau0, eps, nb, nq,  ns,  vx,  vy,  veta );
      }
    }
  }
  
  cout<<"[Info] total deposited energy : " << total_deposited_entropy * IDB->dx * IDB->dy * IDB->deta << endl;
  cout<<"[Info] total deposited nb : " << total_nb * IDB->dx * IDB->dy * IDB->deta << endl;
  out_file.close();
  
}





double opt_glau::baryon_density_eta_envelop_profile_0(double eta){ 
                                                      // as taken in arxiv:1804.10557

  if( fabs(eta - IDB->etamin) < 0.0001 ) NORM_BARYON_ENVELOP = integrate_baryon_density_eta_envelop_profile_0_over_eta() ;  
  double eta_0_nb = get_eta_0_nb() ;
  double THETA_ARG = eta - eta_0_nb ; 
  if ( fabs(THETA_ARG) < 1E-6 ) THETA_ARG = 1E-6 ; 
  double sigma_eta_nb_plus = get_sigma_eta_nb_plus() ; 
  double sigma_eta_nb_minus = get_sigma_eta_nb_minus() ; 
  //cout << " > > > > > > > > >  Integral value : " << integrate_baryon_density_eta_envelop_profile_0_over_eta() << endl ; 
  return ( 1.0 / NORM_BARYON_ENVELOP )
                    * ( theta( THETA_ARG  ) 
                    * exp( - pow( THETA_ARG , 2 ) / ( 2 * pow( sigma_eta_nb_plus, 2) ) )
                    + theta( -THETA_ARG )
                    * exp( - pow( THETA_ARG, 2 ) / ( 2 * pow( sigma_eta_nb_minus, 2) ) ) ) ; 
  
}

double opt_glau::baryon_density_eta_envelop_profile_0_function(double* x, double* p){ 
                                                      // as taken in arxiv:1804.10557
  double eta_0_nb = get_eta_0_nb() ;
  double sigma_eta_nb_plus = get_sigma_eta_nb_plus() ; 
  double sigma_eta_nb_minus = get_sigma_eta_nb_minus() ; 
  return p[0]*( theta( x[0] - eta_0_nb ) * exp( - pow( eta_0_nb - x[0], 2 ) / ( 2 * pow( sigma_eta_nb_plus, 2) ) )
                   + theta( eta_0_nb - x[0] ) * exp( - pow( eta_0_nb - x[0], 2 ) / ( 2 * pow( sigma_eta_nb_minus, 2) ) ) ) ; 
}


double opt_glau::integrate_baryon_density_eta_envelop_profile_0_over_eta(){
        TF1* f1;
        f1=new TF1 ("TA",this,&opt_glau::baryon_density_eta_envelop_profile_0_function,-10,10,1);
        f1->SetParameter(0, 1.0); 
        double TA = f1->Integral(-10,10,1.0E-09);
        delete f1 ; 
        return TA;
}











