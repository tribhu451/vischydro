#include "mc_glau.h"

using namespace std;

mc_glau::mc_glau(idb *InData1)
{
  IDB = InData1; 
  set_mc_glau_params();
}

mc_glau::~mc_glau()
{
}

void mc_glau::set_ic(grid* f, EoS* eos)
{

  cout << "\n" << endl;
  cout<<"      *************************     "<<endl;
  cout<<"      Monte Carlo Glauber Model     "<<endl;
  cout<<"      *************************     "<<endl;
  cout << "\n" << endl;
  
  cout<<"[Info] "<<IDB->projectile<<"+"<<IDB->target<<" at "<<IDB->SNN<<"GeV( sigma_nn : "<<sigma<<"fm)"<<endl;
  cout<<"[Info] deformation parameter of projectile | Beta-2 : "<<p_beta2<<" and Beta-4 : "<<p_beta4<<endl;
  cout<<"[Info] deformation parameter of target | Beta-2 : "<<t_beta2<<" and Beta-4 : "<<t_beta4<<endl;
  if(X_hard != 1.0)
    {
      cout<<"[Info] two component energy deposition with X_hard = "<<X_hard<<endl;
    }
  cout << "[Info] Gaussian smearing with std deviation : " << DELTA << endl;
  cout << "[Info] Multiplicity scaling factor eps0 : " << eps0 << endl;
  
  
  // random number for b and theta
  TRandom3* t1=new TRandom3();
  t1->SetSeed(0);
  long kss=t1->GetSeed();
  gRandom->SetSeed(kss);
  TF1* f1= new TF1("f1","x",0.0,25.0);
  TF1* f2= new TF1("f2","sin(x)",0.0,TMath::Pi());
 

  // initialization
  double XA[A];double YA[A];double ZA[A];
  double XB[B];double YB[B];double ZB[B];
  double npart_x[500],npart_y[500];
  double ncoll_x[2000],ncoll_y[2000];
  for(int j=0;j<=A;j++){XA[j]=0.0;YA[j]=0.0;ZA[j]=0.0;}
  for(int j=0;j<=B;j++){XB[j]=0.0;YB[j]=0.0;ZB[j]=0.0;}
  for(int j=0;j<500;j++){npart_x[j]=0.0;npart_y[j]=0.0;}
  for(int j=0;j<2000;j++){ncoll_x[j]=0.0;ncoll_y[j]=0.0;}

  // orientation angle
  double p_ori_theta = f2->GetRandom(0.0,TMath::Pi());  
  double t_ori_theta = f2->GetRandom(0.0,TMath::Pi());  
  double p_ori_phi = (2.0*TMath::Pi())*(t1->Rndm());
  double t_ori_phi = (2.0*TMath::Pi())*(t1->Rndm());

  // print orientation
  cout<<"[Info] (projectile orientation) p_theta: "<<p_ori_theta<<" p_phi: "<<p_ori_phi<<endl;
  cout<<"[Info] (target orientation) t_theta: "<<t_ori_theta<<" t_phi: "<<t_ori_phi<<endl;
  
  //generate nucleus
  generate_nucleus(XA,YA,ZA,A,p_radius,p_dlt,p_beta2,p_beta4,p_ori_theta,p_ori_phi);
  generate_nucleus(XB,YB,ZB,B,t_radius,t_dlt,t_beta2,t_beta4,t_ori_theta,t_ori_phi);

  // generate impact parameter
  double b=f1->GetRandom(bmin,bmax);                      
  cout<<"[Info] b = "<<b<<" (fm)"<<endl;

  // generate the shifting angle of one nucleus
  double zhi=(2.0*TMath::Pi())*(t1->Rndm());
  
  //shifting of nucleus 
  shift_nucleus( XA, YA,  ZA, A, +b/2.0, zhi, XA, YA, ZA );
  shift_nucleus( XB, YB,  ZB, B, -b/2.0, zhi, XB, YB, ZB);

 // calculate npat ncoll
  int NPart,NColl; 
  calculate_npart_ncoll(XA,YA,XB,YB,NPart,NColl,npart_x,npart_y,ncoll_x, ncoll_y);                    
  cout<<"[Info] No. of participants : "<<NPart<<endl;
  cout<<"[Info] No. of binary collisions : "<<NColl<<endl;


  std::ofstream File0;
  File0.open("hydro_output/mc_glauber_ic_dist.dat");

  // output can be a input to music
  File0<<"#"<<"\t"<<"mc_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<IDB->neta<<"\t"<<"nx="<<"\t"<<IDB->nx<<"\t"<<"ny="<<"\t"<<IDB->ny
       <<"\t"<<"deta="<<"\t"<<IDB->deta<<"\t"<<"dx="<<"\t"<<IDB->dx<<"\t"<<"dy="<<"\t"<<IDB->dy<<endl;
  
  
  
  double total_deposited_entropy = 0.0;  // total deposited entropy
  double total_deposited_energy = 0.0;  // total deposited energy
  
  cell* c;
  
  for(int i=0; i<IDB->nx; i++)
    {
      for(int j=0; j<IDB->ny; j++)
	{
      
	  double x_ = IDB->xmin + i*IDB->dx;
	  double y_ = IDB->ymin + j*IDB->dy;
	  
	  double entr=0;
	  for(int ks=0; ks<NPart; ks++)
	    {
	      double temp1=TMath::Power(x_-npart_x[ks],2)+TMath::Power(y_-npart_y[ks],2);
	      double value1=0.5*npp*(1-X_hard);
	      entr=entr+((value1/(TMath::Sqrt(2*TMath::Pi()*DELTA*DELTA)))*(TMath::Exp((-1.0/(2.0*DELTA*DELTA))*(temp1))));
	    }
	  
	  for(int ks=0; ks<NColl; ks++)
	    {
	      double temp1=TMath::Power(x_-ncoll_x[ks],2)+TMath::Power(y_-ncoll_y[ks],2);
	      double value1=npp*X_hard;
	      entr=entr+((value1/(TMath::Sqrt(2*TMath::Pi()*DELTA*DELTA)))*(TMath::Exp((-1.0/(2.0*DELTA*DELTA))*(temp1))));
	    }
	  
	  if(entr < 1e-5 ) { entr = 0.0; }
	  
	  
	  for(int k=0; k<IDB->neta; k++)
	    {
	      
	      c = f->get_cell(i,j,k);
	      
	      
	      double eta = IDB->etamin + k*IDB->deta;
	      if ( IDB->neta == 1 ) eta = 0.0 ; // for 2+1D hydro
	      

	      // rapidity distribution 
	      //https://arxiv.org/pdf/0902.4121.pdf  (eqn_2.12)	   
	      double H_eta = exp(  - pow( fabs(eta) - IDB->eta_platue / 2.0, 2 )  /  
				   ( 2 * pow(IDB->eta_fall,2) ) *  theta(fabs(eta)-IDB->eta_platue/2) );
	      
	      
 	      double nb= 0; double nq = 0; double ns =0; 
	      
	      double eps = eos->entr_2_eps(eps0*entr*H_eta,nb,nq,ns);  // entropy converted to energy density
	      
	      double vx=0; double vy=0; double vz= 0;
	      double utau = 1.0 / sqrt( 1.0 - vx*vx - vy*vy - vz*vz );
	      double ux = utau*vx;
	      double uy = utau*vy;
	      double uz = utau*vz; 
	      
	      total_deposited_entropy +=  eps0*entr*H_eta ;
	      total_deposited_energy  += eps ;
	      
	      if(fabs(eta)<0.0001)
		{  // output will be a input to music
		  //cout<<eta<<endl;
		  File0<<eta<<"\t"<<x_<<"\t"<<y_<<"\t"<<entr<<"\t"<<utau<<"\t"<<ux
		       <<"\t"<<uy<<"\t"<<uz<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<endl; // writing the dist in a file
		}
	      
	      
	      c->set_prim_var( eos, IDB->tau0, eps, nb, nq,  ns,  vx,  vy,  vz );
	      
	    }
	}
    }
     cout<<"[Info] total amount of deposited entropy is = "<<total_deposited_entropy<<endl;
     cout<<"[Info] total amount of deposited energy is = "<<total_deposited_energy<<endl;
     cout<<"\n";
}



double mc_glau::theta(double _a)
{
if(_a > 0){return 1;}else{return 0;}
}


void mc_glau::generate_nucleus(double* X1, double* Y1,double* Z1,int A,
                double R, double dlt, double BETA2, double BETA4, double etaA, double psiA)
{    
  double X[500];double Y[500];double Z[500];

  TRandom3* tr1 = new TRandom3();
  tr1->SetSeed(0);

  double CMx=0.0;double CMy=0.0;double CMz=0.0;
  int count=0;
  
  do
    {
      double r=(11.0)*(tr1->Rndm());
      double Theta=(TMath::Pi())*(tr1->Rndm());
      double Phi=((2.0)*TMath::Pi())*(tr1->Rndm());
      double test=tr1->Rndm();
      
      
      double Y20=0.25*TMath::Sqrt(5.0/TMath::Pi())*
	(3*TMath::Cos(Theta)*TMath::Cos(Theta)-1.0);
      double Y40=(3.0/(16.0*TMath::Sqrt(TMath::Pi())))* 
	((35*TMath::Power(TMath::Cos(Theta),4))-
	 (30*TMath::Power(TMath::Cos(Theta),2))+3);
      double RAT= R*(1+(BETA2*Y20)+(BETA4*Y40));
      double rho=(1.0/60.0)*(r*r*(TMath::Sin(Theta)))/(1.0+(TMath::Exp((r-RAT)/dlt)));
      
      
      if(test < rho )
	{      
	  
	  X[count]= (r*TMath::Sin(Theta)*TMath::Cos(Phi));
	  Y[count]=(r*TMath::Sin(Theta)*TMath::Sin(Phi));
	  Z[count]=(r*TMath::Cos(Theta));
	  CMx=CMx+X[count]; CMy=CMy+Y[count] ;CMz=CMz+Z[count];    
	  count=count+1;
	}   
    }   
  while(count<A);
  
  CMx=CMx/A;CMy=CMy/A;CMz=CMz/A;
  
  for(int j=0;j<A;j++){ X[j]=X[j]+(-CMx); Y[j]=Y[j]+(-CMy);Z[j]=Z[j]+(-CMz);}
  
  //etaA - nucleus orientaton angle (theta)
  //psiA - nucleus orientation angle (phi)
  
  for(int j=0;j<A;j++)
    {
      X1[j]=(TMath::Cos(psiA)*TMath::Cos(etaA)*X[j])+(-TMath::Sin(psiA)*Y[j])+(-TMath::Cos(psiA)*TMath::Sin(etaA)*Z[j]);
      Y1[j]=(TMath::Sin(psiA)*TMath::Cos(etaA)*X[j])+(TMath::Cos(psiA)*Y[j])+(-TMath::Sin(psiA)*TMath::Sin(etaA)*Z[j]);
      Z1[j]=(TMath::Sin(etaA)*X[j])+(TMath::Cos(etaA)*Z[j]);
    }
  
}




void mc_glau::calculate_npart_ncoll(double* vxA,double* vyA,double* vxB,double* vyB, int &Npart, 
       int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y)
{
  
   Ncoll=0;
   Npart=0;

  double occA[1000];double occB[1000];         //flag during calc of Npart
  //double Ncoll_x[2000]; double Ncoll_y[2000];  // x & y co-ordinate of binary collision sources
  //double Npart_x[1000]; double Npart_y[1000];  // x & y co-ordinate of participant sources

  
  for(int i=0;i<A;i++){occA[i]=0;}
  for(int i=0;i<B;i++){occB[i]=0;}
  
  for (int i=0; i<A; i++)
    {
      for (int j=0; j<B; j++)
	{  
	  double d=TMath::Sqrt( TMath::Power((vxB[j]-vxA[i]),2) + 
				TMath::Power ( (vyB[j]-vyA[i]),2));
	  double D=TMath::Sqrt(sigma/ (TMath::Pi())); 
	  
	  if( d <= D)
	    { 
	      Ncoll_x[Ncoll]=(vxA[i]+vxB[j])/2;
	      Ncoll_y[Ncoll]=(vyA[i]+vyB[j])/2;
	      Ncoll=Ncoll+1;
	      
	      if(occA[i]==0)
		{ 
		  occA[i]=1;Npart_x[Npart]=vxA[i]; Npart_y[Npart]=vyA[i];Npart=Npart+1;
		} 
	      
              
	      if(occB[j]==0)
		{
		  occB[j]=1;Npart_x[Npart]=vxB[j]; Npart_y[Npart]=vyB[j];Npart=Npart+1;
		}
	      
	      
	    }                                                           
	}                                                          
    }                                                         
      // [Info]  shifting the energy distributions center to (0,0,0)   
      double xref1=0.0;
      double yref1=0.0;
      double wref1=0.0;
      for(int k=0;k<Npart;k++)
	{ 
	  xref1=xref1+(Npart_x[k]*(0.5*npp*X_hard));
	  yref1=yref1+(Npart_y[k]*(0.5*npp*X_hard));
	  wref1=wref1+(0.5*npp*X_hard);
	}
      
      double xref2=0.0;
      double yref2=0.0;
      double wref2=0.0;
      for(int k=0;k<Ncoll;k++)
	{ 
	  xref2=xref2+(Ncoll_x[k]*(npp*(1-X_hard)));
	  yref2=yref2+(Ncoll_y[k]*(npp*(1-X_hard)));
	  wref2=wref2+(npp*(1-X_hard));
	}
      
      double xAverage=((xref1+xref2)/(wref1+wref2));
      double yAverage=((yref1+yref2)/(wref1+wref2));
      
      // cout<<xAverage<<"  "<<yAverage<<"\n";
      
      for(int k=0;k<Npart;k++)
	{ 
          Npart_x[k] = Npart_x[k]-xAverage;
          Npart_y[k] = Npart_y[k]-yAverage;
	}
     for(int k=0;k<Ncoll;k++)
	{ 
          Ncoll_x[k] = Ncoll_x[k]-xAverage;
          Ncoll_y[k] = Ncoll_y[k]-yAverage;
	}

}

void mc_glau::shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
                     double zhi,double* X2, double* Y2, double* Z2 )
{
  for(int j=0;j<A;j++){ X2[j]=X1[j]+((b)*TMath::Cos(zhi)); Y2[j]=Y1[j]+((b)*TMath::Sin(zhi));}
}


