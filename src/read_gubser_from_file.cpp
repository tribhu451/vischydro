#include "read_gubser_from_file.h"

using namespace std;

read_gubser_from_file::read_gubser_from_file(InputData *InData1)
{
  InData = InData1; 
}



read_gubser_from_file::~read_gubser_from_file(){}



void read_gubser_from_file::set_ic(fluid* f, EoS* eos, double tau_0)
{
  cout<<"*** [IC] ***"<<endl;
  cout<<"[Info] reading from file for gubser test ... \n"<<endl;
  double x,y,eps,ux,uy,pixx,piyy;
  double pixy,pitt,pitx,pity,piee;
  
  cell* c;
  fstream ic_file;
  ic_file.open("tests/gubser/Initial_Profile.dat");
  if (!ic_file) 
    {cout<<"couldn't find ic file."<<endl; exit(1);}
  
  while(!ic_file.eof())
    { // while start


  ic_file.getline(buff,200);
  iss = new istringstream(buff);
  *iss>>x>>y>>eps>>ux>>uy>>pixx>>piyy>>pixy>>pitt>>pitx>>pity>>piee;
  delete iss;
      
      for(int i=0; i<f->get_nx(); i++)
	{ // i loop
	  double x_ = f->get_x(i);
	  if (abs(x-x_)<0.001)
	    { // if x
	      for(int j=0; j<f->get_ny(); j++)
		{
		  double y_ = f->get_y(j);
		  if(abs(y - y_)<0.001)
		    { //if y
		      
		      c = f->get_cell(i,j,0);         
		      double q =1.0;
		      double rt = TMath::Sqrt(x_*x_+y_*y_);
		      double k_ = TMath::ATanH((2.0*q*q*tau_0*rt)/(1+q*q*tau_0*tau_0+q*q*rt*rt));
		      double vx =  (x_/rt)*(TMath::TanH(k_)); 
		      double vy =  (y_/rt)*(TMath::TanH(k_)); 
		      double vz= 0.0; double nb= 0; double ns=0; double nq=0;
		      c->set_prim_var(eos,tau_0,eps/5.068, nb, nq,  ns,  vx,  vy,  vz);	      
                      c->set_pi(1,1, pixx/5.068);
                      c->set_pi(2,2, piyy/5.068);
                      c->set_pi(2,1, pixy/5.068);
                      c->set_pi(1,2, pixy/5.068);
                      c->set_pi(0,0, pitt/5.068);
                      c->set_pi(0,1, pitx/5.068);
                      c->set_pi(1,0, pitx/5.068);
                      c->set_pi(0,2, pity/5.068);
                      c->set_pi(2,0, pity/5.068);
                      c->set_pi(3,3, piee/5.068);
		      
		    } // if y
		  else
		    {
		      continue ; 
		    }
		  
		} 
	    } // if x
          else
	    {
	      continue ;
	    }
	} // i loop
    }//while end
 ic_file.close();
  cout<<"\n";
}




