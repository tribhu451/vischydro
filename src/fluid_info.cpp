#include "fluid_info.h"

fluid_info::fluid_info(grid* _f, EoS* _eos, idb *_IDB)
{
  f = _f; eos = _eos ; IDB = _IDB;
}

fluid_info::~fluid_info(){}

void fluid_info::midslice_xy_dist(double time)
{
      
      cout<<"[Info] after time "<<time<<" data is recorded"<<endl;
      char name[200];
      
      FILE* File2;
      FILE* File3;
      
      sprintf(name,"hydro_output/midslice_xy_dist_at_%f_fm.dat",time);	
      File2 = fopen(name, "w");
      
      sprintf(name,"hydro_output/rapidity_dist_at_%f_fm.dat",time);	
      File3 = fopen(name, "w");
      
      double erg,vx,vy,vz,nb,nq,ns,pressure;
      
      for(int i=0; i<IDB->nx; i++)
        for(int j=0; j<IDB->ny; j++)
          for(int k=0; k<IDB->neta; k++)
	    {
	      double x_ = IDB->xmin + i*IDB->dx;
	      double y_ = IDB->ymin + j*IDB->dy;
              double eta_ = IDB->etamin + k*IDB->deta;
	      double rt = TMath::Sqrt(x_*x_+y_*y_);
	      f->get_cell(i,j,k)->get_physical_var(eos, time, erg, pressure, nb, nq, ns, vx,vy, vz); 
              double pi00 = f->get_cell(i,j,k)->get_pi(0,0);         
              double pi10 = f->get_cell(i,j,k)->get_pi(1,0);         
              double pi20 = f->get_cell(i,j,k)->get_pi(2,0);         
              double pi30 = f->get_cell(i,j,k)->get_pi(3,0)*time;         
              double pi11 = f->get_cell(i,j,k)->get_pi(1,1);         
              double pi21 = f->get_cell(i,j,k)->get_pi(2,1);         
              double pi31 = f->get_cell(i,j,k)->get_pi(3,1)*time;         
              double pi22 = f->get_cell(i,j,k)->get_pi(2,2);         
              double pi32 = f->get_cell(i,j,k)->get_pi(3,2)*time;         
              double pi33 = f->get_cell(i,j,k)->get_pi(3,3)*pow(time,2);         
              double Pi = f->get_cell(i,j,k)->get_Pi();         
	      double T = eos->temperature(erg,nb, nq,ns);
              double s = eos->entropy(erg,nb,nq,ns);

	      if( IDB->neta==1 || abs(eta_)<0.00001 )  
                 {
                   fprintf(File2,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
			   x_, y_,rt,eta_,erg,vx,vy,vz,pi00,pi10,pi20,pi30,pi11,pi21,pi31,pi22,pi32,pi33,Pi,T,s); 
		 }
	      
	      if( abs(x_) < 0.00001 && abs(y_) < 0.00001 )  
		{
		  fprintf(File3,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
			  x_, y_,rt,eta_,erg,vx,vy,vz,pi00,pi10,pi20,pi30,pi11,pi21,pi31,pi22,pi32,pi33,Pi,T,s); 
                }
	    }
      
      fclose(File2);
      fclose(File3);
        
}


void fluid_info::anisotropy_in_xy_plane(double time)
{
      // set up file name
      const string out_name_ani = "hydro_output/anisotropy_in_xy_plane.dat";
      string out_open_mode_ani;
      FILE *out_file_ani;
      
      
      // If it's the first timestep, overwrite the previous file
      if (time - IDB->tau0 < 0.0001) {
	out_open_mode_ani = "w";
      }
      else {
	out_open_mode_ani = "a";	
      }
      
      out_file_ani =
	fopen(out_name_ani.c_str(), out_open_mode_ani.c_str());
      
      
      
      cout<<"[Info] anisotropy stored at "<<time<<" fm."<<endl;
      double erg,vx,vy,vz,nb,nq,ns,pressure;
      double y2mx2_avg =0.0;     // y^{2} - y^{2} average
      double y2px2_avg =0.0;     // y^{2} + x^{2} average
      double txxmtyy =0.0;       // T_{id}^{xx} - T_{id}^{yy}
      double txxptyy =0.0;       // T_{id}^{xx} + T_{id}^{yy}
      
      double TXXmTYY =0.0;
      double TXXpTYY = 0.0;
      
      int nx  = IDB->nx;
      int ny  = IDB->ny;
      int nz  = IDB->neta;
      
      
      for(int i=0; i<nx ; i++)
	for(int j=0; j<ny; j++)
	  for(int k=(nz/2); k<(nz/2)+1; k++)  // only at the mid rapidity slice
	    {
	      double x_ = IDB->xmin + i*IDB->dx;
	      double y_ = IDB->ymin + j*IDB->dy;
	      //double z_ = f->get_z(k);
	      f->get_cell(i,j,k)->get_physical_var(eos, time, erg, pressure, nb, nq, ns, vx,vy, vz);          
	      y2mx2_avg = y2mx2_avg + ((y_*y_ - x_*x_)*erg);
	      y2px2_avg = y2px2_avg + ((y_*y_ + x_*x_)*erg);
	      
	      double gamma2 = 1./ (1.-vx*vx-vy*vy-vz*vz);
	      txxmtyy = txxmtyy + ((erg+pressure)*gamma2*(vx*vx-vy*vy)); // T_{id}^{xx} - T_{id}^{yy}
	      txxptyy = txxptyy + (((erg+pressure)*gamma2*(vx*vx+vy*vy))+(2*pressure)); // T_{id}^{xx} + T_{id}^{yy}
	      
	      double Pi = f->get_cell(i,j,k)->get_Pi();
	      double pixx = f->get_cell(i,j,k)->get_pi(1,1);
	      double piyy = f->get_cell(i,j,k)->get_pi(2,2);
	      
	      TXXmTYY +=  ( (erg+pressure+Pi)*gamma2*(vx*vx-vy*vy) + (pixx-piyy) ); // T_{vis}^{xx} - T_{vis}^{yy}
	      TXXpTYY +=  ( (erg+pressure+Pi)*gamma2*(vx*vx+vy*vy) + (pixx+piyy) + (2.0*pressure) + (2.0*Pi) ); // T_{vis}^{xx} + T_{vis}^{yy}
	      
	    }
      double ex =  y2mx2_avg/ y2px2_avg;
      double ep = txxmtyy/ txxptyy;
      double epsp = TXXmTYY/ TXXpTYY;
      
      
      fprintf(out_file_ani,"%e %e %e %e\n",
	      time-IDB->tau0,ex,ep,epsp); 
      
      
      fclose(out_file_ani);
        
}


void fluid_info::conservation_check(double tau)
{

      // set up file name
      const string out_name_cons = "hydro_output/total_conserved_quantity.dat";
      string out_open_mode_cons;
      FILE *out_file_cons;
      
      
      // If it's the first timestep, overwrite the previous file
      if (tau - IDB->tau0 < 0.0001) {
	out_open_mode_cons = "w";
      }
      else {
	out_open_mode_cons = "a";	
      }
      
      out_file_cons =
	fopen(out_name_cons.c_str(), out_open_mode_cons.c_str());
      
      cout<<"[Info] total energy and entropy stored at "<<tau<<" fm."<<endl;
      
      double e_total = 0.0; // total energy
      double s_total = 0.0; // total entropy
      
      double eps,p,nb,nq,ns,vx,vy,vz;
      
      for(int ix=0; ix<IDB->nx; ix++)
	for(int iy=0; iy<IDB->ny; iy++)
	  for(int iz=0; iz<IDB->neta; iz++){
	    
	    double eta = IDB->etamin + iz*IDB->deta;
	    if(IDB->neta == 1 ) {eta = 0.0;}
	    f->get_cell(ix,iy,iz)->get_physical_var(eos, tau, eps, p, nb, nq, ns, vx,vy, vz);
	    double pitt = f->get_cell(ix,iy,iz)->get_pi(0,0);
	    double pite = f->get_cell(ix,iy,iz)->get_pi(0,3);
            double Pi = f->get_cell(ix,iy,iz)->get_Pi();         
	    double s = eos->entropy(eps,nb,nq,ns);
	    
	    double utau = 1.0 / sqrt( 1.0 - vx*vx - vy*vy - vz*vz ) ;
	    e_total +=  tau * (  
                       cosh(eta) * ( ( eps + p + Pi ) * pow( utau , 2.0 ) - p - Pi + pitt )
                       + sinh(eta) * ( ( eps + p + Pi ) * pow( utau , 2.0 ) * vz + tau * pite )
                              ) ;
	    s_total += tau*s*utau ;
	    
	  } // for loop of ix,iy,iz end .
      
      
      fprintf(out_file_cons,"%e %e %e\n",
	      tau-IDB->tau0,e_total,s_total); 
      
      fclose(out_file_cons);
        
}


