#include "fluid.h"


fluid::fluid(EoS* _eos, trancoeff* _tr,InputData* _ind)
{
  DATA = _ind;
  eos = _eos;

  double eCrit = eos->temp_2_eps(DATA->Tfreeze,0,0,0);
  trcoef= _tr;

  dt = DATA->dtau;
  nx = DATA->nx;
  ny = DATA->ny;
  nz = DATA->neta;
  xmin = DATA->xmin;
  xmax = DATA->xmax;
  ymin = DATA->ymin;
  ymax = DATA->ymax;
  zmin = DATA->etamin;
  zmax = DATA->etamax;
  dx = (xmax - xmin) / (nx - 1);
  dy = (ymax - ymin) / (ny - 1);
  if(nz==1){dz = zmax-zmin;} else {dz = (zmax - zmin) / (nz - 1);}


  if (nz > 1)
    {
     if(dt > dx/2. || dt > dy/2. || dt > dz/2.){cout<<"[Error]  too big delta_tau"<<endl; exit(1);}
     }
 else
    { 
      if(dt > dx/2. || dt > dy/2.){cout<<"[Error]  too big delta_tau"<<endl; exit(1);}
    }

  
  Cell = new cell[nx * ny * nz];

  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++)
	{
	  get_cell(ix, iy, iz)->set_prev_cell(X_, get_cell(ix - 1, iy, iz));
	  get_cell(ix, iy, iz)->set_next_cell(X_, get_cell(ix + 1, iy, iz));
	  get_cell(ix, iy, iz)->set_prev_cell(Y_, get_cell(ix, iy - 1, iz));
	  get_cell(ix, iy, iz)->set_next_cell(Y_, get_cell(ix, iy + 1, iz));
	  get_cell(ix, iy, iz)->set_prev_cell(Z_, get_cell(ix, iy, iz - 1));
	  get_cell(ix, iy, iz)->set_next_cell(Z_, get_cell(ix, iy, iz + 1));
	  get_cell(ix, iy, iz)->set_pos(ix, iy, iz);
	}


  cout<<"[Info] minimum X : "<<get_xmin()<<"\t   maximum X : "<<get_xmax()<<"\t   nx : "<<get_nx()<<"\t   dx : "<<get_dx()<<endl;
  cout<<"[Info] minimum Y : "<<get_ymin()<<"\t   maximum Y : "<<get_ymax()<<"\t   ny : "<<get_ny()<<"\t   dy : "<<get_dy()<<endl;
  cout<<"[Info] minimum Z : "<<get_zmin()<<"\t   maximum Z : "<<get_zmax()<<"\t   neta : "<<get_nz()<<"\t   deta : "<<get_dz()<<endl;
  cout<<"[Info] dt : "<<dt<<" (fm)"<<endl;
  cout<<"[Info] dt/dx : "<<dt/get_dx()<<"\tdt/dy : "<<dt/get_dy()<<"\tdt/deta : "<<dt/get_dz()<<"\n"<<endl;
    

  //----- Cornelius init
  double arrayDx[3] = {dt, dx, dy};
  cornelius = new Cornelius;
  cornelius->init(3, eCrit, arrayDx);
  
  freeze_file.open("./hydro_output/surface.dat",ios::trunc | ios::out );
freeze_file<<std::setprecision(6)<<std::scientific;
}


fluid::~fluid(){ delete[] Cell;  freeze_file.close();}


void fluid::getCMFvariablesUmunu(cell *c, double tau, double &e, double &nb, 
				 double &nq, double &ns, double &vx, double &vy, double &Y){
  double p; double vz;
  c->get_physical_var(eos, tau, e, p, nb, nq, ns,vx, vy, vz);
  double eta = get_z(c->get_iz());
  Y = eta + 0.5*log((1.+vz)/(1.-vz));
  double gf=1./sqrt(1.-vx*vx-vy*vy-vz*vz);
  vx = vx*gf;
  vy = vy*gf;
}
           //check and use


void fluid::correct_imaginary_cells()
{

  double Q[7]={0.0};
  // Z
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   // left boundary
   get_cell(ix, iy, 2)->get_Q(Q);
   get_cell(ix, iy, 1)->set_Q(Q);
   get_cell(ix, iy, 0)->set_Q(Q);
   // right boundary
   get_cell(ix, iy, nz - 3)->get_Q(Q);
   get_cell(ix, iy, nz - 2)->set_Q(Q);
   get_cell(ix, iy, nz - 1)->set_Q(Q);
  }
 // Y
 for (int ix = 0; ix < nx; ix++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   get_cell(ix, 2, iz)->get_Q(Q);
   get_cell(ix, 1, iz)->set_Q(Q);
   get_cell(ix, 0, iz)->set_Q(Q);
   // right boundary
   get_cell(ix, ny - 3, iz)->get_Q(Q);
   get_cell(ix, ny - 2, iz)->set_Q(Q);
   get_cell(ix, ny - 1, iz)->set_Q(Q);
  }
 // X
 for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   get_cell(2, iy, iz)->get_Q(Q);
   get_cell(1, iy, iz)->set_Q(Q);
   get_cell(0, iy, iz)->set_Q(Q);
   // right boundary
   get_cell(nx - 3, iy, iz)->get_Q(Q);
   get_cell(nx - 2, iy, iz)->set_Q(Q);
   get_cell(nx - 1, iy, iz)->set_Q(Q);
  }
 
}



void fluid::correct_imaginary_cells_full() {
 double Q[7], _pi[4][4], _Pi;
 // Z
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   // left boundary
   get_cell(ix, iy, 2)->get_Q(Q);
   get_cell(ix, iy, 1)->set_Q(Q);
   get_cell(ix, iy, 0)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = get_cell(ix, iy, 2)->get_pi(i, j);
   _Pi = get_cell(ix, iy, 2)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(ix, iy, 0)->set_pi(i, j, _pi[i][j]);
     get_cell(ix, iy, 1)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(ix, iy, 0)->set_Pi(_Pi);
   get_cell(ix, iy, 1)->set_Pi(_Pi);
   // right boundary
   get_cell(ix, iy, nz - 3)->get_Q(Q);
   get_cell(ix, iy, nz - 2)->set_Q(Q);
   get_cell(ix, iy, nz - 1)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = get_cell(ix, iy, nz - 3)->get_pi(i, j);
   _Pi = get_cell(ix, iy, nz - 3)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(ix, iy, nz - 2)->set_pi(i, j, _pi[i][j]);
     get_cell(ix, iy, nz - 1)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(ix, iy, nz - 2)->set_Pi(_Pi);
   get_cell(ix, iy, nz - 1)->set_Pi(_Pi);
  }
 // Y
 for (int ix = 0; ix < nx; ix++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   get_cell(ix, 2, iz)->get_Q(Q);
   get_cell(ix, 1, iz)->set_Q(Q);
   get_cell(ix, 0, iz)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = get_cell(ix, 2, iz)->get_pi(i, j);
   _Pi = get_cell(ix, 2, iz)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(ix, 0, iz)->set_pi(i, j, _pi[i][j]);
     get_cell(ix, 1, iz)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(ix, 0, iz)->set_Pi(_Pi);
   get_cell(ix, 1, iz)->set_Pi(_Pi);
   // right boundary
   get_cell(ix, ny - 3, iz)->get_Q(Q);
   get_cell(ix, ny - 2, iz)->set_Q(Q);
   get_cell(ix, ny - 1, iz)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = get_cell(ix, ny - 3, iz)->get_pi(i, j);
   _Pi = get_cell(ix, ny - 3, iz)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(ix, ny - 2, iz)->set_pi(i, j, _pi[i][j]);
     get_cell(ix, ny - 1, iz)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(ix, ny - 2, iz)->set_Pi(_Pi);
   get_cell(ix, ny - 1, iz)->set_Pi(_Pi);
  }
 // X
 for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   get_cell(2, iy, iz)->get_Q(Q);
   get_cell(1, iy, iz)->set_Q(Q);
   get_cell(0, iy, iz)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = get_cell(2, iy, iz)->get_pi(i, j);
   _Pi = get_cell(2, iy, iz)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(0, iy, iz)->set_pi(i, j, _pi[i][j]);
     get_cell(1, iy, iz)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(0, iy, iz)->set_Pi(_Pi);
   get_cell(1, iy, iz)->set_Pi(_Pi);
   // right boundary
   get_cell(nx - 3, iy, iz)->get_Q(Q);
   get_cell(nx - 2, iy, iz)->set_Q(Q);
   get_cell(nx - 1, iy, iz)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = get_cell(nx - 3, iy, iz)->get_pi(i, j);
   _Pi = get_cell(nx - 3, iy, iz)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(nx - 2, iy, iz)->set_pi(i, j, _pi[i][j]);
     get_cell(nx - 1, iy, iz)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(nx - 2, iy, iz)->set_Pi(_Pi);
   get_cell(nx - 1, iy, iz)->set_Pi(_Pi);
  }
}



/*
void fluid::get_surface3d(double tau)
{
  // cornelius corner poits
  double ****ccube = new double ***[2];
  for (int i1 = 0; i1 < 2; i1++) {
    ccube[i1] = new double **[2];
    for (int i2 = 0; i2 < 2; i2++) {
      ccube[i1][i2] = new double *[2];
      for (int i3 = 0; i3 < 2; i3++) {
	ccube[i1][i2][i3] = new double[2];
      }
    }
  }

 for (int ix = 2; ix < nx - 2; ix++)
   for (int iy = 2; iy < ny - 2; iy++)
     for (int iz = 2; iz < nz - 2; iz++)
       {
	 
	 for (int jx = 0; jx < 2; jx++)
	   for (int jy = 0; jy < 2; jy++)
	     for (int jz = 0; jz < 2; jz++) 
	       {
		 double _p, _nb, _nq, _ns, _vx, _vy, _vz;
		 cell *cc = getCell(ix + jx, iy + jy, iz + jz);
		 cc->get_physical_var(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
		 ccube[1][jx][jy][jz] = e;
		 cc->get_center_var_prev(eos, tau - dt, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
		 ccube[0][jx][jy][jz] = e;
	       } // jz loop
	 
	 
	 cornelius->find_surface_4d(ccube);
	 const int Nsegm = cornelius->get_Nelements();
	 for (int isegm = 0; isegm < Nsegm; isegm++)
	   {
	     freeze_file.precision(15);
	     freeze_file  << tau + cornelius->get_centroid_elem(isegm, 0)
			     << "\t" << getX(ix) + cornelius->get_centroid_elem(isegm, 1)
			     << "\t" << getY(iy) + cornelius->get_centroid_elem(isegm, 2)
			     << "\t" << getZ(iz) + cornelius->get_centroid_elem(isegm, 3)<<endl;
	     
	   }
 
       } // iz loop
 
}
*/

void fluid::get_surface2d(double tau)
{
  double ***ccube = new double**[2];
  for (int i1=0; i1 < 2; i1++) {
    ccube[i1] = new double*[2];
    for (int i2=0; i2 < 2; i2++) {
      ccube[i1][i2] = new double[2];
    }
  }


 for (int ix = 2; ix < nx - 2; ix++)
   for (int iy = 2; iy < ny - 2; iy++)
       { // ix,iy loop
	 
         double Qcube[2][2][2][7] = {0.0}; // tau, x, y and 7-conserved quantities.
         double picube[2][2][10] = {0.0}; //  x, y and 10-components.
         double Picube[2][2] = {0.0}; //  x and y only

         cell* c;
	 for (int jx = 0; jx < 2; jx++)
	   for (int jy = 0; jy < 2; jy++)
	       {
		 double e,_p, _nb, _nq, _ns, _vx, _vy, _vz;
		 c = get_cell(ix + jx, iy + jy, 0); // 0 for eta
		 c->get_physical_var(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
		 ccube[1][jx][jy] = e;
                 c->get_Q(Qcube[1][jx][jy]);
		 c->get_center_var_prev(eos, tau - dt, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
                 c->get_Q(Qcube[0][jx][jy]);
		 ccube[0][jx][jy] = e;

	         for (int ii = 0; ii < 4; ii++)
		  for (int jj = 0; jj <= ii; jj++){
                    int anum = c->indexpi(ii,jj) ;
		    picube[jx][jy][anum] = c->get_pi(ii, jj);}
	         Picube[jx][jy] = c->get_Pi();
	       }
	 
	 
	 cornelius->find_surface_3d(ccube);
	 const int Nsegm = cornelius->get_Nelements();
	 for (int isegm = 0; isegm < Nsegm; isegm++)
	   { // nelement loop
             
             double etad = 0.0 ; // for 2+1D it's always zero 
             double td = tau - dt + cornelius->get_centroid_elem(isegm, 0) ;
             double xd = get_x(ix) + cornelius->get_centroid_elem(isegm, 1) ;
             double yd = get_y(iy) + cornelius->get_centroid_elem(isegm, 2) ;
	     freeze_file << td << "\t" << xd << "\t" << yd << "\t" << etad << "\t" ;   //  [cout-1]

            // Linear interpolation below
             double wt[2] = { 1. - cornelius->get_centroid_elem(isegm, 0) / dt,
			      cornelius->get_centroid_elem(isegm, 0) / dt } ;
             double wx[2] = { 1. - cornelius->get_centroid_elem(isegm, 1) / dx,
			      cornelius->get_centroid_elem(isegm, 1) / dx} ;
             double wy[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / dy,
			      cornelius->get_centroid_elem(isegm, 2) / dy } ;

           double Qc[7] = {0.0} ; // value at center
           double ec,pc,nbc,nqc,nsc,vxc,vyc,vzc,Tc ; // value at center
	   for (int jt = 0; jt < 2; jt++)
	     for (int jx = 0; jx < 2; jx++)
	       for (int jy = 0; jy < 2; jy++)
		   for (int i = 0; i < 7; i++)
		     {
		       Qc[i] += Qcube[jt][jx][jy][i] * wt[jt] * wx[jx] * wy[jy];
		     }

	   for (int i = 0; i < 7; i++) { Qc[i] = Qc[i] / td ; }

	   CALC_2_LRF(eos, Qc, ec, pc, nbc, nqc, nsc , vxc, vyc, vzc);
           Tc = eos->temperature(ec,nbc,nqc,nsc);

	   freeze_file << ec << "\t" << Tc <<"\t";                                     //  [cout-2]

           if ( Tc > 0.175 || Tc < 0.125 )
              {cout<<"Wrong interpolation in hypersurface"
                   <<" !!! Tc is "<<Tc<<" // ec is "<<ec<< " // nbc  = "<<nbc<<endl;exit(1);
                 }

	   double v2 = vxc * vxc + vyc * vyc + vzc * vzc;
	   if (v2 > 1.)
	     {
	       vxc *= sqrt(0.999 / v2);
	       vyc *= sqrt(0.999 / v2);
	       vzc *= sqrt(0.999 / v2);
	       v2 = 0.999;
	     }

	   freeze_file << vxc << "\t" << vyc << "\t" << vzc << "\t" ;                //  [cout-3]
  
           double pic[10] = {0.0} ; // value at center
           double Pic = 0.0 ;

	   for (int jx = 0; jx < 2; jx++)
	     for (int jy = 0; jy < 2; jy++)
		 {
		   for (int ii = 0; ii < 10; ii++)
		     pic[ii] += picube[jx][jy][ii] * wx[jx] * wy[jy] ;

		   Pic += Picube[jx][jy] * wx[jx] * wy[jy] ;
		 }



	         for (int ii = 0; ii < 4; ii++)
		  for (int jj = 0; jj <= ii; jj++){
                    int anum = c->indexpi(ii,jj) ;
                     freeze_file << pic[anum] << "\t" ;  }                         //  [cout-4]
  
          freeze_file << Pic << endl ;                                             //  [cout-5]

	   } // nelement loop
 
       } // ix,iy loop
 
}


