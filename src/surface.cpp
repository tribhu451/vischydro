#include "surface.h"

surf::surf(EoS* _eos,idb* _IDB, grid* _f, cnvrt* _CN)
{
  eos = _eos;
  IDB = _IDB ;
  f = _f ;
  CN = _CN ;

  double eCrit;

  if (IDB->eps_freeze_flag == 0 ){
    cout << "[Info] Constant temperature hypersurface." << endl ; 
    cout<<"[Info] Freezeout temperature = "<<IDB->Tfreeze<<" GeV"<<endl;
    eCrit=eos->temp_2_eps(IDB->Tfreeze,0,0,0);
    cout<<"[Info] Freeze-out hypersurface at Ef = "<<eCrit<<" GeV/fm^3\n"<<endl;
  }
 else{
    cout << "[Info] Constant energy density hypersurface." << endl ;
    cout << "[Info] Freezeout energy density = " << IDB->eps_freeze <<" GeV." << endl; 
    eCrit = IDB->eps_freeze ;
    cout << "[Info] Freeze-out hypersurface at Ef = " << eCrit << " GeV/fm^3\n" << endl;
 }
  
  //----- Cornelius init
  
  if(IDB->neta == 1)
    {
      double arrayDx[3] = {IDB->skip_fo_tau*f->get_dtau(),IDB->skip_fo_x* f->get_dx(),IDB->skip_fo_y* f->get_dy()};
      cornelius = new Cornelius;
      cornelius->init(3, eCrit, arrayDx);   
      freeze_file.open("./hydro_output/surface.dat", ios::out );
      freeze_file<<std::setprecision(10)<<std::scientific;
    }
  else
   {
     double arrayDx[4] = {IDB->skip_fo_tau*f->get_dtau(), IDB->skip_fo_x*f->get_dx(), IDB->skip_fo_y*f->get_dy(), IDB->skip_fo_eta*f->get_deta()};
     cornelius = new Cornelius;
     cornelius->init(4, eCrit, arrayDx);
     freeze_file.open("./hydro_output/surface.dat", ios::out );
     freeze_file<<std::setprecision(10)<<std::scientific;
   }
  
}


surf::~surf()
{
  freeze_file.close();
  delete cornelius;
}


void surf::get_surface2p1d(double tau)
{

  double ***ccube = new double**[2];
  for (int i1=0; i1 < 2; i1++)
    {
      ccube[i1] = new double*[2];
      for (int i2=0; i2 < 2; i2++)
	{
	  ccube[i1][i2] = new double[2];
	}
    }
  

  double _e,_p, _nb, _nq, _ns, _vx, _vy, _vz;  
  double Qcube[2][2][2][7] = {0.0}; // tau, x, y and 7-conserved quantities.
  double picube[2][2][2][10] = {0.0}; //  tau ,x, y and 10-components.
  double Picube[2][2][2] = {0.0}; //   tau, x and y .
  
  cell* c;
  
  for (int ix = 2; ix < f->get_nx() - 2; ix += IDB->skip_fo_x )
    for (int iy = 2; iy < f->get_ny() - 2; iy += IDB->skip_fo_y )
      { // ix,iy loop
	
	// put the values in 8 corner points of the cube
        // (1)
	c = f->get_cell(ix, iy , 0); // 0 for eta
        c->get_Q(Qcube[1][0][0]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][0][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][0][0][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][0][0] = c->get_Pi();
	
	// (2)
	c = f->get_cell(ix +  IDB->skip_fo_x, iy , 0);
        c->get_Q(Qcube[1][1][0]); 
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][1][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][1][0][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][1][0] = c->get_Pi();
	
        // (3)
	c = f->get_cell(ix , iy +  IDB->skip_fo_y  , 0); 
        c->get_Q(Qcube[1][0][1]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][0][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][0][1][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][0][1] = c->get_Pi();
	
	// (4)	
	c = f->get_cell(ix +  IDB->skip_fo_x  , iy +  IDB->skip_fo_y  , 0); 
        c->get_Q(Qcube[1][1][1]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][1][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][1][1][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][1][1] = c->get_Pi();
	
	//(5)
	c = f->get_cell(ix, iy , 0); 
        c->get_Q_fo_prev(Qcube[0][0][0]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][0][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][0][0][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][0][0] = c->get_Pi_fo_prev();
	
	// (6)	
	c = f->get_cell(ix +  IDB->skip_fo_x , iy , 0);
        c->get_Q_fo_prev(Qcube[0][1][0]); 
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][1][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][1][0][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][1][0] = c->get_Pi_fo_prev();
	
	// (7)	
	c = f->get_cell(ix  , iy +  IDB->skip_fo_y, 0); 
        c->get_Q_fo_prev(Qcube[0][0][1]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][0][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][0][1][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][0][1] = c->get_Pi_fo_prev();
	
	// (8)	
	c = f->get_cell(ix +  IDB->skip_fo_x , iy +  IDB->skip_fo_y, 0); 
        c->get_Q_fo_prev(Qcube[0][1][1]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][1][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][1][1][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][1][1] = c->get_Pi_fo_prev();
	
	// cornelius magic	 
	cornelius->find_surface_3d(ccube);
	
	for (int isegm = 0; isegm < cornelius->get_Nelements(); isegm++)
	  { // element loop
	    
	    if (ix <= 4 || ix >= f->get_nx() - 4
		|| iy <= 4 || iy >= f->get_ny() - 4 )
	      {
		cout << "[Error] Freeze-out cell at the boundary !!! "
		     << "The grid is too small !!! exit(1)" << endl;
		exit(1);
	      }
	    
	    // fo position
	    double tau_fo = tau - IDB->skip_fo_tau*f->get_dtau() + cornelius->get_centroid_elem(isegm, 0) ;
	    double x_fo = f->get_x(ix) + cornelius->get_centroid_elem(isegm, 1) ;
	    double y_fo = f->get_y(iy) + cornelius->get_centroid_elem(isegm, 2) ;
	    double eta_fo = 0.0 ;   // for 2+1D it's zero 

	    
	    // normal vector
	    double norm_tau = cornelius->get_normal_elem(isegm, 0);
	    double norm_x = cornelius->get_normal_elem(isegm, 1);
	    double norm_y = cornelius->get_normal_elem(isegm, 2);
	    double norm_eta =  0.0 ; // for 2+1D it's zero

	    
	    // Linear interpolation below
	    double wt[2] = { 1. - cornelius->get_centroid_elem(isegm, 0) / IDB->skip_fo_tau*f->get_dtau(),
			     cornelius->get_centroid_elem(isegm, 0) / IDB->skip_fo_tau*f->get_dtau() } ;
	    double wx[2] = { 1. - cornelius->get_centroid_elem(isegm, 1) / (f->get_dx()* IDB->skip_fo_x),
			     cornelius->get_centroid_elem(isegm, 1) / (f->get_dx()* IDB->skip_fo_x)} ;
	    double wy[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / (f->get_dy()* IDB->skip_fo_y),
			    cornelius->get_centroid_elem(isegm, 2) / (f->get_dy()* IDB->skip_fo_y) } ;
	    
	    double Qc[7] = {0.0} ; // value at C.M. of the element
	    
	    for (int jt = 0; jt < 2; jt++)
	      for (int jx = 0; jx < 2; jx++)
		for (int jy = 0; jy < 2; jy++)
		  for (int i = 0; i < 7; i++)
		    {
		      Qc[i] += Qcube[jt][jx][jy][i] * wt[jt] * wx[jx] * wy[jy];
		    }
	    
	    // e_fo, T_fo, nb_fo, nq_fo, ns_fo
	    for (int i = 0; i < 7; i++) { Qc[i] = Qc[i] / tau_fo ; } 
	    CN->CALC_2_LRF(eos, Qc, _e, _p, _nb, _nq, _ns , _vx, _vy, _vz);
	    double T_fo = eos->temperature(_e,_nb,_nq,_ns);
	    
	    
	    if ( T_fo < IDB->Tfreeze - 0.010  || T_fo > IDB->Tfreeze + 0.010 )
	      {
                cout << "[Error] Wrong interpolation in hypersurface"
		     << " !!! T_fo = " << T_fo <<", e_fo = " << _e
                     << ", nb_fo  = " << _nb << endl;
                cout << "Actual Tfreeze = " << IDB->Tfreeze << endl;
		exit(1);
	      }
	  

            double mub_fo = 0; double muq_fo =0; double mus_fo = 0;
            double e_fo = eos->temp_2_eps(IDB->Tfreeze,0,0,0);
                   T_fo = IDB->Tfreeze ;
            double p_fo = eos->pressure(e_fo,0,0,0);
	      
	    
	    // u^{\mu}	     
	    double utau_fo = 1.0 / sqrt ( 1.0 - _vx * _vx - _vy * _vy - _vz * _vz) ;

	    // if ( sqrt( _vx * _vx + _vy * _vy + _vz * _vz ) > 1.)
	    //   {
	    // 	_vx *= sqrt(0.999 / utau_fo);
	    // 	_vy *= sqrt(0.999 / utau_fo);
	    // 	_vz *= sqrt(0.999 / utau_fo);
	    // 	utau_fo = 0.999;
	    //   }   

	    

	    // viscous quantities	     
	    double pi_fo[10] = {0.0} ; 
	    double Pi_fo = 0.0 ;
	    for (int jt = 0; jt < 2; jt++)	    
	     for (int jx = 0; jx < 2; jx++)
	      for (int jy = 0; jy < 2; jy++)
		{
		  for (int i = 0; i < 10; i++)
		    pi_fo[i] += picube[jt][jx][jy][i] * wt[jt] * wx[jx] * wy[jy] ;
		  
		  Pi_fo += Picube[jt][jx][jy] * wt[jt] * wx[jx] * wy[jy] ;
		}
	    
	        
	    
	    freeze_file << tau_fo << " " << x_fo <<  " " << y_fo << " " << eta_fo << " " ;
 
	    freeze_file << norm_tau << " " << norm_x 
			<< " " << norm_y << " " << norm_eta << " " ;  

	    freeze_file << utau_fo << " " << utau_fo*_vx << " " << utau_fo*_vy << " " << utau_fo*_vz << " " ; 

	    freeze_file << e_fo*5.067653 << " " << T_fo*5.067653  << " " << mub_fo*5.067653 
                        << " " << mus_fo*5.067653  << " " << muq_fo*5.067653  << " " << (e_fo + p_fo) / T_fo << " " ; 

            freeze_file <<        pi_fo[c->indexpi(0,0)] * 5.067653  << " ";
            freeze_file <<        pi_fo[c->indexpi(0,1)] * 5.067653  << " ";
            freeze_file <<        pi_fo[c->indexpi(0,2)] * 5.067653  << " ";
            freeze_file << tau_fo*pi_fo[c->indexpi(0,3)] * 5.067653  << " ";
            freeze_file <<        pi_fo[c->indexpi(1,1)] * 5.067653  << " ";
            freeze_file <<        pi_fo[c->indexpi(1,2)] * 5.067653  << " ";
            freeze_file << tau_fo*pi_fo[c->indexpi(1,3)] * 5.067653  << " ";
            freeze_file <<        pi_fo[c->indexpi(2,2)] * 5.067653  << " ";
            freeze_file << tau_fo*pi_fo[c->indexpi(2,3)] * 5.067653  << " ";
            if (IDB->zetas_flag == 0 ) { 
                      freeze_file << pow(tau_fo,2)*pi_fo[c->indexpi(3,3)] * 5.067653  << endl ;
            }
            else{
                      freeze_file << pow(tau_fo,2)*pi_fo[c->indexpi(3,3)] * 5.067653 << " ";
	              freeze_file << Pi_fo*5.067653 << endl ; 
            }                                          
 

	  } // element loop
	
      } // ix,iy loop
  
}






void surf::get_surface3p1d(double tau)
{
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

  double _e,_p, _nb, _nq, _ns, _vx, _vy, _vz;  
  double Qcube[2][2][2][2][7] = {0.0}; // tau, x, y, eta and 7-conserved quantities.
  double picube[2][2][2][2][10] = {0.0}; // tau, x, y , eta and 10-components.
  double Picube[2][2][2][2] = {0.0}; //  tau,x , y and eta .

  cell* c;
  
  for (int ix = 2; ix < f->get_nx() - 2; ix += IDB->skip_fo_x )
   for (int iy = 2; iy < f->get_ny() - 2; iy += IDB->skip_fo_y )
    for (int ieta = 2; ieta < f->get_neta() - 2; ieta += IDB->skip_fo_eta )
      { // ix,iy,ieta loop

	// put the values in 16 corner points of the cube

        // (1)
	c = f->get_cell(ix, iy , ieta);
        c->get_Q(Qcube[1][0][0][0]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][0][0][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][0][0][0][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][0][0][0] = c->get_Pi();

        // (2)
	c = f->get_cell(ix , iy , ieta +  IDB->skip_fo_eta);
        c->get_Q(Qcube[1][0][0][1]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][0][0][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][0][0][1][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][0][0][1] = c->get_Pi();

        // (3)
	c = f->get_cell(ix , iy+  IDB->skip_fo_y , ieta );
        c->get_Q(Qcube[1][0][1][0]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][0][1][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][0][1][0][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][0][1][0] = c->get_Pi();

        // (4)
	c = f->get_cell(ix , iy+  IDB->skip_fo_y , ieta +  IDB->skip_fo_eta );
        c->get_Q(Qcube[1][0][1][1]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][0][1][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][0][1][1][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][0][1][1] = c->get_Pi();

        // (5)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy , ieta );
        c->get_Q(Qcube[1][1][0][0]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][1][0][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][1][0][0][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][1][0][0] = c->get_Pi();

        // (6)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy , ieta +  IDB->skip_fo_eta);
        c->get_Q(Qcube[1][1][0][1]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][1][0][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][1][0][1][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][1][0][1] = c->get_Pi();

        // (7)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy +  IDB->skip_fo_y , ieta );
        c->get_Q(Qcube[1][1][1][0]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][1][1][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][1][1][0][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][1][1][0] = c->get_Pi();

        // (8)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy +  IDB->skip_fo_y , ieta +  IDB->skip_fo_eta );
        c->get_Q(Qcube[1][1][1][1]);
	c->get_physical_var(eos, tau, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[1][1][1][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[1][1][1][1][c->indexpi(i,j)] = c->get_pi(i, j);}
	Picube[1][1][1][1] = c->get_Pi();

        // (9)  (at previous time ...)
	c = f->get_cell(ix, iy , ieta);
        c->get_Q_fo_prev(Qcube[0][0][0][0]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau() , _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][0][0][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][0][0][0][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][0][0][0] = c->get_Pi_fo_prev();

        // (10)
	c = f->get_cell(ix , iy , ieta +  IDB->skip_fo_eta);
        c->get_Q_fo_prev(Qcube[0][0][0][1]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][0][0][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][0][0][1][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][0][0][1] = c->get_Pi_fo_prev();

        // (11)
	c = f->get_cell(ix , iy+  IDB->skip_fo_y , ieta );
        c->get_Q_fo_prev(Qcube[0][0][1][0]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][0][1][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][0][1][0][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][0][1][0] = c->get_Pi_fo_prev();

        // (12)
	c = f->get_cell(ix , iy+  IDB->skip_fo_y , ieta +  IDB->skip_fo_eta );
        c->get_Q_fo_prev(Qcube[0][0][1][1]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][0][1][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][0][1][1][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][0][1][1] = c->get_Pi_fo_prev();

        // (13)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy , ieta );
        c->get_Q_fo_prev(Qcube[0][1][0][0]);
	c->get_center_var_fo_prev(eos,  tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][1][0][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][1][0][0][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][1][0][0] = c->get_Pi_fo_prev();

        // (14)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy , ieta +  IDB->skip_fo_eta);
        c->get_Q_fo_prev(Qcube[0][1][0][1]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][1][0][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][1][0][1][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][1][0][1] = c->get_Pi_fo_prev();

        // (15)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy +  IDB->skip_fo_y , ieta );
        c->get_Q_fo_prev(Qcube[0][1][1][0]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][1][1][0] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][1][1][0][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][1][1][0] = c->get_Pi_fo_prev();

        // (16)
	c = f->get_cell(ix +  IDB->skip_fo_x , iy +  IDB->skip_fo_y , ieta +  IDB->skip_fo_eta );
        c->get_Q_fo_prev(Qcube[0][1][1][1]);
	c->get_center_var_fo_prev(eos, tau - IDB->skip_fo_tau*f->get_dtau(), _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
	ccube[0][1][1][1] = _e;
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j <= i; j++){
	    picube[0][1][1][1][c->indexpi(i,j)] = c->get_pi_fo_prev(i, j);}
	Picube[0][1][1][1] = c->get_Pi_fo_prev();


	// cornelius magic	 
	cornelius->find_surface_4d(ccube);
	
	for (int isegm = 0; isegm < cornelius->get_Nelements(); isegm++)
	  { // element loop
	    
	    if (ix <= 4 || ix >= f->get_nx() - 4
		|| iy <= 4 || iy >= f->get_ny() - 4 || ieta <= 3 || ieta >= f->get_neta()-3)
	      {
		cout << "[Error] Freeze-out cell at the boundary !!! "
		     << "The grid is too small !!! exit(1)" << endl;
		exit(1);
	      }
	    
	    // fo position
	    double tau_fo = tau - IDB->skip_fo_tau*f->get_dtau() + cornelius->get_centroid_elem(isegm, 0) ;
	    double x_fo = f->get_x(ix) + cornelius->get_centroid_elem(isegm, 1) ;
	    double y_fo = f->get_y(iy) + cornelius->get_centroid_elem(isegm, 2) ;
	    double eta_fo =  f->get_eta(ieta) + cornelius->get_centroid_elem(isegm, 3) ;

	    
	    // normal vector
	    double norm_tau = cornelius->get_normal_elem(isegm, 0);
	    double norm_x = cornelius->get_normal_elem(isegm, 1);
	    double norm_y = cornelius->get_normal_elem(isegm, 2);
	    double norm_eta = cornelius->get_normal_elem(isegm, 3); 

	    // Linear interpolation below
	    double wt[2] = { 1. - cornelius->get_centroid_elem(isegm, 0) / IDB->skip_fo_tau*f->get_dtau(),
			     cornelius->get_centroid_elem(isegm, 0) / IDB->skip_fo_tau*f->get_dtau() } ;
	    double wx[2] = { 1. - cornelius->get_centroid_elem(isegm, 1) / (f->get_dx()* IDB->skip_fo_x),
			     cornelius->get_centroid_elem(isegm, 1) / (f->get_dx()* IDB->skip_fo_x)} ;
	    double wy[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / (f->get_dy()* IDB->skip_fo_y),
			    cornelius->get_centroid_elem(isegm, 2) / (f->get_dy()* IDB->skip_fo_y) } ;
	    double weta[2] = {1. - cornelius->get_centroid_elem(isegm, 3) / (f->get_deta()* IDB->skip_fo_eta),
			    cornelius->get_centroid_elem(isegm, 3) / (f->get_deta()* IDB->skip_fo_eta) } ;
	    
	    double Qc[7] = {0.0} ; // value at C.M. of the element
	    
	    for (int jt = 0; jt < 2; jt++)
	      for (int jx = 0; jx < 2; jx++)
		for (int jy = 0; jy < 2; jy++)
		 for (int jeta = 0; jeta < 2; jeta++)
		  for (int i = 0; i < 7; i++)
		    {
		      Qc[i] += Qcube[jt][jx][jy][jeta][i] * wt[jt] * wx[jx] * wy[jy] * weta[jeta];
		    }
	    
	    // e_fo, T_fo, nb_fo, nq_fo, ns_fo
	    for (int i = 0; i < 7; i++) { Qc[i] = Qc[i] / tau_fo ; } 
	    CN->CALC_2_LRF(eos, Qc, _e, _p, _nb, _nq, _ns , _vx, _vy, _vz);
	    double T_fo = eos->temperature(_e,_nb,_nq,_ns);
	    
	    
	    if ( T_fo < IDB->Tfreeze - 0.010  || T_fo > IDB->Tfreeze + 0.010 )
	      {
                cout << "[Error] Wrong interpolation in hypersurface"
		     << " !!! T_fo = " << T_fo <<", e_fo = " << _e
                     << ", nb_fo  = " << _nb << endl;
                cout << "Actual Tfreeze = " << IDB->Tfreeze << endl;
		exit(1);
	      }
	    	    
	    
            double mub_fo = 0; double muq_fo =0; double mus_fo = 0;
            double e_fo = eos->temp_2_eps(IDB->Tfreeze,0,0,0);
                   T_fo = IDB->Tfreeze ;
            double p_fo = eos->pressure(e_fo,0,0,0);
	      

	    // u^{\mu}	     
	    double utau_fo = 1.0 / sqrt ( 1.0 - _vx * _vx - _vy * _vy - _vz * _vz) ;

	    // if ( sqrt( _vx * _vx + _vy * _vy + _vz * _vz ) > 1.)
	    //   {
	    // 	_vx *= sqrt(0.999 / utau_fo);
	    // 	_vy *= sqrt(0.999 / utau_fo);
	    // 	_vz *= sqrt(0.999 / utau_fo);
	    // 	utau_fo = 0.999;
	    //   }   

	    
	    
	    // viscous quantities	     
	    double pi_fo[10] = {0.0} ; 
	    double Pi_fo = 0.0 ;
	    for (int jt = 0; jt < 2; jt++)	    
	     for (int jx = 0; jx < 2; jx++)
	      for (int jy = 0; jy < 2; jy++)
		for (int jeta = 0; jeta < 2; jeta++)
		{
		  for (int i = 0; i < 10; i++)
		    pi_fo[i] += picube[jt][jx][jy][jeta][i] * wt[jt] * wx[jx] * wy[jy] * weta[jeta];
		  
		  Pi_fo += Picube[jt][jx][jy][jeta] * wt[jt] * wx[jx] * wy[jy] * weta[jeta];
		}
	    
	    
	    freeze_file << tau_fo << " " << x_fo << " " << y_fo << " " << eta_fo << " " ;
 
	    freeze_file << norm_tau << " " << norm_x 
			<< " " << norm_y << " " << norm_eta << " " ;  

	    freeze_file << utau_fo << " " << utau_fo*_vx << " " << utau_fo*_vy << " " << utau_fo*_vz << " " ; 

	    freeze_file << e_fo*5.067653 << " " << T_fo*5.067653 << " " << mub_fo*5.067653
                        << " " << mus_fo*5.067653 << " " << muq_fo*5.067653 << " " << (e_fo + p_fo) / T_fo << " " ; 

            freeze_file <<        pi_fo[c->indexpi(0,0)] * 5.067653 << " ";
            freeze_file <<        pi_fo[c->indexpi(0,1)] * 5.067653 << " ";
            freeze_file <<        pi_fo[c->indexpi(0,2)] * 5.067653 << " ";
            freeze_file << tau_fo*pi_fo[c->indexpi(0,3)] * 5.067653 << " ";
            freeze_file <<        pi_fo[c->indexpi(1,1)] * 5.067653 << " ";
            freeze_file <<        pi_fo[c->indexpi(1,2)] * 5.067653 << " ";
            freeze_file << tau_fo*pi_fo[c->indexpi(1,3)] * 5.067653 << " ";
            freeze_file <<        pi_fo[c->indexpi(2,2)] * 5.067653 << " ";
            freeze_file << tau_fo*pi_fo[c->indexpi(2,3)] * 5.067653 << " ";
            if (IDB->zetas_flag == 0 ) { 
                      freeze_file << pow(tau_fo,2)*pi_fo[c->indexpi(3,3)] * 5.067653  << endl ;
            }
            else{
                      freeze_file << pow(tau_fo,2)*pi_fo[c->indexpi(3,3)] * 5.067653 << " ";
	              freeze_file << Pi_fo*5.067653 << endl ; 
            }                                     
                   	    
	  } // element loop
	


      } // ix,iy,ieta loop
	
  
}



























