#include "master.h"

master::master(idb* _IDB)
{
  // Creating the "hydro_output" directory if it doesn't exist
  if (mkdir("hydro_output", 0777) == -1)
    cerr << "[Info] " << strerror(errno) << endl;
  else
    cout << "[Info] hydro_output directory is created" << endl;

  IDB = _IDB;

  if(IDB->eos == 0)
   {
     eos = new EoS0();
   }
  if(IDB->eos == 1)
    {
     eos = new EoS1();
    }
  else
    {
      eos = new EoS2();
    }
  
  CN = new cnvrt();
  
}

master::~master()
{
  delete IDB;
  delete eos;
  delete out;
  delete h;
  delete g;
  delete CN;
}

void master::init()
{
  trcoef = new trancoeff( eos, IDB);

  // make the grid
  g = new grid(IDB,CN, trcoef,eos);
  g->make_grid();


  // set IC on the grid
  if(IDB->ic_mode == 0)
    {
      opt_glau* ic = new opt_glau(IDB);
      ic->set_ic(g,eos);
    }

  if(IDB->ic_mode == 1)
    {
      mc_glau* ic = new mc_glau(IDB);
      ic->set_ic(g,eos);
    }

  if(IDB->ic_mode == 2)
    {
      read_ic* ic = new read_ic(IDB);
      ic->set_ic(g,eos);
    }

 
 // wanna the hydro field output ?
 out = new fluid_info(g,eos,IDB);


}

void master::run_hydro()
{

  int stop_flag ; 
 
    // freezeout hypersurface initialization for therminator ...
    evolve* map =new evolve();

  if( IDB->eps_freeze_flag == 0 ){
    map->ini(g,IDB->tau0,IDB->tauMax,IDB->dtau);
    map->put(0,g,IDB->tau0,eos);
  }


  sf = new surf(eos, IDB , g , CN );
  h = new hydro(eos, g , IDB , IDB->tau0 , IDB->dtau, CN, trcoef);


  int nstep = (IDB->tauMax - IDB->tau0)/IDB->dtau;
  auto start = high_resolution_clock::now();  
  for(int istep=0; istep<nstep ; istep++)
    { // step loop


      // save the fluis info on grid
      if( istep*1.0/IDB->save_every_N_steps == int (istep/IDB->save_every_N_steps) )
        {
          out->midslice_xy_dist(h->get_tau());
          out->conservation_check(h->get_tau());
          out->anisotropy_in_xy_plane(h->get_tau());
        }


      // Befor the 1st step save for freezeout.
      if( istep == 0 ) { 
        save_for_fo(istep);
        stop_flag = check_to_stop( g, eos, h->get_tau() ) ;
      }

      // evolution step.
      auto stop = high_resolution_clock::now();
      auto sdxduration = duration_cast<seconds>(stop - start);
      cout <<"[ "<< sdxduration.count()/60 << " mins. ] ";
      printf( "\t\033[93m[ tau : %2.2f fm ] \033[0m",h->get_tau());
      cout << "\t" << "[ "<<istep<<"/"<<nstep<<" ] completed." << endl;
     


      // entire hydro evolution here
      h->evolve();

      // find the hypersurface after the evolution and save for freezeout.
      find_hyper_surface(istep+1, h->get_tau()); 
      save_for_fo(istep+1); 


      if( IDB->eps_freeze_flag == 0 ){
        int gg=(istep+1)/map->stepsave;
         if ((gg*map->stepsave)==(istep+1))
	  {
	    map->put(gg,g,h->get_tau(),eos);
	    cout<<"\033[32m [Info] Saving for time="<<h->get_tau()<<"\tEvolution step="<<istep+1<<
	      ", "<<gg+1<<"step saved. \033[0m "<<endl;
	  }
      }


     // check wheather all cells are below eps freeze or not ?
      stop_flag = check_to_stop( g, eos, h->get_tau() ) ;
          if (stop_flag > 0) 
            { continue; }
	  else 
            { cout<<"[complete] all cells freezed..."<<endl; break ;}
	     
      
    }// step loop


  if( IDB->eps_freeze_flag == 0 ){
    // Getting the freezeout hypersurface in .xml file .
    cout << "$ $ $ HYPERSURFACE FOR THERMINATOR $ $ $" << endl ; 
    double Tfreeze = IDB->Tfreeze  ;
    cout<<"\n\n freezeout temperature : "<<Tfreeze<<" GeV"<<endl;
    map->gethypersurface(g,eos,Tfreeze);
    map->hypersurface("hydro_output/",1);
  }
  

}



int master::check_to_stop(grid* f, EoS* eos, double tau )
{
  double e,p,nb,nq,ns,vx,vy,vz;
  int flag = 0;
  double max_temp = 0;
  double max_eps = 0;
  double max_nb = 0;

  for(int ix = 0; ix < IDB->nx ; ix++)
    for(int iy = 0; iy <IDB->ny; iy++)
      for(int iz = 0; iz <IDB->neta ; iz++)
       {
	 cell* c = f->get_cell(ix,iy,iz); 
	 c->get_physical_var(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
	 double temperature = eos->temperature(e,nb,nq,ns);

         if(ix < 4 && temperature > IDB->Tfreeze + 0.01 ) {cout << "SMALL GRID [EXIT]" << endl; exit(1) ; }
         if(ix > IDB->nx-4 && temperature > IDB->Tfreeze + 0.01 ) {cout << "SMALL GRID [EXIT]" << endl; exit(1) ; }
         if(iy < 4 && temperature > IDB->Tfreeze + 0.01 ) {cout << "SMALL GRID [EXIT]" << endl; exit(1) ; }
         if(iy > IDB->ny-4 && temperature > IDB->Tfreeze + 0.01 ) {cout << "SMALL GRID [EXIT]" << endl; exit(1) ; }
         if(IDB->neta != 1 && iz < 4 && temperature > IDB->Tfreeze + 0.01 ) {cout << "SMALL GRID [EXIT]" << endl; exit(1) ; }
         if(IDB->neta != 1 && iz > IDB->neta-4 && temperature > IDB->Tfreeze + 0.01 ) {cout << "SMALL GRID [EXIT]" << endl; exit(1) ; }

         if (temperature > max_temp ) { max_temp = temperature ; }
         if (e > max_eps) { max_eps = e ; }
         if (nb > max_nb) { max_nb = nb ; }

         if(IDB->eps_freeze_flag == 0 && temperature > IDB->Tfreeze - 0.001 ) { flag += 1; }
         if(IDB->eps_freeze_flag == 1 && e > IDB->eps_freeze - 0.001 ) { flag += 1; }
         
      }

       printf( "Max e : %2.3f GeV/fm^3   Max T : %3.3f MeV.  ", max_eps, max_temp*1000.0);
  
  return flag;
}



void master::save_for_fo(int istep)
{
// save for freezeout in every [skip_fo_tau] steps
if( fabs( (istep*1.0) / IDB->skip_fo_tau - 
        int ( (istep*1.0) / IDB->skip_fo_tau ) ) < 0.0001  )
{
  cell *c;
   for(int ix = 0; ix<g->get_nx(); ix++)
    for(int iy = 0; iy<g->get_ny(); iy++)
     for(int iz = 0; iz<g->get_neta(); iz++)
	{      
	  c = g->get_cell(ix, iy, iz );
	  c->save_for_fo();

	}

}
else
{
return ;
}
  
}

void master::find_hyper_surface(int istep, double tau)
{
// finding the freezeout hypersurface in every [skip_fo_tau] steps
if( fabs ( (istep*1.0) / IDB->skip_fo_tau - 
        int ( (istep*1.0) / IDB->skip_fo_tau ) ) < 0.0001 && istep >= 1)
{
  if(g->get_neta() == 1 )
    { 
     sf->get_surface2p1d(tau);
    }
  else
      {
     sf->get_surface3p1d(tau);
      }
}



}



