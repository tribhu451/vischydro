#include "hydro.h"

#define del_eps 1e-7

hydro::hydro(EoS* _eos, grid* _f,idb* _IDB,double _tau0, double _dt, cnvrt* _CN, trancoeff* _tr)
{
 eos = _eos;
 f = _f;
 CN = _CN;
 tau = _tau0;
 IDB = _IDB ;
 dt =_dt;
 hr = new hrr();
 trcoef = _tr;

}

hydro::~hydro()
{
  delete trcoef;
  delete hr;
}


void hydro::set_dtau(double deltaTau)
{
  dt = deltaTau;
  if (dt > IDB->dx / 2. ||
      dt > IDB->dy / 2. /*|| dt > tau*IDB->deta */) {
    cout << "too big delta_tau " << dt << "  " << IDB->dx << "  " << IDB->dy
	 << "  " << tau * IDB->deta << endl;
    exit(1);
  }
}


void hydro::evolve()
{

  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix<f->get_nx(); ix++)
	{      
	  cell *c = f->get_cell(ix, iy, iz );
	  c->save_Q_prev();
	  c->clear_flux();
	}

  // X-direction flux
  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix<f->get_nx()-1; ix++)
	{
	  hlle_flux(f->get_cell(ix,iy, iz), f->get_cell(ix+1,iy,iz), X_, PREDICT);
	  // cout<<iy<<"\t"<<iz<<"\t"<<ix<<endl;
	}
  
  // Y-direction flux
  for(int iz = 0; iz<f->get_neta(); iz++)
    for(int ix = 0; ix<f->get_nx(); ix++)
      for(int iy = 0; iy<f->get_ny()-1; iy++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy+1,iz), Y_, PREDICT);
	}
  
if(f->get_neta() > 1) // don,t calculate z-flux in 2+1D hydro
{  
  // Z-direction flux
  for(int ix = 0; ix<f->get_nx(); ix++)
    for(int iy = 0; iy<f->get_ny(); iy++)
      for(int iz = 0; iz<f->get_neta() -1; iz++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy,iz+1), Z_, PREDICT);
	}
}  
  
  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix< f->get_nx(); ix++)
        {
	  cell *c = f->get_cell(ix, iy,iz);
	  sourcestep( PREDICT,  ix,  iy, iz,  tau);
	  c->update_Q_to_Qh_by_flux();
	  c->clear_flux();
	}
  
  // X-direction flux
  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix<f->get_nx()-1; ix++)
	{
	  hlle_flux(f->get_cell(ix,iy, iz), f->get_cell(ix+1,iy,iz), X_, CORRECT);
	}
  
  
  // Y-direction flux
  for(int iz = 0; iz<f->get_neta(); iz++)
    for(int ix = 0; ix<f->get_nx(); ix++)
      for(int iy = 0; iy<f->get_ny()-1; iy++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy+1,iz), Y_, CORRECT);
	}
  
if(f->get_neta()> 1) // don,t calculate z-flux in 2+1D hydro
{  
  // Z-direction flux
  for(int ix = 0; ix<f->get_nx(); ix++)
    for(int iy = 0; iy<f->get_ny(); iy++)
      for(int iz = 0; iz<f->get_neta()-1; iz++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy,iz+1), Z_, CORRECT);
	}
  
}  
  for(int iy = 0; iy <f->get_ny(); iy++)
    for(int iz = 0; iz< f->get_neta(); iz++)
      for(int ix = 0; ix< f->get_nx(); ix++)
	{
	  cell *c = f->get_cell(ix, iy,iz); 
	  sourcestep( CORRECT,  ix,  iy, iz, tau);
	  c->update_by_flux();
	  c->clear_flux();	  
	}
  
  tau = tau+dt;                 // tau increased
  f->correct_imaginary_cells(); //boundary condition
  
  
  if (trcoef->isViscous())
    {
      //cout<<"viscous calculation ... "<<endl;
      ISformal();  // evolution of viscous quantities according to IS equations
      
      // X dir
      for (int iy = 0; iy < f->get_ny(); iy++)
	for (int iz = 0; iz <f->get_neta(); iz++)
	  for (int ix = 0; ix <f->get_nx()- 1; ix++)
	    {
	      visc_flux(f->get_cell(ix, iy, iz), f->get_cell(ix + 1, iy, iz), X_);
	    }
      
      // Y dir
      for (int iz = 0; iz < f->get_neta(); iz++)
	for (int ix = 0; ix < f->get_nx(); ix++)
	  for (int iy = 0; iy <f->get_ny() - 1; iy++)
	    {
	      visc_flux(f->get_cell(ix, iy, iz), f->get_cell(ix, iy + 1, iz), Y_);
	    }
if(f->get_neta() > 1) // don,t calculate z-flux in 2+1D hydro
{ 
      // Z dir
      for (int ix = 0; ix < f->get_nx(); ix++)
	for (int iy = 0; iy <f->get_ny(); iy++)
	  for (int iz = 0; iz < f->get_neta() - 1; iz++)
	    {
	      visc_flux(f->get_cell(ix, iy, iz), f->get_cell(ix, iy, iz + 1), Z_);
	    }
      
}
      
      for (int iy = 0; iy < f->get_ny(); iy++)
	for (int iz = 0; iz < f->get_neta(); iz++)
	  for (int ix = 0; ix < f->get_nx(); ix++)
	    {
	      visc_source_step(ix, iy, iz);
	      f->get_cell(ix, iy, iz)->update_by_visc_flux();
	      f->get_cell(ix, iy, iz)->clear_flux();
	    }
      
    }
  else
    {  // end viscous part
    }
 f->correct_imaginary_cells_full();
  
}




void hydro::hlle_flux(cell* left, cell* right, int direction, int mode)
{
  const double dta = mode == 0 ? dt/2. : dt;
  double el,vxl,vyl,vzl,nbl,nql,nsl,pl;
  double er,vxr,vyr,vzr,nbr,nqr,nsr,pr;
  double El,Er;
  double Utl=0.,Uxl=0.,Uyl=0.,Uzl=0.,Ubl=0.,Uql=0.,Usl=0.;
  double Utr=0.,Uxr=0.,Uyr=0.,Uzr=0.,Ubr=0.,Uqr=0.,Usr=0.;
  double Ftl=0.,Fxl=0.,Fyl=0.,Fzl =0,Fbl=0., Fql =0, Fsl =0;
  double Ftr=0.,Fxr=0.,Fyr=0.,Fzr =0,Fbr=0., Fqr =0, Fsr =0;
  double csb,vb,bl=0.,br=0.;
  double flux[7]={0.0};
  double tauFactor;  // fluxes are also multiplied by tau
  double dx =0;
  
  if(mode == PREDICT)
    {
      left->get_right_var(eos, tau, el, pl, nbl,nql,nsl, vxl, vyl, vzl, direction);
      right->get_left_var(eos, tau, er, pr, nbr,nqr,nsr, vxr, vyr, vzr, direction);
      El = (el+pl)*(1.0/(1.0-vxl*vxl-vyl*vyl-vzl*vzl));
      Er = (er+pr)*(1.0/(1.0-vxr*vxr-vyr*vyr-vzr*vzr));
      tauFactor = tau ;
    }
  else
    {
      left->get_right_varH(eos, tau, el, pl, nbl,nql,nsl, vxl, vyl, vzl, direction );
      right->get_left_varH( eos, tau, er, pr, nbr,nqr,nsr, vxr, vyr, vzr, direction);
      El = (el+pl)*(1.0/(1.0-vxl*vxl-vyl*vyl-vzl*vzl));
      Er = (er+pr)*(1.0/(1.0-vxr*vxr-vyr*vyr-vzr*vzr));
      tauFactor = tau +(0.5*dt);
    }
  
  if (el < 0.)
    {
      el = 0.;
      pl = 0.;
    }
  if (er < 0.)
    {
      er = 0.;
      pr = 0.;
    }
  
  if (el < del_eps && er < del_eps ) return;  // *1) no flux calculation if both sides are empty cells
  
  double gammal = 1.0 / sqrt(1 - vxl * vxl - vyl * vyl - vzl*vzl );
  double gammar = 1.0 / sqrt(1 - vxr * vxr - vyr * vyr - vzr*vzr );
  Utl = gammal * gammal * (el + pl) - pl;
  Uxl = gammal * gammal * (el + pl) * vxl;
  Uyl = gammal * gammal * (el + pl) * vyl;
  Uzl = gammal * gammal * (el + pl) * vzl;
  Ubl = gammal * nbl;
  Uql = gammal * nql;
  Usl = gammal * nsl;
  
  
  Utr = gammar * gammar * (er + pr) - pr;
  Uxr = gammar * gammar * (er + pr) * vxr;
  Uyr = gammar * gammar * (er + pr) * vyr;
  Uzr = gammar * gammar * (er + pr) * vzr;
  Ubr = gammar * nbr;
  Uqr = gammar * nqr;
  Usr = gammar * nsr;
  
  
  if(direction == X_)
    {
      Ftl = Utl * vxl + pl * vxl;
      Fxl = Uxl * vxl + pl;
      Fyl = Uyl * vxl;
      Fzl = Uzl * vxl;
      Fbl = Ubl * vxl;
      Fql = Uql * vxl;
      Fsl = Usl * vxl;
      
      Ftr = Utr * vxr + pr * vxr;
      Fxr = Uxr * vxr + pr;
      Fyr = Uyr * vxr;
      Fzr = Uzr * vxr;
      Fbr = Ubr * vxr;
      Fqr = Uqr * vxr;
      Fsr = Usr * vxr;
      
      
      // for the case of constant c_s only
      csb = sqrt(eos->cs2_() +
		 0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vxl - vxr, 2));

      vb = (sqrt(El) * vxl + sqrt(Er) * vxr) / (sqrt(El) + sqrt(Er));
      bl = min(0., min((vb - csb) / (1 - vb * csb),
		       (vxl - eos->cs_()) / (1 - vxl * eos->cs_())));
      br = max(0., max((vb + csb) / (1 + vb * csb),
		       (vxr + eos->cs_()) / (1 + vxr * eos->cs_())));
      
      dx = f->get_dx();
      
      
      if (el == 0.) bl = -1.; 
      if (er == 0.) br = 1.;
    }
  
  
  
  if(direction == Y_)
    {
      Ftl = Utl * vyl + pl * vyl;
      Fxl = Uxl * vyl ;
      Fyl = Uyl * vyl + pl;
      Fzl = Uzl * vyl;
      Fbl = Ubl * vyl;
      Fql = Uql * vyl;
      Fsl = Usl * vyl;
      
      Ftr = Utr * vyr + pr * vyr;
      Fxr = Uxr * vyr ;
      Fyr = Uyr * vyr + pr;
      Fzr = Uzr * vyr;
      Fbr = Ubr * vyr;
      Fqr = Uqr * vyr;
      Fsr = Usr * vyr;
      
      


      // for the case of constant c_s only
      
     csb = sqrt(eos->cs2_() +
		 0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vyl - vyr, 2));
     vb = (sqrt(El) * vyl + sqrt(Er) * vyr) / (sqrt(El) + sqrt(Er));
     bl = min(0., min((vb - csb) / (1 - vb * csb),
		      (vyl - eos->cs_()) / (1 - vyl * eos->cs_())));
     br = max(0., max((vb + csb) / (1 + vb * csb),
                   (vyr + eos->cs_()) / (1 + vyr * eos->cs_())));

   
     dx = f->get_dy();
     
     
     if (el == 0.) bl = -1.;
     if (er == 0.) br = 1.;
    }
  
  
  if(direction == Z_)
    {
      double tau1 =tauFactor;
      Ftl = Utl * vzl/tau1 + pl * vzl/tau1;
      Fxl = Uxl * vzl/tau1 ;
      Fyl = Uyl * vzl/tau1;
      Fzl = Uzl * vzl/tau1 + pl/tau1;
      Fbl = Ubl * vzl/tau1;
      Fql = Uql * vzl/tau1;
      Fsl = Usl * vzl/tau1;
      
      Ftr = Utr * vzr/tau1 + pr * vzr/tau1;
      Fxr = Uxr * vzr/tau1 ;
      Fyr = Uyr * vzr/tau1 ;
      Fzr = Uzr * vzr/tau1 + pr/tau1 ;
      Fbr = Ubr * vzr/tau1;
      Fqr = Uqr * vzr/tau1;
      Fsr = Usr * vzr/tau1;
      
      
      
    // for the case of constant c_s only
     
     csb = sqrt(eos->cs2_() +
		 0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vzl - vzr, 2));

     vb = (sqrt(El) * vzl + sqrt(Er) * vzr) / (sqrt(El) + sqrt(Er));
     bl = (1.0/tau)*min(0., min((vb - csb) / (1 - vb * csb),
				(vzl - eos->cs_()) / (1 - vzl * eos->cs_())));
     br = (1.0/tau)*max(0., max((vb + csb) / (1 + vb * csb),
				 (vzr + eos->cs_()) / (1 + vzr * eos->cs_())));

    


     dx =f->get_deta();
     
     
    if (el == 0.) bl = -1./tau; 
    if (er == 0.) br = 1./tau;
    }
  

  if(bl == 0. && br == 0.) return;
  
  flux[T_] = tauFactor * dta / dx *
    (-bl * br * (Utl - Utr) + br * Ftl - bl * Ftr) / (-bl + br);
  flux[X_] = tauFactor * dta / dx *
    (-bl * br * (Uxl - Uxr) + br * Fxl - bl * Fxr) / (-bl + br);
  flux[Y_] = tauFactor * dta / dx *
    (-bl * br * (Uyl - Uyr) + br * Fyl - bl * Fyr) / (-bl + br);
  flux[Z_] = tauFactor * dta / dx *
    (-bl * br * (Uzl - Uzr) + br * Fzl - bl * Fzr) / (-bl + br);
  flux[NB_] = tauFactor * dta / dx *
    (-bl * br * (Ubl - Ubr) + br * Fbl - bl * Fbr) / (-bl + br);
  flux[NQ_] = tauFactor * dta / dx *
    (-bl * br * (Uql - Uqr) + br * Fql - bl * Fqr) / (-bl + br);
  flux[NS_] = tauFactor * dta / dx *
    (-bl * br * (Usl - Usr) + br * Fsl - bl * Fsr) / (-bl + br);
  
  
  left->add_flux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
		-flux[NQ_], -flux[NS_]);
  right->add_flux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
		 flux[NS_]);

}




void hydro::sourcestep(int mode, int ix, int iy, int iz, double _tau)
{ 
  double S[7];
  double e,p,vx,vy,vz,nb,nq,ns;
  if( mode == PREDICT )
    {
      double _dt  = 0.5*dt; // (1/2)*dtau -> for predictor step 
      f->get_cell(ix,iy,iz)->get_Q(S);
      for(int i =0; i<7; i++){S[i] = S[i]/_tau; }
      CN->CALC_2_LRF(eos,  S, e, p, nb, nq, ns, vx, vy, vz);
      if (e < del_eps ) return;
      f->get_cell(ix,iy,iz)->add_flux( (-S[T_] * vz * vz - p * (1. + vz * vz))*_dt, 0.0, 0.0, -S[Z_]*_dt, 0.0,0.0, 0.0);  
    }
  else
    {
      _tau = _tau+0.5*dt;
      double _dt  = dt;
      f->get_cell(ix,iy,iz)->get_Qh(S);
      for(int i =0; i<7; i++){S[i] = S[i]/_tau; }
      CN->CALC_2_LRF(eos,  S, e, p, nb, nq, ns, vx, vy, vz);
      if (e < del_eps) return;
      f->get_cell(ix,iy,iz)->add_flux(  (-S[T_] * vz * vz - p * (1. + vz * vz))*_dt, 0.0, 0.0, (-S[Z_])*_dt, 0.0,0.0, 0.0);     
    }
  
}





         ////////////////////
        ///// Viscous //////
       ////////////////////





void hydro::NSquant(int ix,int iy,int iz, double dmunu[4][4], double &PiNS, double piNS[4][4] )
{



  
  double e1,p1,nb1,nq1,ns1,vx1,vy1,vz1,e0,p0,vx0,vy0,vz0,nb0,nq0,ns0;
  double ut0=1.,ux0=0.,uy0=0.,uz0=0.,ut1=1.,ux1=0.,uy1=0.,uz1=0.; //u^{mu}
  dmunu[4][4]={0.0}; //\partial_{mu} u^{\nu}
  double Delta[4][4] = {0.0};
  //g^{\mu \nu} 
  double gmunu[4] = { 1.0, -1.0 , -1.0, -1.0/pow(tau-0.5*dt,2) }; 

  // current cells variable at half time step
  double eh,ph,nbh,nqh,nsh,vxh,vyh,vzh ;   
  f->get_cell(ix,iy,iz)->get_center_varH(eos, tau-0.5*dt, eh, ph, nbh, nqh, nsh, vxh, vyh, vzh);  

  // current cells variable at full time step
  double e_loc,p_loc,nb_loc,nq_loc,ns_loc,vx_loc,vy_loc,vz_loc;
  f->get_cell(ix,iy,iz)->get_physical_var(eos, tau, e_loc, p_loc, nb_loc, nq_loc, ns_loc, vx_loc, vy_loc, vz_loc);
  

  // those cells which are situated at extreme end ( +/-x or +/-y side only)
  if(ix == 0 || ix == f->get_nx()-1 || iy == 0 || iy == f->get_ny()-1 )
    {
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  {
	    dmunu[i][j] = 0.0;
	    Delta[i][j] = 0.0;
	    piNS[i][j] = 0.0;
	  }
      PiNS=0.0;
      
      return ;
    }
  

  if(eh < del_eps ) 
{
      for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  {
	    dmunu[i][j] = 0.0;
	    Delta[i][j] = 0.0;
	    piNS[i][j] = 0.0;
	  }
      PiNS=0.0;
      
      return ;
    
}
  
  cell* c = f->get_cell(ix,iy,iz);
  
  // \partial_{;\tau} u^{\tau,x,y,\eta}
  c->get_center_var_prev(eos, tau-dt, e0, p0, nb0, nq0, ns0, vx0, vy0, vz0);
  if(e0 > del_eps && e1 > del_eps)
    {

      ut0 = 1.0 / sqrt( 1.0 - vx0*vx0 - vy0*vy0 - vz0*vz0 ) ;
      ux0 = ut0 * vx0 ;
      uy0 = ut0 * vy0 ;
      uz0 = ut0 * vz0 / (tau-dt) ;

      ut1 = 1.0 / sqrt( 1.0 - vx_loc*vx_loc - vy_loc*vy_loc - vz_loc*vz_loc );
      ux1 = ut1 * vx_loc ;
      uy1 = ut1 * vy_loc ;
      uz1 = ut1 * vz_loc / tau ;

      dmunu[0][0] = (ut1-ut0) / dt ;
      dmunu[0][1] = (ux1-ux0) / dt ;
      dmunu[0][2] = (uy1-uy0) / dt ;
      dmunu[0][3] = (uz1-uz0) / dt ;
    }
  else
    {
      dmunu[0][0] = 0. ;
      dmunu[0][1] = 0. ;
      dmunu[0][2] = 0. ;
      dmunu[0][3] = 0. ;
    }
  
  
  //\partial_{;x} u^{\tau,x,y,\eta}
  f->get_cell(ix-1,iy,iz)->get_center_varH(eos, tau-0.5*dt , e0, p0, nb0, nq0, ns0, vx0, vy0, vz0);
  f->get_cell(ix+1,iy,iz)->get_center_varH(eos, tau-0.5*dt , e1, p1, nb1, nq1, ns1, vx1, vy1, vz1);
  if(e0 > del_eps && e1 > del_eps )
    {
      ut0 = 1.0 / sqrt( 1.0 - vx0*vx0 - vy0*vy0 -vz0*vz0 );
      ux0 = ut0 * vx0 ;
      uy0 = ut0 * vy0 ;
      uz0 = ut0 * vz0 / (tau-0.5*dt) ; 
      ut1 = 1.0 / sqrt( 1.0 - vx1*vx1 - vy1*vy1 -vz1*vz1 );
      ux1 = ut1 * vx1 ;
      uy1 = ut1 * vy1 ;
      uz1 = ut1 * vz1 / (tau-0.5*dt) ;

      dmunu[1][0] = ( ut1 - ut0 ) / ( 2.0 * f->get_dx() ) ;
      dmunu[1][1] = ( ux1 - ux0 ) / ( 2.0 * f->get_dx() ) ;
      dmunu[1][2] = ( uy1 - uy0 ) / ( 2.0 * f->get_dx() ) ;
      dmunu[1][3] = ( uz1 - uz0 ) / ( 2.0 * f->get_dx() ) ;

    }
  else
    {
      dmunu[1][0] = 0. ;
      dmunu[1][1] = 0. ;
      dmunu[1][2] = 0. ;
      dmunu[1][3] = 0. ;
    }
  
  //\partial_{;y} u^{\tau,x,y,\eta}
  f->get_cell(ix,iy-1,iz)->get_center_varH(eos, tau-0.5*dt, e0, p0, nb0, nq0, ns0, vx0, vy0, vz0);
  f->get_cell(ix,iy+1,iz)->get_center_varH(eos, tau-0.5*dt, e1, p1, nb1, nq1, ns1, vx1, vy1, vz1);
  if(e0 > del_eps && e1 > del_eps)
    {
      ut0 = 1.0 / sqrt( 1.0 - vx0*vx0 - vy0*vy0 -vz0*vz0 ) ;
      ux0 = ut0 * vx0 ;
      uy0 = ut0 * vy0 ;
      uz0 = ut0 * vz0 / (tau - 0.5*dt) ; 
      ut1 = 1.0 / sqrt( 1.0 - vx1*vx1 - vy1*vy1 - vz1*vz1 );
      ux1 = ut1 * vx1 ;
      uy1 = ut1 * vy1 ;
      uz1 = ut1 * vz1 / (tau-0.5*dt) ;


      dmunu[2][0] = (ut1-ut0) / ( 2.0 * f->get_dy() ) ;
      dmunu[2][1] = (ux1-ux0) / ( 2.0 * f->get_dy() ) ;
      dmunu[2][2] = (uy1-uy0) / ( 2.0 * f->get_dy() ) ;
      dmunu[2][3] = (uz1-uz0) / ( 2.0 * f->get_dy() ) ;

    }
  else
    {
      dmunu[2][0] = 0.;
      dmunu[2][1] = 0. ;
      dmunu[2][2] = 0.;
      dmunu[2][3] = 0. ;
  }
  
  
if(f->get_neta() > 1)
{ 

      //\partial_{;\eta} u^{\tau,x,y,\eta}
      f->get_cell(ix,iy,iz-1)->get_center_varH(eos, tau-0.5*dt, e0, p0, nb0, nq0, ns0, vx0, vy0, vz0);
      f->get_cell(ix,iy,iz+1)->get_center_varH(eos, tau-0.5*dt, e1, p1, nb1, nq1, ns1, vx1, vy1, vz1);
      if(e0 > del_eps && e1 > del_eps)
	{
	  ut0 = 1.0 / sqrt( 1.0 - vx0*vx0 - vy0*vy0 - vz0*vz0 ) ;
	  ux0 = ut0 * vx0 ;
	  uy0 = ut0 * vy0 ;
	  uz0 = ut0 * vz0 / (tau-0.5*dt) ; 
	  ut1 = 1.0 / sqrt( 1.0 - vx1*vx1 - vy1*vy1 - vz1*vz1 ) ;
	  ux1 = ut1 * vx1 ;
	  uy1 = ut1 * vy1 ;
	  uz1 = ut1 * vz1 / (tau-0.5*dt) ;

	  dmunu[3][0] = (ut1-ut0) / ( 2.0 * f->get_deta() ) ;
	  dmunu[3][1] = (ux1-ux0) / ( 2.0 * f->get_deta() ) ;
	  dmunu[3][2] = (uy1-uy0) / ( 2.0 * f->get_deta() ) ;
	  dmunu[3][3] = (uz1-uz0) / ( 2.0 * f->get_deta() ) ;

	}
      else
	{
	  dmunu[3][0] = 0. ;
	  dmunu[3][1] = 0. ;
	  dmunu[3][2] = 0. ;
	  dmunu[3][3] = 0. ;
	}
      
}

else {
	  dmunu[3][0] = 0. ;
	  dmunu[3][1] = 0. ;
	  dmunu[3][2] = 0. ;
	  dmunu[3][3] = 0. ;
}
 

  // u^{\mu} to calculate \Delta^{\mu \nu}
  double u[4] = {1.0,0.0,0.0,0.0};

  //some extra from christoffel symbols
  //simultaneously u^{\mu} will be calculated
  if(eh > del_eps )
    {
      u[0] = 1.0 / sqrt( 1.0 - vxh*vxh - vyh*vyh - vzh*vzh ) ;
      u[1] = u[0] * vxh;
      u[2] = u[0] * vyh;
      u[3] = u[0] * vzh / (tau-0.5*dt) ; 
      dmunu[3][0] += (tau-0.5*dt) * u[3] ;
      dmunu[3][3] += ( 1.0 / (tau-0.5*dt) ) * u[0] ;
      dmunu[0][3] += ( 1.0 / (tau-0.5*dt) ) * u[3] ;
    }

 
  /////////////////////////////////////////////////////
  //    calculate \Pi_{NS} and \pi_{NS}^{\mu \nu}    //
  /////////////////////////////////////////////////////
  
  
  // \Delta^{\mu \nu}
  for(int mu=0; mu<4; mu++) 
    for(int nu=0; nu<4; nu++)
      {
        if (mu == nu) 
          Delta[mu][nu] = gmunu[mu] - u[mu] * u[nu];
        else
          Delta[mu][nu] = - u[mu] * u[nu];
      }
  

  // \partial_{;\lambda} u^{\lambda}
  double theta = dmunu[0][0] + dmunu[1][1] + dmunu[2][2] + dmunu[3][3] ;

  
  double etaS,zetaS;
  //get \eta & \zeta
  double temp = eos->temperature(e_loc,nb_loc,nq_loc,ns_loc);
  trcoef->getEta(e_loc, temp, etaS, zetaS);
  double ent = eos->entropy(e_loc,nb_loc,nq_loc,ns_loc) ;
  double eta = etaS * ent  ;
  double zeta = zetaS * ent ;
 

  // -\zeta \partial_{\lambda} u^{\lambda}     
  PiNS = - zeta * theta / 5.068; // fm^{-4} --> GeV/fm^3
  
  double z[4][4] = {0.0};
  for(int mu=0; mu<4; mu++) 
    for(int nu=0; nu<4; nu++)
      for(int lambda=0; lambda<4; lambda++)
	{
	  z[mu][nu] +=  0.5 * Delta[mu][lambda] * dmunu[lambda][nu] 
	    + 0.5 * Delta[nu][lambda] * dmunu[lambda][mu] ;	  
	}
  
  for(int mu=0; mu<4; mu++) 
    for(int nu=0; nu<4; nu++)
       piNS[mu][nu] = 2 * eta * ( z[mu][nu] - (1.0/3.0) * Delta[mu][nu] * theta ) / 5.068;
       // hbarc = 5.068, fm^{-4} --> GeV/fm^3 
  

 //--------- debug part: NaN/inf check, trace check, diag check, transversality
 // check
 for (int i = 0; i < 4; i++)
  for (int j = 0; j < 4; j++) {
   if (piNS[i][j] != 0. && fabs(1.0 - piNS[j][i] / piNS[i][j]) > 1e-5)
    cout << "non-diag: " << piNS[i][j] << "  " << piNS[j][i] << endl;
   if (std::isinf(piNS[i][j]) || std::isnan(piNS[i][j])) {
    cout << "hydro:NSquant: inf/nan i " << i << " j " << j << endl;
    exit(1);
   }
  }
  
  
}






// function to solve IS equation
void hydro::ISformal()
{
  
  double eps,prs,nb,ns,nq,vx,vy,vz;
  double u[4] = {1,0,0,0};
  double flux[7]={0.0};
  double PiNS = 0.0; 
  double piNS[4][4] = {0.0};
  double dmunu[4][4] = {0.0};
  double sigNS[4][4] = {0.0};
  double Delta[4][4] = {0.0};  
  
  double gmunu[4] = { 1.0, -1.0 , -1.0, -1.0/pow(tau-0.5*dt,2) };     
  
  for(int ix = 0; ix < f->get_nx() ; ix++)
    for(int iy = 0; iy <  f->get_ny() ; iy++)
      for(int iz = 0; iz <  f->get_neta() ; iz++)
	{ //loop for cell
	  
	  cell* c = f->get_cell(ix,iy,iz); // get that cell

	  //get energy density, pressure etc. at the cell center at 1/2 time step.
	  c->get_center_varH(eos, tau - 0.5*dt , eps, prs, nb, nq, ns, vx, vy, vz);

	  if (eps < del_eps)
	    {             // empty cell?
	      for (int i = 0; i < 4; i++)
		for (int j = 0; j <= i; j++)
		  {
		    c->set_piH0(i, j, 0.0);
		    c->set_pi0(i, j, 0.0);
		  }
	      c->set_PiH0(0.0);
	      c->set_Pi0(0.0);
	    }
	  
          else

	    {	      
	      
	      //calculates u^{\mu}
	      u[0] = 1.0 / sqrt( 1.0 - vx*vx - vy*vy - vz*vz );
	      u[1] = u[0] * vx;
	      u[2] = u[0] * vy;
	      u[3] = u[0] * vz / ( tau - 0.5*dt );  
	      
	      // \Delta^{\mu \nu}
  for(int mu=0; mu<4; mu++) 
    for(int nu=0; nu<4; nu++)
      {
        if (mu == nu) 
          Delta[mu][nu] = gmunu[mu] - u[mu] * u[nu];
        else
          Delta[mu][nu] = - u[mu] * u[nu];
      }
	      
	      /////////////////////////////////////////////////////
	      //    source term  -> +[Q_{vis}]_{i}^{n}           //
	      /////////////////////////////////////////////////////
	      
	      
	      for(int i =0 ; i<4; i++)
		flux[i] = (tau-dt) * ( c->get_Pi() * u[0] * u[i] + c->get_pi(0,i) );
	      flux[0] -=  (tau-dt) * c->get_Pi();
	      flux[3] *= (tau-dt) ; // [ NOTHING WRONG !!! Q[3] = \tau^{2} T_{vis}^{\eta \tau} ]
	      
	      for(int i=0; i<4; i++)
		{
		  if(std::isinf(flux[i]) or std::isnan(flux[i]))
		    {cout<<"[Error] flux inf in Q_{vis}^n "<<endl; exit(1); }
		}
	      
	      c->add_flux(flux[0],flux[1],flux[2],flux[3],0.,0.,0.); 
	      
	      
	      /////////////////////////////////////////////////////	  
	      //            get relaxation times                 //
	      /////////////////////////////////////////////////////
	      
	      
	      double taupi, tauPi; double T;
	      T = eos->temperature(eps,nb,nq,ns);
	      trcoef->getTau(eps, T, taupi, tauPi);
	      double deltapipi, taupipi, lambdapiPi, phi7, delPiPi, lamPipi;
	      trcoef->getOther(eps, nb, nq, ns, deltapipi, taupipi, lambdapiPi, phi7);
	      trcoef->getOtherBulk(eps, nb, nq, ns, delPiPi, lamPipi);

              if(tauPi<dt){cout<<"[Error] tauPi < dt "<<endl;}
              if(taupi<dt){cout<<"[Error] taupi < dt "<<endl;}
	      
	      ////////////////////////////////////////////////////////////////  
	      // get \partial_{\mu} u^{\nu} , \pi_{NS}^{\mu \nu} & \Pi_{NS} //
	      ////////////////////////////////////////////////////////////////
	      
	      
	      NSquant( ix, iy, iz,dmunu, PiNS, piNS);
	      double theta = dmunu[0][0] + dmunu[1][1] + dmunu[2][2] + dmunu[3][3] ;
	      
	      
	      ////////////////////////////////////////////////////////////////  
	      //     sigmaNS = piNS / (2*eta),                              // 
	      //     protect against division by zero in the eta=0 case.    //
	      ////////////////////////////////////////////////////////////////
	      
	      double etaS, zetaS;
	      trcoef->getEta(eps, T, etaS, zetaS);
	      const double s = eos->entropy(eps, nb, nq, ns);
	      const double eta = etaS * s;
	      for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		  {
		    sigNS[i][j] = 0.5 * piNS[i][j] / eta ;
		    if( eta <= 0.0 ) sigNS[i][j] = 0.0 ;
		  }	  
	      
	      double cov_sigNS[4][4]={0.0};
	      for(int i = 0; i < 4; i++)  
		for(int j = 0; j < 4; j++) 
                  cov_sigNS[i][j] = (1.0/gmunu[i]) * (1.0/gmunu[j]) * sigNS[i][j] ;

	      
	      /////////////////////////////////////////////
	      //     [Let's solve the evolution equation]//
	      /////////////////////////////////////////////
	      
	      
	      // no formal solution is taken [that is not recommended]
	      for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++)
		  {
		    c->set_piH0( i ,j , c->get_pi(i,j) -  (0.5*dt) / ( u[0] * taupi ) * ( c->get_pi(i,j) - piNS[i][j] )  );//pi set
		  } 	
	      c->set_PiH0( c->get_Pi() - (0.5 * dt ) / ( u[0] * tauPi ) * ( c->get_Pi() - PiNS )   ); //Pi set
	      
	      
	      //[source term] -> \delta_{\Pi \Pi} / \tau_{\Pi} \Pi \partial_{\gamma} u^{\gamma}
	      c->add_PiH0( -delPiPi * c->get_Pi() * theta * (0.5*dt/u[0]) );

	      
	      // [source term] -> 8/5 (1/3 - c_s^2 ) \pi^{\mu \nu} \sig_{\mu \nu}	
	      //[after converting contravariant \sigma^{\mu \nu} to covariant one.
	      for(int mu = 0; mu<4; mu++)  
		for(int nu= 0; nu<4; nu++)
		  c->add_PiH0( lamPipi * c->get_pi(mu,nu) * cov_sigNS[mu][nu] * (0.5*dt/u[0]) ) ;
		 
	      
	      
	      // [source term] -> [ u^{\nu} \pi^{\mu \beta} + u^{\mu} \pi^{\nu \beta} ] \times
	      //                    g_{\beta \rho} u^{\gamma} \partial{\gamma} u^{\rho}
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++)
		  for(int beta=0; beta<4; beta++)
		    for(int gamma=0; gamma<4; gamma++)
			{
                            c->add_piH0(mu,nu,-( ( u[nu]*c->get_pi(mu,beta) + 
						 u[mu]*c->get_pi(nu,beta) )* (1.0/gmunu[beta])  //g_{\mu \mu} = 1.0/g^{\mu \mu}
						 * u[gamma]*dmunu[gamma][beta] ) * (0.5*dt/u[0]) );
			}
	      
	      
	      // [source term] -> u^{\gamma} \Gamma^{\mu}_{\gamma \lambda} \pi^{\lambda \nu}
	      //                 -u^{\gamma} \Gamma^{\nu}_{\gamma \lambda} \pi^{ \mu \lambda}  
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++)
		  for(int gamma=0; gamma<4; gamma++)
		    for(int lambda=0; lambda<4; lambda++){
		      c->add_piH0(mu,nu,( -u[gamma] * Gamma(mu,gamma,lambda,tau-0.75*dt) * c->get_pi(lambda,nu) 
					  -u[gamma] * Gamma(nu,gamma,lambda,tau-0.75*dt) * c->get_pi(mu,lambda) ) * (0.5*dt/u[0]) );
		      
		    }
	      
	      //[source term] -> \delta_{\pi \pi} / \tau_{\pi} pi^{\mu \nu} \partial_{\gamma} u^{\gamma}
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++){
		  c->add_piH0(mu, nu, - deltapipi* c->get_pi(mu,nu) * theta * (0.5*dt/u[0]) );
		}
	      
	      
	      
	      // [source term] -> \frac{\phi_7}{\tau_{\pi}} \pi_{\alpha}^{\langle mu} \pi^{\nu \rangle \alpha}
	      //           and -> \frac{\tau_{\pi \pi}}{\tau_{\pi}} \pi_{\alpha}^{\langle \mu} \sigma^{\nu \rangle \alpha }	  
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++)
		  for(int alpha=0; alpha<4; alpha++)
		    for(int lambda=0; lambda<4; lambda++){
		      double pi_sig = 0.5* (  c->get_pi(alpha,mu) * sigNS[nu][alpha] * (1.0/gmunu[alpha])
					      + c->get_pi(alpha,nu) * sigNS[mu][alpha] * (1.0/gmunu[alpha])  )
			- (1.0/3.0) * Delta[mu][nu] * c->get_pi(alpha,lambda) * sigNS[alpha][lambda]
			* (1.0/gmunu[alpha]) * (1.0/gmunu[lambda]) ;  
		      double pi_pi = c->get_pi(alpha,mu) * c->get_pi(nu,alpha) * (1.0/gmunu[alpha])
			- (1.0/3.0) * Delta[mu][nu] * c->get_pi(alpha,lambda) * c->get_pi(alpha,lambda)
			* (1.0/gmunu[alpha]) * (1.0/gmunu[lambda]) ;  
		      c->add_piH0(mu,nu, -taupipi* pi_sig * (0.5*dt/u[0]) ) ;
		      c->add_piH0(mu,nu, (phi7/taupi)* pi_pi * (0.5*dt/u[0]) ) ;         
		    }
	      
	      // [source term] -> \frac{\lambda_{\pi \Pi}}{\tau_{\pi}} \Pi \sigma^{\mu \nu}
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++){
		  c->add_piH0(mu,nu, lambdapiPi * c->get_Pi() * sigNS[mu][nu] * (0.5*dt/u[0]) );
		}	
	     
	      
	      ///////////////////////////////////////////////////////////////////////////////////////////////
	      //                      t+(1/2 * Dt)  ->   t+Dt                                              //
	      ///////////////////////////////////////////////////////////////////////////////////////////////
	      
	      
	      for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++)
		  {
		    c->set_pi0( i ,j , c->get_pi(i,j) -  (dt) / ( u[0] * taupi ) * ( c->get_piH0(i,j) - piNS[i][j] )  );
		  } 	
	      c->set_Pi0( c->get_Pi() - (dt) / ( u[0] * tauPi ) * ( c->get_PiH0() - PiNS )  );  
	      
	      // \Pi's sources
	      c->add_Pi0( -delPiPi * c->get_PiH0() * theta * (dt/u[0]) );
	      
	      for(int mu = 0; mu<4; mu++)  
		for(int nu= 0; nu<4; nu++){
		  c->add_Pi0( lamPipi * c->get_piH0(mu,nu) * cov_sigNS[mu][nu] * (dt/u[0]) ) ;
		}  
	      
	      // \pi's sources
	      
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++)
		  for(int beta=0; beta<4; beta++)
		    for(int gamma=0; gamma<4; gamma++)
			{
                          c->add_pi0(mu,nu,-( ( u[nu]*c->get_piH0(mu,beta) + 
		                          u[mu]*c->get_piH0(nu,beta) )* (1.0/gmunu[beta])  
						 * u[gamma]*dmunu[gamma][beta] ) * (dt/u[0]) );
			}
	      
	      
	      
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++)
		  for(int gamma=0; gamma<4; gamma++)
		    for(int lambda=0; lambda<4; lambda++){
		      c->add_pi0(mu,nu,( -u[gamma] * Gamma(mu,gamma,lambda,tau-0.5*dt) * c->get_piH0(lambda,nu) 
					 -u[gamma] * Gamma(nu,gamma,lambda,tau-0.5*dt) * c->get_piH0(mu,lambda) ) * (dt/u[0]) );
		      
		    }
	      
	     
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++){
		  c->add_pi0(mu, nu, - deltapipi * c->get_piH0(mu,nu) * theta * (dt/u[0]) );
		}
	      
	      
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++)
		  for(int alpha=0; alpha<4; alpha++)
		    for(int lambda=0; lambda<4; lambda++){
		      double pi_sig = 0.5* (  c->get_piH0(alpha,mu) * sigNS[nu][alpha] * (1.0/gmunu[alpha])
					      + c->get_piH0(alpha,nu) * sigNS[mu][alpha] * (1.0/gmunu[alpha])  )
			- (1.0/3.0) * Delta[mu][nu] * c->get_piH0(alpha,lambda) * sigNS[alpha][lambda]
			* (1.0/gmunu[alpha]) * (1.0/gmunu[lambda]) ;  
		      double pi_pi = c->get_pi(alpha,mu) * c->get_piH0(nu,alpha) * (1.0/gmunu[alpha])
			- (1.0/3.0) * Delta[mu][nu] * c->get_piH0(alpha,lambda) * c->get_piH0(alpha,lambda)
			* (1.0/gmunu[alpha]) * (1.0/gmunu[lambda]) ;  
		      c->add_pi0(mu,nu, -taupipi* pi_sig * (dt/u[0]) ) ;
		      c->add_pi0(mu,nu, (phi7/taupi)* pi_pi * (dt/u[0]) ) ;         
		    }
	      
	      for(int mu=0; mu<4; mu++)
		for(int nu=0; nu<=mu; nu++){
		  c->add_pi0(mu,nu, lambdapiPi * c->get_PiH0() * sigNS[mu][nu] * (dt/u[0]) );
		}	  
	      
	    }	
	} //loop for cell
  
  
  
	  /////////////////////////////////////////
	  //             advection               //
	  /////////////////////////////////////////
  
  
  double pi[4][4], piH[4][4], Pi, PiH;
  
  for (int ix = 0; ix < f->get_nx() ; ix++)
    for (int iy = 0; iy <  f->get_ny() ; iy++)
      for (int iz = 0; iz <  f->get_neta() ; iz++)
	{ //advection loop (all cells)
	  cell *c = f->get_cell(ix, iy, iz);
	  c->get_center_varH(eos, tau - 0.5 * dt, eps, prs, nb, nq, ns, vx, vy,
			     vz);  // getPrimVar() before
	  if (eps < del_eps) continue;
	  double xm = -vx * dt / f->get_dx();
	  double ym = -vy * dt / f->get_dy();
	  double zm = -vz * dt /f->get_deta() /  (tau - 0.5 * dt);
	  double xmH = -vx * dt / f->get_dx() / 2.0;
	  double ymH = -vy * dt / f->get_dy() / 2.0;
	  double zmH = -vz * dt / f->get_deta() / 2.0/  (tau - 0.5 * dt);
	  double wx[2] = {(1. - fabs(xm)), fabs(xm)};
	  double wy[2] = {(1. - fabs(ym)), fabs(ym)};
	  double wz[2] = {(1. - fabs(zm)), fabs(zm)};
	  double wxH[2] = {(1. - fabs(xmH)), fabs(xmH)};
	  double wyH[2] = {(1. - fabs(ymH)), fabs(ymH)};
	  double wzH[2] = {(1. - fabs(zmH)), fabs(zmH)};
	  for (int i = 0; i < 4; i++)
	    for (int j = 0; j < 4; j++)
	      {
		pi[i][j] = piH[i][j] = 0.0;
	      }
	  Pi = PiH = 0.0;
	  for (int jx = 0; jx < 2; jx++)
	    for (int jy = 0; jy < 2; jy++)
	      for (int jz = 0; jz < 2; jz++)
		{
		  // pi,Pi-->full step, piH,PiH-->half-step
		  cell *c1 = f->get_cell(ix + jx * sign(xm), iy + jy * sign(ym),
					iz + jz * sign(zm));
		  for (int i = 0; i < 4; i++)
		    for (int j = 0; j < 4; j++) {
		      pi[i][j] += wx[jx] * wy[jy] * wz[jz] * c1->get_pi0(i, j);
		      piH[i][j] += wxH[jx] * wyH[jy] * wzH[jz] * c1->get_piH0(i, j);
		    }
		  Pi += wx[jx] * wy[jy] * wz[jz] * c1->get_Pi0();
		  PiH += wxH[jx] * wyH[jy] * wzH[jz] * c1->get_PiH0();
		}




	  /////////////////////////////////////////////////////////////////////
	  //======= hydro applicability check (viscous corrections limiter): //
          //======= hydro regulation routine (hrr)                           // 
          /////////////////////////////////////////////////////////////////////

          bool rescaled  = false ;
          double quant[6] = {tau-0.5*dt,eps,prs,vx,vy,vz};

          /*
          hr->hrr_vhlle(piH,  PiH, quant, piH,  PiH, rescaled);
          if(rescaled == true ) {c->set_visc_correct_flag(1); }
          hr->hrr_vhlle(pi,  Pi, quant, pi,  Pi, rescaled);
          if(rescaled == true ) {c->set_visc_correct_flag(1); }
          */

          /*
          hr->hrr_music(piH,  PiH, quant, piH,  PiH, rescaled);
          if(rescaled == true ) {c->set_visc_correct_flag(1); }
          hr->hrr_music(pi,  Pi, quant, pi,  Pi, rescaled);
          if(rescaled == true ) {c->set_visc_correct_flag(1); }
          */

          hr->hrr_music2(piH,  PiH, quant, piH,  PiH, rescaled);
          if(rescaled == true ) {c->set_visc_correct_flag(1); }
          hr->hrr_music2(pi,  Pi, quant, pi,  Pi, rescaled);
          if(rescaled == true ) {c->set_visc_correct_flag(1); }



       
   
	  /////////////////////////////////////////
	  // updating to the new values          //
	  /////////////////////////////////////////



	  for (int i = 0; i < 4; i++)
	    for (int j = 0; j <= i; j++)
	      {
		c->set_pi(i, j, pi[i][j]);
		c->set_piH(i, j, piH[i][j]);
	      }
	  c->set_Pi(Pi);
	  c->set_PiH(PiH);


	  /////////////////////////////////////////////////////
	  // source term  - (tau+dt)*delta_Q_(i+1)/delta_tau //
	  /////////////////////////////////////////////////////


	  u[0]= 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
	  u[1] = u[0] * vx;
	  u[2] = u[0] * vy;
	  u[3] = u[0] * vz / (tau - 0.5*dt);


	  for (int i = 0; i < 4; i++)
	    flux[i] = -tau * (c->get_pi(0, i) + c->get_Pi() * u[0] * u[i]); //-[Q_{vis}]_{i}^{n+1}
	  flux[0] += tau * c->get_Pi();
          flux[3] *= tau ;

           for(int i=0; i<4; i++)
              {
               if(std::isinf(flux[i]) or std::isnan(flux[i]))
                 {
                  cout<<"[Error] flux inf in Q_{vis}^{n+1} "<<endl;
                  exit(1); 
                 }
               } 

	  c->add_flux(flux[0], flux[1], flux[2], flux[3], 0., 0., 0.);
	}  // advection loop (all cells)
	  

}





void hydro::visc_flux(cell *left, cell *right, int direction)
{
  int index; double dxa=1.0;
  double el,er,p,nb,nq,ns,vxl,vyl,vzl,vxr,vyr,vzr,u[4],Delta[4][4],flux[7] ;
  double tau1 = tau - 0.5*dt ;
  left->get_center_varH(eos, tau1 , el, p, nb, nq, ns, vxl, vyl, vzl);
  right->get_center_varH(eos, tau1 , er, p, nb, nq, ns, vxr, vyr, vzr);

  if(el < del_eps && er < del_eps) return;

  vxl = 0.5 * (vxl + vxr);
  vyl = 0.5 * (vyl + vyr);
  vzl = 0.5 * (vzl + vzr);
  u[0] = 1.0 / sqrt(1.0 - vxl * vxl - vyl * vyl - vzl * vzl);
  u[1] = u[0] * vxl;
  u[2] = u[0] * vyl;
  u[3] = u[0] * vzl / (tau1);
  
  double g[4][4];    //g^{\mu \nu}
  for(int mu=0; mu<4; mu++)
    for(int nu=0; nu<4; nu++)
      {if (mu != nu) g[mu][nu] = 0.0;}
  g[0][0] = 1.0; g[1][1] = -1.0; g[2][2] = -1.0;
  g[3][3] = -1.0/((tau1)*(tau1));
  
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      Delta[i][j] = g[i][j] - u[i]*u[j] ;
  
  if(direction == X_) {index = 1; dxa = f->get_dx();} 
  if(direction == Y_) {index = 2; dxa = f->get_dy();} 
  if(direction == Z_) {index = 3; dxa = f->get_deta();} 
  
  flux[T_] = (dt/dxa)*tau1*( -0.5*(left->get_PiH()+right->get_PiH())*Delta[T_][index] 
			     + 0.5*(left->get_piH(T_,index)+right->get_piH(T_,index))  );
  flux[X_] = (dt/dxa)*tau1*( -0.5*(left->get_PiH()+right->get_PiH())*Delta[X_][index]
			     + 0.5*(left->get_piH(X_,index)+right->get_piH(X_,index))  );
  flux[Y_] = (dt/dxa)*tau1*( -0.5*(left->get_PiH()+right->get_PiH())*Delta[Y_][index] 
			     + 0.5*(left->get_piH(Y_,index)+right->get_piH(Y_,index))  );
  flux[Z_] = (dt/dxa)*tau1*tau1*( -0.5*(left->get_PiH()+right->get_PiH())*Delta[Z_][index] 
				  + 0.5*(left->get_piH(Z_,index)+right->get_piH(Z_,index))  ); // tau^2 is multiplied here
  
  for(int i=0; i<4; i++)
    {
      if(std::isinf(flux[i]) or std::isnan(flux[i]))
	{cout<<"[Error] Visc flux inf "<<endl; exit(1); }
    } 
  
  left->add_flux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.);
  right->add_flux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.);
}






void hydro::visc_source_step(int ix, int iy, int iz) {
 double e, p, nb, nq, ns, vx, vy, vz;
 double uuu[4];
 double k[7];
 double tau1 = tau - dt / 2.;
	
 cell *c = f->get_cell(ix, iy, iz);

 c->get_center_varH(eos, tau - dt / 2., e, p, nb, nq, ns, vx, vy, vz);
 if (e < del_eps) return;
 uuu[0] = 1. / sqrt(1. - vx * vx - vy * vy - vz * vz);
 uuu[1] = uuu[0] * vx;
 uuu[2] = uuu[0] * vy;
 uuu[3] = uuu[0] * vz/tau1;

 double delta_etaeta = (-1.0/(tau1*tau1))-uuu[3]*uuu[3] ;
 double delta_etatau = -uuu[3]*uuu[0];

 k[T_] = (-tau1*tau1)*(-delta_etaeta*c->get_PiH() + c->get_piH(3,3));
 k[X_] = 0.;
 k[Y_] = 0.;
 k[Z_] = (-tau1)*(-delta_etatau*c->get_PiH() + c->get_piH(3,0));
 for (int i = 0; i < 4; i++) k[i] *= dt;
 c->add_flux(k[T_], k[X_], k[Y_], k[Z_], 0., 0., 0.);
}



void hydro::setNSvalues() {
 double e, p, nb, nq, ns, vx, vy, vz, piNS[4][4], PiNS;
 for (int ix = 0; ix < f->get_nx(); ix++)
  for (int iy = 0; iy < f->get_ny(); iy++)
   for (int iz = 0; iz < f->get_neta(); iz++) {
    cell *c = f->get_cell(ix, iy, iz);
    c->get_physical_var(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
    if (e <= 0.) continue;
    // NSquant(ix, iy, iz, piNS, PiNS, dmu, du) ;
    //############## set NS values assuming initial zero flow + Bjorken z
    // flow
    double T;
    double etaS, zetaS;
    double s = eos->entropy(e, nb, nq, ns);  // entropy density in the current cell
    //eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
    T = eos->temperature(e,nb,nq,ns);
    trcoef->getEta(e, T, etaS, zetaS);
    for (int i = 0; i < 4; i++)
     for (int j = 0; j < 4; j++) piNS[i][j] = 0.0;  // reset piNS
    piNS[1][1] = piNS[2][2] = 2.0 / 3.0 * etaS * s / tau / 5.068;
    piNS[3][3] = -2.0/pow(tau,2) * piNS[1][1];
    PiNS = 0.0;
    for (int i = 0; i < 4; i++)
     for (int j = 0; j <= i; j++) c->set_pi(i, j, piNS[i][j]);
    c->set_Pi(PiNS);
   }
 cout << "setNS done\n";
}








