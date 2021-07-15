#include "grid.h"
#include "global.h"

grid::grid(idb* _IDB, cnvrt* _CN, trancoeff* _tr, EoS* _eos)
{
  IDB = _IDB ;
  CN = _CN;
  trcoef = _tr ;
  eos = _eos;
}

grid::~grid()
{
 delete[] Cell; 
}

void grid::make_grid()
{
  
  if(IDB->ic_mode == 2)
    {
      double dummy2;
      string dummy;
      double temp_nx, temp_ny ; // x and y from file.
      double xf, yf;
      int neta_;
      double deta_,dx_,dy_;
      
      
      std::fstream ic_file;
      ic_file.open(IDB->init_file_name.c_str(),ios::in);
      if (!ic_file) 
	{
	  cout<<"couldn't find ic file."<<endl;
	  exit(1);
	}
      
      ic_file.getline(buff,200);
      
      iss = new istringstream(buff);
      
      *iss >> dummy >> dummy >> dummy2 
	   >> dummy >> neta_ >> dummy >> temp_nx >> dummy >> temp_ny 
	   >> dummy >> deta_ >> dummy >> dx_ >> dummy >> dy_ ;
      
      IDB->nx = temp_nx ;
      IDB->ny = temp_ny ;
      
      delete iss;
      
      ic_file.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> dummy2 >> xf >> yf 
	   >> dummy2 >> dummy2 >> dummy2 >> dummy2
	   >> dummy2  >> dummy2  >> dummy2  >> dummy2;
      
      delete iss;
      IDB->xmin = xf ;
      IDB->xmax = fabs(xf) ;
      IDB->ymin = yf ;
      IDB->ymax = fabs(yf) ;
      
      IDB->dx =  ( IDB->xmax- IDB->xmin )
	         /  ( IDB->nx - 1 ) ;
      IDB->dy =  ( IDB->ymax- IDB->ymin )
	         / ( IDB->ny - 1 ) ;
      
      cout<<"[Info] grid resized."<<endl;
      cout<<"[Info] xmin : "<<IDB->xmin<<"  xmax : "<<IDB->xmax<<endl;
      cout<<"[Info] ymin : "<<IDB->ymin<<"  ymax : "<<IDB->ymax<<endl;
      cout<<"[Info] dx : "<<IDB->dx<<endl;
      cout<<"[Info] dy : "<<IDB->dy<<endl;
      
    }


  if (IDB->dtau > IDB->dx / 2. ||
      IDB->dtau > IDB->dy / 2. )
     {
      cout << "[Error]  too big delta_tau : " << IDB->dtau <<
              "  " << IDB->dx << "  " << IDB->dtau << endl;
      exit(1);
     }

 if(IDB->neta > 1 && 
    IDB->dtau > IDB->tau0*IDB->deta/2. )
   {
     cout << "[Error]  too big delta_tau : " << IDB->dtau <<
          "  tau*deta/2. = "<< IDB->tau0*IDB->deta/2. << endl; 
     exit(1);
   }
  
  
  cout<<"[Info] nx = "<<IDB->nx<<endl;
  cout<<"[Info] ny = "<<IDB->ny<<endl;
  cout<<"[info] neta = "<<IDB->neta<<"\n"<<endl;



  nx = IDB->nx ;
  ny = IDB->ny ;
  neta = IDB->neta ;

  dx = IDB->dx ;
  dy = IDB->dy ;
  deta = IDB->deta ;

  xmin = IDB->xmin;
  xmax = IDB->xmax;
  ymin = IDB->ymin;
  ymax = IDB->ymax;
  etamin = IDB->etamin;
  etamax = IDB->etamax;
  
  dtau = IDB->dtau ;
  
  
  
  // now grid properties are fixed and it's ready to be made
  Cell = new cell[nx * ny * neta];



/*  
  for(int i = 0; i < IDB->nx; i++)
    {
      vector < vector < cell > > w;
      tube.push_back( w );
      for(int j = 0; j < IDB->ny; j++)
	{
	  vector <cell> v;
	  tube[i].push_back( v );
	  for(int k = 0; k < IDB->neta; k++)
	    {
	      //cout<<"i = "<<i << " j = " << j << " k = " << k << endl;
	      tube[i][j].push_back(*new_cell(i,j,k));
	      //c1->set_pos(i,j,k);
	      //cout<<new_cell(i,j,k)<<endl;
	    }
	}
    }
*/ 
  
  
  for (int iz = 0; iz < IDB->neta; iz++)
    for (int iy = 0; iy < IDB->ny; iy++)
      for (int ix = 0; ix < IDB->nx; ix++)
	{
          //cout<<"iy = "<<iy << " ix = " << ix <<endl;
          get_cell(ix, iy, iz)->set_pos(ix,iy,iz);
          get_cell(ix, iy, iz)->set_cnvrt(CN);
	  get_cell(ix, iy, iz)->set_prev_cell(X_, get_cell(ix - 1, iy, iz));
	  get_cell(ix, iy, iz)->set_next_cell(X_, get_cell(ix + 1, iy, iz));
	  get_cell(ix, iy, iz)->set_prev_cell(Y_, get_cell(ix, iy - 1, iz));
	  get_cell(ix, iy, iz)->set_next_cell(Y_, get_cell(ix, iy + 1, iz));
	  get_cell(ix, iy, iz)->set_prev_cell(Z_, get_cell(ix, iy, iz - 1));
	  get_cell(ix, iy, iz)->set_next_cell(Z_, get_cell(ix, iy, iz + 1));    
	}
  
}


cell* grid::new_cell(int ix,int iy,int iz)
{
  cell* c1 = new cell();
  c1->set_pos(ix,iy,iz);
  return c1;
}





void grid::correct_imaginary_cells()
{

  double Q[7]={0.0};
  // Z
 for (int ix = 0; ix < IDB->nx; ix++)
  for (int iy = 0; iy < IDB->ny; iy++) {
   // left boundary
   get_cell(ix, iy, 2)->get_Q(Q);
   get_cell(ix, iy, 1)->set_Q(Q);
   get_cell(ix, iy, 0)->set_Q(Q);
   // right boundary
   get_cell(ix, iy, IDB->neta - 3)->get_Q(Q);
   get_cell(ix, iy, IDB->neta - 2)->set_Q(Q);
   get_cell(ix, iy, IDB->neta - 1)->set_Q(Q);
  }
 // Y
 for (int ix = 0; ix < IDB->nx; ix++)
  for (int iz = 0; iz < IDB->neta; iz++) {
   // left boundary
   get_cell(ix, 2, iz)->get_Q(Q);
   get_cell(ix, 1, iz)->set_Q(Q);
   get_cell(ix, 0, iz)->set_Q(Q);
   // right boundary
   get_cell(ix, IDB->ny - 3, iz)->get_Q(Q);
   get_cell(ix, IDB->ny - 2, iz)->set_Q(Q);
   get_cell(ix, IDB->ny - 1, iz)->set_Q(Q);
  }
 // X
 for (int iy = 0; iy < IDB->ny; iy++)
  for (int iz = 0; iz < IDB->neta; iz++) {
   // left boundary
   get_cell(2, iy, iz)->get_Q(Q);
   get_cell(1, iy, iz)->set_Q(Q);
   get_cell(0, iy, iz)->set_Q(Q);
   // right boundary
   get_cell(IDB->nx - 3, iy, iz)->get_Q(Q);
   get_cell(IDB->nx - 2, iy, iz)->set_Q(Q);
   get_cell(IDB->nx - 1, iy, iz)->set_Q(Q);
  }
 
}



void grid::correct_imaginary_cells_full() {
 double Q[7], _pi[4][4], _Pi;
 // Z
 for (int ix = 0; ix < IDB->nx; ix++)
  for (int iy = 0; iy < IDB->ny; iy++) {
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
   get_cell(ix, iy, IDB->neta - 3)->get_Q(Q);
   get_cell(ix, iy, IDB->neta - 2)->set_Q(Q);
   get_cell(ix, iy, IDB->neta - 1)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = get_cell(ix, iy, IDB->neta - 3)->get_pi(i, j);
   _Pi = get_cell(ix, iy, IDB->neta - 3)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(ix, iy, IDB->neta - 2)->set_pi(i, j, _pi[i][j]);
     get_cell(ix, iy, IDB->neta - 1)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(ix, iy, IDB->neta - 2)->set_Pi(_Pi);
   get_cell(ix, iy, IDB->neta - 1)->set_Pi(_Pi);
  }
 // Y
 for (int ix = 0; ix < IDB->nx; ix++)
  for (int iz = 0; iz < IDB->neta; iz++) {
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
   get_cell(ix, IDB->ny - 3, iz)->get_Q(Q);
   get_cell(ix, IDB->ny - 2, iz)->set_Q(Q);
   get_cell(ix, IDB->ny - 1, iz)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = get_cell(ix, IDB->ny - 3, iz)->get_pi(i, j);
   _Pi = get_cell(ix, IDB->ny - 3, iz)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(ix, IDB->ny - 2, iz)->set_pi(i, j, _pi[i][j]);
     get_cell(ix, IDB->ny - 1, iz)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(ix, IDB->ny - 2, iz)->set_Pi(_Pi);
   get_cell(ix, IDB->ny - 1, iz)->set_Pi(_Pi);
  }
 // X
 for (int iy = 0; iy < IDB->ny; iy++)
  for (int iz = 0; iz < IDB->neta; iz++) {
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
   get_cell(IDB->nx - 3, iy, iz)->get_Q(Q);
   get_cell(IDB->nx - 2, iy, iz)->set_Q(Q);
   get_cell(IDB->nx - 1, iy, iz)->set_Q(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = get_cell(IDB->nx - 3, iy, iz)->get_pi(i, j);
   _Pi = get_cell(IDB->nx - 3, iy, iz)->get_Pi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     get_cell(IDB->nx - 2, iy, iz)->set_pi(i, j, _pi[i][j]);
     get_cell(IDB->nx - 1, iy, iz)->set_pi(i, j, _pi[i][j]);
    }
   get_cell(IDB->nx - 2, iy, iz)->set_Pi(_Pi);
   get_cell(IDB->nx - 1, iy, iz)->set_Pi(_Pi);
  }
}






