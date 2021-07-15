#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include "idb.h"

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::cin;
using std::fstream;
using std::ios;
using std::istringstream;

// read_id : read 'I'nput 'D'ata

class read_id{

public :

read_id(){};
~read_id(){};

// This functions reads the input.dat file and sets the input parameters in the code 
void read_id_from_file(idb *input_parameter_list, string input_file_name)
{
  string a_; char a[50];
  
  istringstream* iss;
  char   buff[200];
  
  fstream File0;
  File0.open(input_file_name,ios::in);
  if(!File0){cout<<"No input file, exit..."; exit(1);}
  
  int number = 0;
  while(!File0.eof())
    {
      File0.getline(buff,200);
      if (!(*buff) || (*buff == '#')) {number ++; continue;}
      iss = new istringstream(buff);
      *iss >> a_ >> a ;
      if(a_ == "nx" )    {input_parameter_list->nx = atof(a);}
      if(a_ == "ny" )    {input_parameter_list-> ny = atof(a);}
      if(a_ == "neta" )  {input_parameter_list->neta = atof(a);}
      if(a_ == "xmin" )  {input_parameter_list->xmin = atof(a);}
      if(a_ == "xmax" )  {input_parameter_list->xmax = atof(a);}
      if(a_ == "ymin" )  {input_parameter_list->ymin = atof(a);}
      if(a_ == "ymax" )  {input_parameter_list->ymax = atof(a);}
      if(a_ == "etamin" )  {input_parameter_list->etamin = atof(a);}
      if(a_ == "etamax" )  {input_parameter_list->etamax = atof(a);}
      if(a_ == "dtau" )  {input_parameter_list->dtau = atof(a);}
      if(a_ == "tau0" )  {input_parameter_list->tau0 = atof(a);}
      if(a_ == "Tfreeze" )  {input_parameter_list->Tfreeze = atof(a);}
      if(a_ == "tauMax" )  {input_parameter_list->tauMax = atof(a);}
      if(a_ == "ic_mode" )  {input_parameter_list->ic_mode = atof(a);}

      if(a_ == "SNN" )  {input_parameter_list->SNN = atof(a);}

      if(a_ == "eos" )  {input_parameter_list->eos = atof(a);}

      if(a_ == "bmin" )  {input_parameter_list->bmin = atof(a);}  
      if(a_ == "bmax" )  {input_parameter_list->bmax = atof(a);}  

      if(a_ == "eta_platue" )  {input_parameter_list->eta_platue = atof(a);} 
      if(a_ == "eta_fall" )  {input_parameter_list->eta_fall = atof(a);} 

      if(a_ == "skip_fo_tau" )  {input_parameter_list->skip_fo_tau = atof(a);}   
      if(a_ == "skip_fo_x" )  {input_parameter_list->skip_fo_x = atof(a);}   
      if(a_ == "skip_fo_y" )  {input_parameter_list->skip_fo_y = atof(a);}   
      if(a_ == "skip_fo_eta" )  {input_parameter_list->skip_fo_eta = atof(a);}   

      if(a_ == "save_every_N_steps" )  {input_parameter_list->save_every_N_steps = atof(a);}   

      if(input_parameter_list->ic_mode == 0) //Optical-Glauber inputs
        {
          if(a_ == "opt_eps0" )  {input_parameter_list->opt_eps0 = atof(a);} 
          if(a_ == "species" )    {input_parameter_list->species = a;}      
        }     
      
      if(input_parameter_list->ic_mode == 1) //MC-Glauber inputs
        {
          if(a_ == "projectile" )    {input_parameter_list->projectile = a;}
          if(a_ == "target" )    {input_parameter_list->target = a;}
          if(a_ == "DELTA" )  {input_parameter_list->DELTA = atof(a);}       
          if(a_ == "mc_eps0" )  {input_parameter_list->mc_eps0 = atof(a);}       
        }     
      
          if(a_ == "init_file_name" )    {input_parameter_list->init_file_name = a;}      

          if(a_ == "etas_flag" )    {input_parameter_list->etas_flag = atof(a);} 
          if(a_ == "etas" )  {input_parameter_list->etas = atof(a);}     
          if(a_ == "zetas_flag" )    {input_parameter_list->zetas_flag = atof(a);}      
          if(a_ == "t_etas_flag" )    {input_parameter_list->t_etas_flag = atof(a);}      

      delete iss;
      number++;
    } 

   // set dx,dy and deta by calculating it from nx, ny, neta, xmin, xmax .... etc.
    input_parameter_list->dx =  ( input_parameter_list->xmax- input_parameter_list->xmin )
                                     /  ( input_parameter_list->nx - 1 ) ;
    input_parameter_list->dy =  ( input_parameter_list->ymax- input_parameter_list->ymin )
                                      /  ( input_parameter_list->ny - 1 ) ;

   if (input_parameter_list->neta == 1) { input_parameter_list->deta = 1e-15; } //[important] for 2+1D deta is a very small no. not zero.
   else  { input_parameter_list->deta =   ( input_parameter_list->etamax- input_parameter_list->etamin )
                                      /  ( input_parameter_list->neta - 1 ); }


  
  File0.close();

} 
};

