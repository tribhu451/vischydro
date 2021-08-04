#include "eos2.h"

EoS2::~EoS2()
{
  
  for (int itable = 0; itable < 7 ; itable++)
    {
      mtx_free(mu_B_tb[itable],
	       nb_length[itable], e_length[itable]);
      mtx_free(pressure_tb[itable],
	       nb_length[itable], e_length[itable]);
      mtx_free(temperature_tb[itable],
	       nb_length[itable], e_length[itable]);
    }
  
  
  delete [] mu_B_tb;
  delete [] pressure_tb;
  delete [] temperature_tb;
}



EoS2::EoS2()
{
  cout << "[Info] Using NEOS_B EoS from Phys. Rev. C100, 024907 (2019) [arXiv:1902.05095] " << endl;
  
  const int ntables = 7;
  string eos_file_string_array[ntables] = {"1", "2", "3", "4", "5", "6", "7"};
  

  pressure_tb    = new double** [ntables];
  temperature_tb = new double** [ntables];
  mu_B_tb        = new double** [ntables];
  
  
  for (int itable = 0; itable < ntables; itable++)
    {
      std::ifstream eos_p("eos/neos_b/neos" + eos_file_string_array[itable]
			  + "_p.dat");
      std::ifstream eos_T("eos/neos_b/neos" + eos_file_string_array[itable]
			  + "_t.dat");
      std::ifstream eos_mub("eos/neos_b/neos" + eos_file_string_array[itable]
			    + "_mub.dat");
      
      
      if(!eos_p){cout << "eos/neos_b/neos"<< eos_file_string_array[itable]
		      <<"_p.dat" <<" not found, code terminated. exiting ...."<<endl; exit(1);}
      if(!eos_T){cout << "eos/neos_b/neos" << eos_file_string_array[itable]
		      << "_t.dat" <<" not found, code terminated. exiting ...."<<endl; exit(1);}
      if(!eos_mub){cout <<  "eos/neos_b/neos" << eos_file_string_array[itable]
			<< "_mub.dat" <<" not found, code terminated. exiting ...."<<endl; exit(1);}
      

      // read the first two lines with general info:
      // first value of rhob, first value of epsilon
      // deltaRhob, deltaE, number of rhob points, number of epsilon points
      // the table size is
      // (number of rhob points + 1, number of epsilon points + 1)
      int N_e, N_rhob;
      eos_p >> nb_bounds[itable] >> e_bounds[itable];
      eos_p >> nb_spacing[itable] >> e_spacing[itable]
	    >> N_rhob >> N_e;
      nb_length[itable] = N_rhob + 1;
      e_length[itable]  = N_e + 1;
      
      // e_bounds[itable]  /= Util::hbarc;   // 1/fm^4
      // e_spacing[itable] /= Util::hbarc;   // 1/fm^4
      
      // skip the header in T and mu_B files
      string dummy;
      std::getline(eos_T, dummy);
      std::getline(eos_T, dummy);
      std::getline(eos_mub, dummy);
      std::getline(eos_mub, dummy);
      
      
      // allocate memory for EOS arrays
      pressure_tb[itable] = mtx_malloc(nb_length[itable],
				       e_length[itable]);
      temperature_tb[itable] = mtx_malloc(nb_length[itable],
					  e_length[itable]);
      mu_B_tb[itable] = mtx_malloc(nb_length[itable],
				   e_length[itable]);
      
      
      // read pressure, temperature and chemical potential values
      for (int j = 0; j < e_length[itable]; j++)
	{
	  for (int i = 0; i < nb_length[itable]; i++)
	    {
	      eos_p >> pressure_tb[itable][i][j];
	      eos_T >> temperature_tb[itable][i][j];
	      eos_mub >> mu_B_tb[itable][i][j];
	    }
	}
    }
  
  
}




//       *  *  *  *  *  *  *  *  *  *  *  *  *  *     //
//       *  *  *  *  *  *  *  *  *  *  *  *  *  *     //
//       *  *  *  *  *  *  *  *  *  *  *  *  *  *     //
//                          ..
//                          ..
//                          ..
//                          ..
//                          ..
//                       .........
//                         .....
//                           .



void EoS2::check_eos()
{
}


double EoS2::pressure( double e,double _nb, double _nq, double _ns)
{
    int table_idx = get_table_idx(e);
    double f = interpolate2D(e, std::abs(_nb), table_idx, pressure_tb);
    f = std::max(1e-15, f);
    return(f);
}


double EoS2::temperature( double e,double _nb, double _nq, double _ns)
{
    int table_idx = get_table_idx(e);
    double T = interpolate2D(e, std::abs(_nb), table_idx,
                             temperature_tb);  
    T = std::max(1e-15, T);
    return(T);
}


double EoS2::mu_b( double e,double _nb, double _nq, double _ns)
{
    int table_idx = get_table_idx(e);
    double sign = _nb/(std::abs(_nb) + 1e-15);
    double mu = sign*interpolate2D(e, std::abs(_nb), table_idx,
                                   mu_B_tb);  
    return(mu);
}



double EoS2::entropy( double e,double _nb, double _nq, double _ns)
{   
    auto P    = pressure(e, _nb, _nq ,_nq );
    auto T    = temperature(e, _nb, _nq, _nq);
    auto muB  = mu_b(e, _nb, _nq, _ns);
    auto muS  = 0 ;
    auto muC  = 0 ;
    auto rhoS = 0 ;
    auto rhoC = 0 ;
    auto f    = (e + P - muB*_nb - muS*rhoS - muC*rhoC)/(T + 1e-15);
    return(std::max(1e-15, f));
}


double EoS2::cs2(double e,double _nb, double _nq, double _ns)
{
   double eLeft = 0.9*e;
   double eRight = 1.1*e;
   double pL = pressure(eLeft, _nb, _nq, _ns);   
   double pR = pressure(eRight, _nb, _nq, _ns); 
   double dpde = (pR - pL)/(eRight - eLeft);

   int table_idx = get_table_idx(e);
   double deltaRhob = nb_spacing[table_idx];    
   double rhobLeft  = _nb - deltaRhob*0.5;
   double rhobRight = _nb + deltaRhob*0.5;
          pL = pressure(e, rhobLeft, _nq, _ns);      
          pR = pressure(e, rhobRight, _nq, _ns);   
   double dpdrho = (pR - pL)/(rhobRight - rhobLeft);

    double v_min = 0.01;
    double v_max = 1./3;
    double P = pressure(e, _nb, _nq, _ns);
    double v_sound = dpde + _nb/(e + P + 1e-15)*dpdrho;
    v_sound = std::max(v_min, std::min(v_max, v_sound));
    return(v_sound);

}


double EoS2::temp_2_eps( double T,double _nb, double _nq, double _ns)
{

    double eps_lower = 1e-15;
    double eps_upper = e_bounds[6] + e_spacing[6]*e_length[6];
    double eps_mid   = (eps_upper + eps_lower)/2.;
    double T_lower   = temperature(eps_lower, _nb, _nq, _ns);
    double T_upper   = temperature(eps_upper, _nb, _nq, _ns);
    int ntol         = 1000;
    if (T < 0.0 || T > T_upper) {
        cout << "get_T2e:: T is out of bound, "
             << "T = " << T << ", T_upper = " << T_upper
             << ", T_lower = " << T_lower << endl;
        exit(1);
    }
    if (T < T_lower) return(eps_lower);

    double rel_accuracy = 1e-8;
    double abs_accuracy = 1e-15;
    double T_mid;
    int iter = 0;
    while (((eps_upper - eps_lower)/eps_mid > rel_accuracy
            && (eps_upper - eps_lower) > abs_accuracy) && iter < ntol) {
        T_mid = temperature(eps_mid, _nb, _nq, _ns);
        if (T < T_mid)
            eps_upper = eps_mid;
        else 
            eps_lower = eps_mid;
        eps_mid = (eps_upper + eps_lower)/2.;
        iter++;
    }
    if (iter == ntol) {
        cout << "get_T2e_finite_rhob:: max iteration reached, "
             << "T = " << T << ", rhob = " << _nb << endl;;
        cout << "T_upper = " << T_upper
             << " , T_lower = " << T_lower << endl;
        cout << "eps_upper = " << eps_upper
             << " , eps_lower = " << eps_lower
             << ", diff = " << (eps_upper - eps_lower) << endl;
        exit(1);
    }
    return (eps_mid);

}



double EoS2::entr_2_eps(double s,double _nb, double _nq, double _ns)
{   
    double eps_lower = 1e-15;
    double eps_upper = e_bounds[6] + e_spacing[6]*e_length[6];
    double eps_mid   = (eps_upper + eps_lower)/2.;
    double s_lower   = entropy(eps_lower, _nb, _nq, _ns);
    double s_upper   = entropy(eps_upper, _nb, _nq, _ns);
    int ntol         = 1000;
    if (s < 0.0 || s > s_upper) {
        cout << "get_s2e_finite_rhob:: s is out of bound, "
             << "s = " << s << ", s_upper = " << s_upper
             << ", s_lower = " << s_lower << endl;
        exit(1);
    }

    if (s < s_lower) return(eps_lower);

    double rel_accuracy = 1e-8;
    double abs_accuracy = 1e-15;
    double s_mid;
    int iter = 0;
    while (((eps_upper - eps_lower)/eps_mid > rel_accuracy
            && (eps_upper - eps_lower) > abs_accuracy) && iter < ntol) {
        s_mid = entropy(eps_mid, _nb, _nq, _ns);
        if (s < s_mid)
            eps_upper = eps_mid;
        else 
            eps_lower = eps_mid;
        eps_mid = (eps_upper + eps_lower)/2.;
        iter++;
    }
    if (iter == ntol) {
        cout << "get_s2e_finite_rhob:: max iteration reached, "
             << "s = " << s << ", rhob = " << _nb << endl;
        cout << "s_upper = " << entropy(eps_upper, _nb,_nq,_ns)
             << " , s_lower = " << entropy(eps_lower, _nb,_nq,_ns) << endl;
        cout << "eps_upper = " << eps_upper
             << " , eps_lower = " << eps_lower
             << ", diff = " << (eps_upper - eps_lower) << endl;
        exit(1);
    }
    return (eps_mid);

}










