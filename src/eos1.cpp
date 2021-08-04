#include "eos1.h"



void EoS1::is_file_exist(fstream& file)
{
  if(!file){cout<<"Files for EoS not found, code terminated. exiting ...."<<endl; exit(1);}
  //else{cout<<"[Get] EoS file "<<endl;}
}

EoS1::~EoS1()
{

    gsl_spline_free (spline_temp_1);
    gsl_spline_free (spline_temp_2);
    gsl_spline_free (spline_temp_3);
    gsl_spline_free (spline_temp_4);
    gsl_spline_free (spline_temp_5);
    gsl_spline_free (spline_temp_6);
    gsl_spline_free (spline_temp_7);


    gsl_spline_free (spline_ntrpy_1);
    gsl_spline_free (spline_ntrpy_2);
    gsl_spline_free (spline_ntrpy_3);
    gsl_spline_free (spline_ntrpy_4);
    gsl_spline_free (spline_ntrpy_5);
    gsl_spline_free (spline_ntrpy_6);
    gsl_spline_free (spline_ntrpy_7);

    gsl_spline_free (spline_prs_1);
    gsl_spline_free (spline_prs_2);
    gsl_spline_free (spline_prs_3);
    gsl_spline_free (spline_prs_4);
    gsl_spline_free (spline_prs_5);
    gsl_spline_free (spline_prs_6);
    gsl_spline_free (spline_prs_7);

    gsl_spline_free (spline_teos_e);
    gsl_spline_free (spline_teos_p);
    gsl_spline_free (spline_teos_s);
    gsl_spline_free (spline_teos_t);

    gsl_interp_accel_free (acc);

}

EoS1::EoS1()
{

    cout<<"\n[Info] allocating s95p-v1 EoS from arXiv:0912.2541\n"<<endl;

  // There are 14(7-7) files . These files contain entropy densty, prs densty etc. in equal energy denst. intervals.
  
  infile1_d.open("eos/EoS_s95p-v1/s95p-v1_dens1.dat",ios::in);
  infile2_d.open("eos/EoS_s95p-v1/s95p-v1_dens2.dat",ios::in);
  infile3_d.open("eos/EoS_s95p-v1/s95p-v1_dens3.dat",ios::in);
  infile4_d.open("eos/EoS_s95p-v1/s95p-v1_dens4.dat",ios::in);
  infile5_d.open("eos/EoS_s95p-v1/s95p-v1_dens5.dat",ios::in);
  infile6_d.open("eos/EoS_s95p-v1/s95p-v1_dens6.dat",ios::in);
  infile7_d.open("eos/EoS_s95p-v1/s95p-v1_dens7.dat",ios::in);

  infile1_t.open("eos/EoS_s95p-v1/s95p-v1_par1.dat",ios::in);    
  infile2_t.open("eos/EoS_s95p-v1/s95p-v1_par2.dat",ios::in);    
  infile3_t.open("eos/EoS_s95p-v1/s95p-v1_par3.dat",ios::in);    
  infile4_t.open("eos/EoS_s95p-v1/s95p-v1_par4.dat",ios::in);    
  infile5_t.open("eos/EoS_s95p-v1/s95p-v1_par5.dat",ios::in);    
  infile6_t.open("eos/EoS_s95p-v1/s95p-v1_par6.dat",ios::in);    
  infile7_t.open("eos/EoS_s95p-v1/s95p-v1_par7.dat",ios::in);    

  file_teosdens.open("eos/EoS_s95p-v1/s95p-v1_Teosdens.dat",ios::in);    

  is_file_exist(infile1_d);
  is_file_exist(infile2_d);
  is_file_exist(infile3_d);
  is_file_exist(infile4_d);
  is_file_exist(infile5_d);
  is_file_exist(infile6_d);
  is_file_exist(infile7_d);

  is_file_exist(infile1_t);
  is_file_exist(infile2_t);
  is_file_exist(infile3_t);
  is_file_exist(infile4_t);
  is_file_exist(infile5_t);
  is_file_exist(infile6_t);
  is_file_exist(infile7_t);

  is_file_exist(file_teosdens);

  for(int i=0; i<750; i++)  //reading v1_Teosdens file
    { 

      file_teosdens.getline(buff,200);
      if (!(*buff) || (*buff == '#')) { continue;}
      iss = new istringstream(buff);
      *iss >> teos_e[i] >> teos_p[i] >> teos_s[i] >> nb >>teos_t[i] 
           >> nb >> nb ;
        teos_t[i] *= 0.001 ;
      delete iss;
    } 



  infile1_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_1;
  delete iss;
  
  infile1_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_1 >> ne_1 ;
  delete iss;

  infile2_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_2;
  delete iss;
  
  infile2_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_2 >> ne_2 ;
  delete iss;

  infile3_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_3;
  delete iss;
  
  infile3_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_3 >> ne_3 ;
  delete iss;

  infile4_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_4;
  delete iss;
  
  infile4_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_4 >> ne_4 ;
  delete iss;

  infile5_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_5;
  delete iss;
  
  infile5_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_5 >> ne_5 ;
  delete iss;

  infile6_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_6;
  delete iss;
  
  infile6_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_6 >> ne_6 ;
  delete iss;

  infile7_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> e0_7;
  delete iss;
  
  infile7_d.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> de_7 >> ne_7 ;
  delete iss;


  
  for(int i=ne_1-1; i>=0; i--)  //reverse counting
    { 
      //if (!(*buff) || (*buff == '#')) {number ++; continue;}
      infile1_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_1[i] >> prs_1[i] >> ent_1[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_2-1; i>=0; i--)
    { 
      infile2_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_2[i] >> prs_2[i] >> ent_2[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_3-1; i>=0; i--)
    { 
      infile3_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_3[i] >> prs_3[i] >> ent_3[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_4-1; i>=0; i--)
    { 
      infile4_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_4[i] >> prs_4[i] >> ent_4[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_5-1; i>=0; i--)
    { 
      infile5_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_5[i] >> prs_5[i] >> ent_5[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_6-1; i>=0; i--)
    { 
      infile6_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_6[i] >> prs_6[i] >> ent_6[i] >> nb >> mub;
      delete iss;
    } 

  for(int i=ne_7-1; i>=0; i--)
    { 
      infile7_d.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> eps_7[i] >> prs_7[i] >> ent_7[i] >> nb >> mub;
      delete iss;
    } 

  infile1_t.getline(buff,200);
  infile2_t.getline(buff,200);
  infile3_t.getline(buff,200);
  infile4_t.getline(buff,200);
  infile5_t.getline(buff,200);
  infile6_t.getline(buff,200);
  infile7_t.getline(buff,200);

  infile1_t.getline(buff,200);
  infile2_t.getline(buff,200);
  infile3_t.getline(buff,200);
  infile4_t.getline(buff,200);
  infile5_t.getline(buff,200);
  infile6_t.getline(buff,200);
  infile7_t.getline(buff,200);
  
  for(int i=ne_1-1; i>=0; i--) // reverse counting
    { 
      infile1_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_1[i]  >> mus >> muq;
      delete iss;
    } 

  for(int i=ne_2-1; i>=0; i--)
    { 
      infile2_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_2[i]  >> mus >> muq;
      delete iss;
    } 

  for(int i=ne_3-1; i>=0; i--)
    { 
      infile3_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_3[i]  >> mus >> muq;
      delete iss;
    } 
  for(int i=ne_4-1; i>=0; i--)
    { 
      infile4_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_4[i]  >> mus >> muq;
      delete iss;
    } 
  for(int i=ne_5-1; i>=0; i--)
    { 
      infile5_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_5[i]  >> mus >> muq;
      delete iss;
    } 

  for(int i=ne_6-1; i>=0; i--)
    { 
      infile6_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_6[i]  >> mus >> muq;
      delete iss;
    } 
  for(int i=ne_7-1; i>=0; i--)
    { 
      infile7_t.getline(buff,200);
      iss = new istringstream(buff);
      *iss >> tmp_7[i]  >> mus >> muq;
      delete iss;
    } 


 infile1_d.close();
 infile2_d.close();
 infile3_d.close();
 infile4_d.close();
 infile5_d.close();
 infile6_d.close();
 infile7_d.close();

 infile1_t.close();
 infile2_t.close();
 infile3_t.close();
 infile4_t.close();
 infile5_t.close();
 infile6_t.close();
 infile7_t.close();

 file_teosdens.close();
 // File reading done .



   acc = gsl_interp_accel_alloc();



   // Initialize temperature
   spline_temp_1 = gsl_spline_alloc (gsl_interp_linear, ne_1);
   gsl_spline_init (spline_temp_1, eps_1, tmp_1, ne_1);

   spline_temp_2 = gsl_spline_alloc (gsl_interp_linear, ne_2);
   gsl_spline_init (spline_temp_2, eps_2, tmp_2, ne_2);

   spline_temp_3 = gsl_spline_alloc (gsl_interp_linear, ne_3);
   gsl_spline_init (spline_temp_3, eps_3, tmp_3, ne_3);

   spline_temp_4 = gsl_spline_alloc (gsl_interp_linear, ne_4);
   gsl_spline_init (spline_temp_4, eps_4, tmp_4, ne_4);

   spline_temp_5 = gsl_spline_alloc (gsl_interp_linear, ne_5);
   gsl_spline_init (spline_temp_5, eps_5, tmp_5, ne_5);

   spline_temp_6 = gsl_spline_alloc (gsl_interp_linear, ne_6);
   gsl_spline_init (spline_temp_6, eps_6, tmp_6, ne_6);

   spline_temp_7 = gsl_spline_alloc (gsl_interp_linear, ne_7);
   gsl_spline_init (spline_temp_7, eps_7, tmp_7, ne_7);


   // Initialize entropy
   spline_ntrpy_1 = gsl_spline_alloc (gsl_interp_linear, ne_1);
   gsl_spline_init (spline_ntrpy_1, eps_1, ent_1, ne_1);

   spline_ntrpy_2 = gsl_spline_alloc (gsl_interp_linear, ne_2);
   gsl_spline_init (spline_ntrpy_2, eps_2, ent_2, ne_2);

   spline_ntrpy_3 = gsl_spline_alloc (gsl_interp_linear, ne_3);
   gsl_spline_init (spline_ntrpy_3, eps_3, ent_3, ne_3);

   spline_ntrpy_4 = gsl_spline_alloc (gsl_interp_linear, ne_4);
   gsl_spline_init (spline_ntrpy_4, eps_4, ent_4, ne_4);

   spline_ntrpy_5 = gsl_spline_alloc (gsl_interp_linear, ne_5);
   gsl_spline_init (spline_ntrpy_5, eps_5, ent_5, ne_5);

   spline_ntrpy_6 = gsl_spline_alloc (gsl_interp_linear, ne_6);
   gsl_spline_init (spline_ntrpy_6, eps_6, ent_6, ne_6);

   spline_ntrpy_7 = gsl_spline_alloc (gsl_interp_linear, ne_7);
   gsl_spline_init (spline_ntrpy_7, eps_7, ent_7, ne_7);


   // Initialize pressure
   spline_prs_1 = gsl_spline_alloc (gsl_interp_linear, ne_1);
   gsl_spline_init (spline_prs_1, eps_1, prs_1, ne_1);

   spline_prs_2 = gsl_spline_alloc (gsl_interp_linear, ne_2);
   gsl_spline_init (spline_prs_2, eps_2, prs_2, ne_2);

   spline_prs_3 = gsl_spline_alloc (gsl_interp_linear, ne_3);
   gsl_spline_init (spline_prs_3, eps_3, prs_3, ne_3);

   spline_prs_4 = gsl_spline_alloc (gsl_interp_linear, ne_4);
   gsl_spline_init (spline_prs_4, eps_4, prs_4, ne_4);

   spline_prs_5 = gsl_spline_alloc (gsl_interp_linear, ne_5);
   gsl_spline_init (spline_prs_5, eps_5, prs_5, ne_5);

   spline_prs_6 = gsl_spline_alloc (gsl_interp_linear, ne_6);
   gsl_spline_init (spline_prs_6, eps_6, prs_6, ne_6);

   spline_prs_7 = gsl_spline_alloc (gsl_interp_linear, ne_7);
   gsl_spline_init (spline_prs_7, eps_7, prs_7, ne_7);


   // Teosdens initialize
   spline_teos_e = gsl_spline_alloc (gsl_interp_linear, 750);
   gsl_spline_init (spline_teos_e, teos_t, teos_e, 750);

   spline_teos_p = gsl_spline_alloc (gsl_interp_linear, 750);
   gsl_spline_init (spline_teos_p, teos_t, teos_p, 750);

   spline_teos_s = gsl_spline_alloc (gsl_interp_linear, 750);
   gsl_spline_init (spline_teos_s, teos_t, teos_s, 750);


}


double EoS1::pressure( double eg,double _nb, double _nq, double _ns) //input- energy densty, number densities and output-pressure
{
        if( eg < e0_1  )             { return eg/3.0;}
        if (eg >= e0_1 && eg < e0_2) { return gsl_spline_eval (spline_prs_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { return gsl_spline_eval (spline_prs_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { return gsl_spline_eval (spline_prs_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { return gsl_spline_eval (spline_prs_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { return gsl_spline_eval (spline_prs_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { return gsl_spline_eval (spline_prs_6, eg, acc);}
        if (eg >= e0_7             ) { return gsl_spline_eval (spline_prs_7, eg, acc);}
        else { return eg/3.0; }

}

double EoS1::temperature( double eg,double _nb, double _nq, double _ns) //input- energy densty, number densities and output-temperature
{

        if (eg < e0_1              ) { return gsl_spline_eval(spline_temp_1, e0_1, acc) / e0_1 * (eg); }
        if (eg >= e0_1 && eg < e0_2) { return gsl_spline_eval (spline_temp_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { return gsl_spline_eval (spline_temp_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { return gsl_spline_eval (spline_temp_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { return gsl_spline_eval (spline_temp_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { return gsl_spline_eval (spline_temp_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { return gsl_spline_eval (spline_temp_6, eg, acc);}
        else                         { return gsl_spline_eval (spline_temp_7, eg, acc);}

 
}


double EoS1::entropy( double eg,double _nb, double _nq, double _ns)  //input- energy densty, number densities and output-entropy
{
        if (eg < e0_1              ) { return  gsl_spline_eval(spline_ntrpy_1, e0_1, acc) / e0_1 * (eg);}
        if (eg >= e0_1 && eg < e0_2) { return gsl_spline_eval (spline_ntrpy_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { return gsl_spline_eval (spline_ntrpy_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { return gsl_spline_eval (spline_ntrpy_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { return gsl_spline_eval (spline_ntrpy_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { return gsl_spline_eval (spline_ntrpy_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { return gsl_spline_eval (spline_ntrpy_6, eg, acc);}
        else                         { return gsl_spline_eval (spline_ntrpy_7, eg, acc);}
 
}


double EoS1::cs2( double eg,double _nb, double _nq, double _ns)//input- energy densty, number densities and output-c_{s}^{2}
{
        if (eg < e0_1              ) { return  gsl_spline_eval_deriv(spline_prs_1, e0_1, acc) / e0_1 * (eg);}
        if (eg >= e0_1 && eg < e0_2) { return gsl_spline_eval_deriv (spline_prs_1, eg, acc);}
        if (eg >= e0_2 && eg < e0_3) { return gsl_spline_eval_deriv (spline_prs_2, eg, acc);}
        if (eg >= e0_3 && eg < e0_4) { return gsl_spline_eval_deriv (spline_prs_3, eg, acc);}
        if (eg >= e0_4 && eg < e0_5) { return gsl_spline_eval_deriv (spline_prs_4, eg, acc);}
        if (eg >= e0_5 && eg < e0_6) { return gsl_spline_eval_deriv (spline_prs_5, eg, acc);}
        if (eg >= e0_6 && eg < e0_7) { return gsl_spline_eval_deriv (spline_prs_6, eg, acc);}
        else                         { return gsl_spline_eval_deriv (spline_prs_7, eg, acc);}

}


double EoS1::temp_2_eps( double T,double _nb, double _nq, double _ns)  
{

        double yi=1E-15;
        yi = gsl_spline_eval (spline_teos_e, T, acc);
        return yi;
 
}



//! This function returns local energy density [1/fm^4] from
//! a given entropy density [1/fm^3] and rhob [1/fm^3]
//! using binary search
// [Taken from MUSIC // Copyright 2018 @ Chun Shen ]
double EoS1::entr_2_eps( double s,double _nb, double _nq, double _ns)  
{
    double eps_lower = 1e-15;
    double eps_upper = e0_7;
    double eps_mid   = (eps_upper + eps_lower)/2.;
    double s_lower   = entropy(eps_lower, 0,0,0);
    double s_upper   = entropy(eps_upper, 0,0,0);
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
        s_mid = entropy(eps_mid, 0,0,0);
        if (s < s_mid)
            eps_upper = eps_mid;
        else 
            eps_lower = eps_mid;
        eps_mid = (eps_upper + eps_lower)/2.;
        iter++;
    }
    if (iter == ntol) {
        cout << "get_s2e_finite_rhob:: max iteration reached, "
             << "s = " << s  << endl;
        cout << "s_upper = " << entropy(eps_upper, 0,0,0)
             << " , s_lower = " << entropy(eps_lower, 0,0,0) << endl;
        cout << "eps_upper = " << eps_upper
             << " , eps_lower = " << eps_lower
             << ", diff = " << (eps_upper - eps_lower) << endl;
        exit(1);
    }
    return (eps_mid);
}





void EoS1::check_eos()
{
  std::ofstream File_eos1;
  File_eos1.open("plot_to_check_eos.dat");
  for(double ee =0.0001; ee<= 1000; ee=ee+0.03)
  {
   File_eos1<<temperature(ee,0,0,0)<<"\t"<<ee/pow(temperature(ee,0,0,0),4.0)<<endl;
  }
File_eos1.close();
}








