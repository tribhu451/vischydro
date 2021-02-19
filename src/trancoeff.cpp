#include "trancoeff.h"

trancoeff::trancoeff(EoS *_eos, idb* _InData)
{
  InData = _InData ;
  etas_flag = InData->etas_flag;
  t_etas_flag = InData->t_etas_flag;
  etaS = InData->etas;
  zetas_flag = InData->zetas_flag;
  eos = _eos;

  if(etas_flag > 0 || zetas_flag > 0 )
    cout<<"[info] viscous hydro running ..."<<endl;
  else
   cout<<"[info] Ideal hydro running ... "<<endl;
  
  if(etas_flag > 0 && t_etas_flag > 0)
    cout<<"[info] eta / s  is temp depedent "<<endl;

  if(etas_flag > 0 && t_etas_flag == 0)
    cout<<"[info] const eta/s = "<<etaS<<endl;

  if (zetas_flag > 0)
    cout<<"[info] bulk on with temp. dependet zeta/s "<<endl;

  if (etas_flag > 0 && zetas_flag == 0)
    cout<<"[info] bulk off zeta/s = 0 "<<endl;

  cout<<"\n"<<endl;
  
}

trancoeff::~trancoeff()
{
}

double trancoeff::EtaS(double e, double temp)
{ 
  if(t_etas_flag > 0)
    {
      double Ttr = 0.18 ;
      double etas;
      if (temp<Ttr)
	etas = 0.681 - 0.0594*temp/Ttr - 0.544*pow(temp/Ttr,2.0) ;
      else
	etas = -0.289 + 0.288*temp/Ttr + 0.0818*pow(temp/Ttr,2.0) ;
      
      return etas;
    }
  else
    {
      return etaS ;
    }
}



double trancoeff::zetaS(double e, double temp)
{
  
  //temperature dependent zeta/s
  // as used in MUSIC [by Gabriel]
  temp /= 0.19733 ;
  double Ttr = 0.18/0.1973;
  double dummy = temp/Ttr ;
  
  double A1 = -13.77;
  double A2 = 27.55 ;
  double A3 = 13.45 ;
  
  double lambda1 = 0.90;
  double lambda2 = 0.25;
  double lambda3 = 0.90;
  double lambda4 = 0.22;
  
  double sigma1 = 0.0250;
  double sigma2 = 0.1300;
  double sigma3 = 0.0025;
  double sigma4 = 0.0220;
  
  
  double zetas = A1*dummy*dummy + A2*dummy - A3 ;
  if(temp < 0.995*Ttr)
    {
      zetas = lambda3*exp( ( dummy -1 ) / sigma3 ) + lambda4*exp( ( dummy -1 ) / sigma4 ) + 0.03000 ;
    } 
  
  if(temp > 1.05*Ttr)
    {
      zetas = lambda1*exp( -( dummy -1 ) / sigma1 ) + lambda2*exp( -( dummy -1 ) / sigma2 ) + 0.001 ;
    } 
  
  return zetas ;
}



void trancoeff::getEta(double e, double T, double &_etaS, double &_zetaS)
{
  _etaS = EtaS(e,T);
  
  if(zetas_flag > 0)
    {_zetaS = zetaS(e,T);}
  else {_zetaS = 0 ;}
}


void trancoeff::getTau(double e, double T, double &_taupi, double &_tauPi)
{
  if (T > 0.)
    {
      _taupi = 5. / 5.068 * EtaS(e,T) / (T + 1e-15); 
    }
  else
    {
      cout<<"temperature was negative : "<<T<<" -> making it positiive. 1E-15 "<<endl;
      T = 1E-15;
      _taupi = 5. / 5.068 * EtaS(e,T) / (T + 1e-15); 
    }
  
  if (T > 0.)
    { _tauPi = 1. / 5.068 * zetaS(e,T) / (14.55 * pow(1./3. - eos->cs2(e,0,0,0),2) * (T+1e-15) );}
  else
    {
      cout<<"temperature was negative : "<<T<<" -> making it positiive. 1E-15"<<endl;
      T = 1E-15;
      _tauPi = 1. / 5.068 * zetaS(e,T) / (14.55 * pow(1./3. - eos->cs2(e,0,0,0),2) * (T+1e-15) );
    }
  
}
