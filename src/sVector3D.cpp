/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

//#include "THGlobal.h"
#include <iostream>
//#include "nag.h"
#include "sVector3D.h"


// STATIC VARIABLES

double sVector3D::mXin;
double sVector3D::mXi0;
double sVector3D::mXip;
double sVector3D::mX0;
double sVector3D::mXd;

inline double sign(double x){ return 0;}//if (x==0) return 0;return x>0?1:-1;}

// OPERATORS

float& sVector3D::operator()(unsigned int i, unsigned int j, unsigned int k)
{
  return mVec[i][j][k];
}

// CLASS FUNCTIONS

sVector3D::sVector3D()
: mVecName(""), mVec(0),
mXmin(0.0), mXmax(1.0), mXpts(1),
mYmin(0.0), mYmax(1.0), mYpts(1),
mZmin(0.0), mZmax(1.0), mZpts(1)
{
  mVec		= new float** [1];
  mVec[0]	= new float*  [1];
  mVec[0][0]	= new float   [1];
  mVec[0][0][0]	= 0.0;
}

sVector3D::sVector3D(const char* aName,
		   double aXmin, double aXmax, int aXpts,
		   double aYmin, double aYmax, int aYpts,
		   double aZmin, double aZmax, int aZpts)
: mVecName(aName),
mXmin(aXmin), mXmax(aXmax), mXpts(aXpts),
mYmin(aYmin), mYmax(aYmax), mYpts(aYpts),
mZmin(aZmin), mZmax(aZmax), mZpts(aZpts)
{
  if(mXpts <= 0) mXpts = 1;
  if(mYpts <= 0) mYpts = 1;
  if(mZpts <= 0) mZpts = 1;
  
  mDi = (mXpts - 1) / (mXmax - mXmin);
  mDj = (mYpts - 1) / (mYmax - mYmin);
  mDk = (mZpts - 1) / (mZmax - mZmin);

  mVec = new float** [mXpts];
  for (int i=0; i<mXpts; i++) {
    mVec[i] = new float* [mYpts];
    for (int j=0; j<mYpts; j++) {
      mVec[i][j] = new float [mZpts];
      for (int k=0; k<mZpts; k++)
	mVec[i][j][k] = 0.0;
    }
  }
}

sVector3D::sVector3D(const sVector3D& aVector)
{
  mXmin = aVector.mXmin; mXmax = aVector.mXmax; mXpts = aVector.mXpts; mDi = aVector.mDi;
  mYmin = aVector.mYmin; mYmax = aVector.mYmax; mYpts = aVector.mYpts; mDj = aVector.mDj;
  mZmin = aVector.mZmin; mZmax = aVector.mZmax; mZpts = aVector.mZpts; mDk = aVector.mDk;

  mVecName = "";
  mVec = new float** [mXpts];
  for (int i=0; i<mXpts; i++) {
    mVec[i] = new float* [mYpts];
    for (int j=0; j<mYpts; j++) {
      mVec[i][j] = new float [mZpts];
      for (int k=0; k<mZpts; k++)
	mVec[i][j][k] = 0.0;
    }
  }
}

sVector3D::~sVector3D()
{
  if(mVec) {
    for (int i=0; i<mXpts; i++) {
      for (int j=0; j<mYpts; j++)
        delete[] mVec[i][j];
      delete[] mVec[i];
    }
    delete[] mVec;
  }
}

std::string	sVector3D::GetName() const	{ return mVecName; }
double		sVector3D::GetXMin() const	{ return mXmin; }
double		sVector3D::GetXMax() const	{ return mXmax; }
int		sVector3D::GetXPts() const	{ return mXpts; }
double		sVector3D::GetYMin() const	{ return mYmin; }
double		sVector3D::GetYMax() const	{ return mYmax; }
int		sVector3D::GetYPts() const	{ return mYpts; }
double		sVector3D::GetZMin() const	{ return mZmin; }
double		sVector3D::GetZMax() const	{ return mZmax; }
int		sVector3D::GetZPts() const	{ return mZpts; }

double sVector3D::Interpolate(double aX, double aY, double aZ)
{
  if(mZpts>1)
    return Interpolate3D(aX, aY, aZ);
  else if (mYpts>1)
    return Interpolate2D(aX, aY);
  else if (mXpts>1)
    return Interpolate1D(aX);
  else {
    return mVec[0][0][0];
  }
}

sVector3D* sVector3D::DerivativeX(const char* aName)
{
  int ti;
  sVector3D* tVec;

  mXd = (mXmax - mXmin) / (mXpts - 1);
  tVec = new sVector3D(*this);
  tVec->mVecName = aName;
  for (int i=0; i<mXpts; i++) {
    ti = InitDerivative(i, mXmin, mXmax, mXpts);
    for (int j=0; j<mYpts; j++)
      for (int k=0; k<mZpts; k++)
	(*tVec)(i,j,k) = Derivative(mVec[ti-1][j][k], mVec[ti][j][k], mVec[ti+1][j][k]);
  }
  return tVec;
}

sVector3D* sVector3D::DerivativeY(const char* aName)
{
  int tj;
  sVector3D* tVec;

  mXd = (mYmax - mYmin) / (mYpts - 1);
  tVec = new sVector3D(*this);
  tVec->mVecName = aName;
  for (int i=0; i<mXpts; i++)
    for (int j=0; j<mYpts; j++) {
      tj = InitDerivative(j, mYmin, mYmax, mYpts);
      for (int k=0; k<mZpts; k++)
	(*tVec)(i,j,k) = Derivative(mVec[i][tj-1][k], mVec[i][tj][k], mVec[i][tj+1][k]);
    }
  return tVec;
}

sVector3D* sVector3D::DerivativeZ(const char* aName)
{
  int tk;
  sVector3D* tVec;

  mXd = (mZmax - mZmin) / (mZpts - 1);
  tVec = new sVector3D(*this);
  tVec->mVecName = aName;
  for (int i=0; i<mXpts; i++)
    for (int j=0; j<mYpts; j++)
      for (int k=0; k<mZpts; k++) {
	tk = InitDerivative(k, mZmin, mZmax, mZpts);
	(*tVec)(i,j,k) = Derivative(mVec[i][j][tk-1], mVec[i][j][tk], mVec[i][j][tk+1]);
      }
  return tVec;
}

double sVector3D::Interpolate1D(double aX)
{
  int    i, j, k;
  double ti;

  ti = (aX - mXmin) * mDi;	i = (int) ti;	if(i+1 > mXpts-1) i--;	ti -= i;
  j  = 0;
  k  = 0;
  return
    mVec[i  ][j  ][k  ] * (1-ti) + mVec[i+1][j  ][k  ] * ti;
}

double sVector3D::Interpolate2D(double aX, double aY)
{
  int    i, j, k;
  double ti,tj;

  ti = (aX - mXmin) * mDi;	i = (int) ti;	if(i+1 > mXpts-1) i--;	ti -= i;
  tj = (aY - mYmin) * mDj;	j = (int) tj;	if(j+1 > mYpts-1) j--;	tj -= j;
  k  = 0;
  return
    (mVec[i  ][j  ][k  ] * (1-ti) + mVec[i+1][j  ][k  ] * ti) * (1-tj) +
    (mVec[i  ][j+1][k  ] * (1-ti) + mVec[i+1][j+1][k  ] * ti) *    tj;
}

double sVector3D::Interpolate3D(double aX, double aY, double aZ)
{
  int    i, j, k;
  double ti,tj,tk;

  ti = (aX - mXmin) * mDi;	i = (int) ti;	if(i+1 > mXpts-1) i--;	ti -= i;
  tj = (aY - mYmin) * mDj;	j = (int) tj;	if(j+1 > mYpts-1) j--;	tj -= j;
  tk = (aZ - mZmin) * mDk;	k = (int) tk;	if(k+1 > mZpts-1) k--;	tk -= k;
  return
    (
      (mVec[i  ][j  ][k  ] * (1-ti) + mVec[i+1][j  ][k  ] * ti) * (1-tj) +
      (mVec[i  ][j+1][k  ] * (1-ti) + mVec[i+1][j+1][k  ] * ti) *    tj
    ) * (1-tk) + (
      (mVec[i  ][j  ][k+1] * (1-ti) + mVec[i+1][j  ][k+1] * ti) * (1-tj) +
      (mVec[i  ][j+1][k+1] * (1-ti) + mVec[i+1][j+1][k+1] * ti) *    tj
    ) *    tk;
}

inline int sVector3D::InitDerivative(int aIdx, double aAMin, double aAMax, int aAPts)
{
  if(aIdx == 0) {
    mXin = aAMin;
    mXi0 = aAMin + mXd;
    mXip = aAMin + mXd * 2.0;
    mX0  = aAMin;
    return aIdx + 1;
  } else if (aIdx == aAPts - 1) {
    mXin = aAMax - mXd * 2.0;
    mXi0 = aAMax - mXd;
    mXip = aAMax;
    mX0  = aAMax;
    return aIdx - 1;
  } else {
    mXin = aAMin + mXd * (aIdx - 1);
    mXi0 = aAMin + mXd *  aIdx;
    mXip = aAMin + mXd * (aIdx + 1);
    mX0  = mXi0;
    return aIdx;
  }
}

inline double sVector3D::Derivative(double aFin, double aFi0, double aFip)
{
  return (
    mXin * (mXin - 2.0 * mX0) * (aFip - aFi0) +
    mXi0 * (mXi0 - 2.0 * mX0) * (aFin - aFip) +
    mXip * (mXip - 2.0 * mX0) * (aFi0 - aFin)
  ) / (2.0 * mXd * mXd * mXd);
}


void sVector3D::put(double x,int i,int j,int k)
{
  mVec[i][j][k]=x;
}

// added by Piotr Bozek



double sVector3D::derx(int i, int j, int k){
  //      if (i==0||j==0||k==0||i==(mXpts-1)||j==(mYpts-1)||k==(mZpts-1))
  //       {return 0.;}
  if (i==0){
    double d2=0.5*mDi*(4.*mVec[i+1][j  ][k  ]-3.*mVec[i][j][k]
		       -mVec[i+2][j  ][k  ]);
    return d2; 
  }
  else
    if (i==(mXpts-1)){
      double d2=0.5*mDi*(-4.*mVec[i-1][j  ][k  ]+3.*mVec[i][j][k]
			+mVec[i-2][j  ][k  ]);
      return d2;
    }
    else{
      double d3=0.5*mDi*(mVec[i+1][j  ][k  ]-mVec[i-1][j  ][k  ]);
      return d3;
    }
}


double sVector3D::dery(int i, int j, int k){
  //    if (i==0||j==0||k==0||i==(mXpts-1)||j==(mYpts-1)||k==(mZpts-1))
  //     {return 0.;}
  if (j==0){
    double d2=0.5*mDj*(4.*mVec[i  ][j+1][k  ]-3.*mVec[i][j][k]
		      -mVec[i  ][j+2][k  ]);
    return d2;
  }
  else
    if (j==(mYpts-1)){
      double d2=0.5*mDj*(-4.*mVec[i  ][j-1][k  ]+3.*mVec[i][j][k]
			+mVec[i  ][j-2][k  ]);
      return d2;
    }
    else{
      double d3=0.5*mDj*(mVec[i  ][j+1][k  ]-mVec[i  ][j-1][k  ]);
      return d3;
    }
}


double sVector3D::dere(int i, int j, int k){
  //   if (i==0||j==0||k==0||i==(mXpts-1)||j==(mYpts-1)||k==(mZpts-1))
  //     {return 0.;}
  if (k==0){
    double d2=0.5*mDk*(4.*mVec[i  ][j  ][k+1]-3.*mVec[i][j][k]
		      -mVec[i  ][j  ][k+2]);
    return d2;
  }
  else
    if (k==(mZpts-1)){
      double d2=0.5*mDk*(-4.*mVec[i  ][j  ][k-1]+3.*mVec[i][j][k]
			+mVec[i  ][j  ][k-2]);
      return d2;
    }
    else{
      double d3=0.5*mDk*(mVec[i  ][j  ][k+1]-mVec[i  ][j  ][k-1]);
      return d3;
    }
}

/*
double sVector3D::derx(int i, int j, int k){
  //      if (i==0||j==0||k==0||i==(mXpts-1)||j==(mYpts-1)||k==(mZpts-1))
  //       {return 0.;}
  if (i==0){
    double d1=(mVec[i+1][j  ][k  ]-mVec[i][j][k])*mDi*par::thetaflux;
    double d2=0.5*mDi*(4.*mVec[i+1][j  ][k  ]-3.*mVec[i][j][k]
		       -mVec[i+2][j  ][k  ]);
    return minmod2(d1,d2);
  }
  else
    if (i==(mXpts-1)){
      double d1=(mVec[i][j][k]-mVec[i-1][j  ][k  ])*mDi*par::thetaflux;
      double d2=0.5*mDi*(-4.*mVec[i-1][j  ][k  ]+3.*mVec[i][j][k]
			+mVec[i-2][j  ][k  ]);
      return minmod2(d1,d2);
    }
    else{
      double d1=(mVec[i][j][k]-mVec[i-1][j  ][k  ])*mDi*par::thetaflux;
      double d2=(mVec[i+1][j  ][k  ]-mVec[i][j][k])*mDi*par::thetaflux;
      double d3=0.5*mDi*(mVec[i+1][j  ][k  ]-mVec[i-1][j  ][k  ]);
      return minmod3(d1,d2,d3);
    }
}


double sVector3D::dery(int i, int j, int k){
  //    if (i==0||j==0||k==0||i==(mXpts-1)||j==(mYpts-1)||k==(mZpts-1))
  //     {return 0.;}
  if (j==0){
    double d1=(mVec[i  ][j+1][k  ]-mVec[i][j][k])*mDj*par::thetaflux;
    double d2=0.5*mDj*(4.*mVec[i  ][j+1][k  ]-3.*mVec[i][j][k]
		      -mVec[i  ][j+2][k  ]);
    return minmod2(d1,d2);
  }
  else
    if (j==(mYpts-1)){
      double d1=(mVec[i][j][k]-mVec[i  ][j-1][k  ])*mDj*par::thetaflux;
      double d2=0.5*mDj*(-4.*mVec[i  ][j-1][k  ]+3.*mVec[i][j][k]
			+mVec[i  ][j-2][k  ]);
      return minmod2(d1,d2);
    }
    else{
     double d1=(mVec[i][j][k]-mVec[i  ][j-1][k  ])*mDj*par::thetaflux;
      double d2=(mVec[i  ][j+1][k  ]-mVec[i][j][k])*mDj*par::thetaflux;
      double d3=0.5*mDj*(mVec[i  ][j+1][k  ]-mVec[i  ][j-1][k  ]);
      return minmod3(d1,d2,d3);
    }
}


double sVector3D::dere(int i, int j, int k){
  //   if (i==0||j==0||k==0||i==(mXpts-1)||j==(mYpts-1)||k==(mZpts-1))
  //     {return 0.;}
  if (k==0){
    double d1=(mVec[i  ][j  ][k+1]-mVec[i][j][k])*mDk*par::thetaflux;
    double d2=0.5*mDk*(4.*mVec[i  ][j  ][k+1]-3.*mVec[i][j][k]
		      -mVec[i  ][j  ][k+2]);
    return minmod2(d1,d2);
  }
  else
    if (k==(mZpts-1)){
      double d1=(mVec[i][j][k]-mVec[i  ][j  ][k-1])*mDk*par::thetaflux;
      double d2=0.5*mDk*(-4.*mVec[i  ][j  ][k-1]+3.*mVec[i][j][k]
			+mVec[i  ][j  ][k-2]);
      return minmod2(d1,d2);
    }
    else{
      double d1=(mVec[i][j][k]-mVec[i  ][j  ][k-1])*mDk*par::thetaflux;
      double d2=(mVec[i  ][j  ][k+1]-mVec[i][j][k])*mDk*par::thetaflux;
      double d3=0.5*mDk*(mVec[i  ][j  ][k+1]-mVec[i  ][j  ][k-1]);
      return minmod3(d1,d2,d3);
    }
}
*/

double sVector3D::minmod2(double x, double y){
  return y;
  //     return 0;
  if ((x>0)&&(y>0)){
    return    (x < y) ? x : y;
  }
  else
    if ((x<0)&&(y<0)){
      return    (x < y) ?  y : x;
    }
    else{
      return 0.;
    }
}

double sVector3D::minmod3(double x, double y, double z){
  //   return 0;
  return z;
  if ((x>0)&&(y>0)&&(z>0)){
        x=(x < y) ? x : y;
	return (x < z) ? x : z;
  }
  else
    if ((x<0)&&(y<0)&&(z<0)){
      x=(x < y) ?  y : x;
      return (x < z) ?  z : x;
    }
    else{
      return 0.;
    }
}
