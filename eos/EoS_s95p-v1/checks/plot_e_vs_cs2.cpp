#include<iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

using std::cin;
using std::cout;
using std::endl;
using std::ios;


int
main (void)
{
  int i;
  double xi, yi;

  std::fstream File0;
  File0.open("initial points", ios::out);

  
  std::fstream File1;
  File1.open("spline interpolation", ios::out);

  std::vector<double> x;
  std::vector<double> y;

  for (i = 0; i < 10; i++)
    {
      x.push_back(i + 0.5 * sin (i));
      y.push_back(i + cos (i * i));
      File0<<x[i]<<" \t "<<y[i]<<endl;
    }

   cout<<x.size()<<endl; 

    double x1[x.size()];
    double y1[y.size()];

    for(int i=0; i<x.size(); i++){x1[i]=x[i];}
    for(int i=0; i<y.size(); i++){y1[i]=y[i];}

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, x.size());

    gsl_spline_init (spline, x1, y1, x.size());

    for (xi = x[0]; xi < x[9]; xi += 0.2)
      {
        yi = gsl_spline_eval (spline, xi, acc);
        File1<<xi<<" \t "<<yi<<endl;
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  
  return 0;
}
