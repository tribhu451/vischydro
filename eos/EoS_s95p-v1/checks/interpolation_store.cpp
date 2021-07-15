#include<iostream>
#include <cmath>
#include <fstream>
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
  double xi, yi, x[10], y[10];

  std::fstream File0;
  File0.open("initial points", ios::out);

  
  std::fstream File1;
  File1.open("spline interpolation", ios::out);


  for (i = 0; i < 10; i++)
    {
      x[i] = i + 0.5 * sin (i);
      y[i] = i + cos (i * i);
      File0<<x[i]<<" \t "<<y[i]<<endl;
    }

 
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, 10);

    gsl_spline_init (spline, x, y, 10);

    for (xi = x[0]; xi < x[9]; xi += 0.2)
      {
        yi = gsl_spline_eval (spline, xi, acc);
        File1<<xi<<" \t "<<yi<<endl;
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  
  return 0;
}
