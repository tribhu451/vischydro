# give 6 arguments 
# semi-analytical 1.2fm , numerical 1.2 fm , semi-analytic 1.5fm, numerical 1.5fm 
# semi-analytical 2.0fm , numerical 2.0 fm


# python3 util/gubser_visc_test.py tests/gubser/y\=0_tau\=1.2_SemiAnalytic.dat hydro_output/midslice_xy_dist_at_1.200000_fm.dat tests/gubser/y\=0_tau\=1.5_SemiAnalytic.dat hydro_output/midslice_xy_dist_at_1.500000_fm.dat tests/gubser/y\=0_tau\=2.0_SemiAnalytic.dat hydro_output/midslice_xy_dist_at_2.000000_fm.dat

import matplotlib.pyplot as plt
import sys


py1 = 7;  # semi - analytical column no.
py2 = 13; # from fluid_info.cpp

  
X1, Y1 = [], []
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  X1.append(values[0])
  Y1.append(values[py1])

X2, Y2 = [], []
for line in open(sys.argv[2], 'r'):
  values = [float(s) for s in line.split()]
  X2.append(values[0])
  Y2.append(values[py2])

X3, Y3 = [], []
for line in open(sys.argv[3], 'r'):
  values = [float(s) for s in line.split()]
  X3.append(values[0])
  Y3.append(values[py1])

X4, Y4 = [], []
for line in open(sys.argv[4], 'r'):
  values = [float(s) for s in line.split()]
  X4.append(values[0])
  Y4.append(values[py2])

X5, Y5 = [], []
for line in open(sys.argv[5], 'r'):
  values = [float(s) for s in line.split()]
  X5.append(values[0])
  Y5.append(values[py1])

X6, Y6 = [], []
for line in open(sys.argv[6], 'r'):
  values = [float(s) for s in line.split()]
  X6.append(values[0])
  Y6.append(values[py2])


fig = plt.figure()
ax = fig.gca()
plt.xlim(-4.8,4.8)
#plt.ylim(0.0,0.155)
plt.plot(X1, Y1, ls = '-', label = r'semi-analytical')
plt.plot(X2, Y2, ls = '--', label = r'numerical ($\tau=$1.2fm)')
plt.plot(X3, Y3, ls = '-', label = r'semi-analytical')
plt.plot(X4, Y4, ls = '--', label = r'numerical ($\tau = $1.5fm)')
plt.plot(X5, Y5, ls = '-', label = r'semi-analytical')
plt.plot(X6, Y6, ls = '--', label = r'numerical ($\tau = $2.0fm)')
plt.legend()
plt.xlabel(r'x')
#plt.ylabel(r'$\tau^2 \pi^{\eta \eta}$ (GeV/fm$^3$)')
plt.ylabel(r'$\pi^{xy}$ (GeV/fm$^3$)')
#plt.yscale('linear')

ax.tick_params(bottom=False, top=True, left=True, right=True)
#ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True)


#plt.show()
plt.savefig('pixy.pdf', format = 'pdf')
                
    

    

