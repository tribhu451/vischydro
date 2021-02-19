import matplotlib.pyplot as plt
import sys


X1, Y1 = [], []
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  X1.append(values[0])
  Y1.append(values[1])


X2, Y2 = [], []
for line in open(sys.argv[2], 'r'):
  values = [float(s) for s in line.split()]
  X2.append(values[0])
  Y2.append(values[2])

X3, Y3 = [], []
for line in open(sys.argv[3], 'r'):
  values = [float(s) for s in line.split()]
  X3.append(values[0])
  Y3.append(values[1])


X4, Y4 = [], []
for line in open(sys.argv[4], 'r'):
  values = [float(s) for s in line.split()]
  X4.append(values[0])
  Y4.append(values[2])

X5, Y5 = [], []
for line in open(sys.argv[5], 'r'):
  values = [float(s) for s in line.split()]
  X5.append(values[0])
  Y5.append(values[1])


X6, Y6 = [], []
for line in open(sys.argv[6], 'r'):
  values = [float(s) for s in line.split()]
  X6.append(values[0])
  Y6.append(values[3])



fig = plt.figure()
ax = fig.gca()
plt.xlim(0.,8.0)
plt.ylim(0.0,0.185)
plt.scatter(X1, Y1, marker = 'o', s =40, label = r'VISHNU 2+1D ($\eta / s = 0.0$) [$\epsilon_p$]')
plt.plot(X2, Y2, ls = '-', label = r'my code ($\eta / s = 0.0$) [$\epsilon_p$] ', linewidth=5)
plt.scatter(X3, Y3, marker = 'o', s = 40, label = r'VISHNU 2+1D ($\eta / s = 0.08$) [$\epsilon_p$]')
plt.plot(X4, Y4, ls = '-', label = r'my code ($\eta / s = 0.08$) [$\epsilon_p$]', linewidth=5)
plt.scatter(X5, Y5, marker = 'o', s =40, label = r'VISHNU 2+1D ($\eta / s = 0.08$) [$\epsilon_{p}^{\prime}$]')
plt.plot(X6, Y6, ls = '-', label = r'my code ($\eta / s = 0.08$) [$\epsilon_{p}^{\prime}$]', linewidth=5)
plt.xlabel(r'$\tau$ (fm)')
plt.ylabel(r'$\frac{ \langle T^{xx} - T^{yy} \rangle } { \langle T^{xx} + T^{yy} \rangle} $')
plt.yscale('linear')
plt.legend()

ax.tick_params(bottom=False, top=True, left=True, right=True)
#ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True)


#plt.show()
plt.savefig('anisotropy.pdf', format = 'pdf')
                
    

    

