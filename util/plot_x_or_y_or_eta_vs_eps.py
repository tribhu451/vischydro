import matplotlib.pyplot as plt
import sys


###################  MUSIC  #######################
X1, Y1  = [],[]
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  if values[0]>12.9 and values[0]<13.1 :             # when time = 1 fm 
     if values[1]>-0.001 and values[1]<0.001 :     # when eta = 0 
       if values[2]>-0.001 and values[2]<0.001 :   # when y = 0 
        X1.append(values[3])                       # put x
        Y1.append(values[4]/5.068)                 # put epsilon in GeV/fm^3


X2, Y2  = [],[]
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  if values[0]>12.9 and values[0]<13.1 :             #when time = 1 fm 
     if values[2]>-0.001 and values[2]<0.001 : # when eta = 0 
       if values[3]>-0.001 and values[3]<0.001 : # when y = 0 
        X2.append(values[1])                        # put x
        Y2.append(values[4]/5.068)                 # put epsilon in GeV/fm^3


X3, Y3  = [],[]
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  if values[0]>12.9 and values[0]<13.1 :             #when time = 1 fm 
     if values[3]>-0.001 and values[3]<0.001 : # when eta = 0 
       if values[1]>-0.001 and values[1]<0.001 : # when y = 0 
        X3.append(values[2])                        # put x
        Y3.append(values[4]/5.068)                 # put epsilon in GeV/fm^3



###################  vischydro  #######################
A1, B1  = [], []
for line in open(sys.argv[2], 'r'):
  values = [float(s) for s in line.split()]
  if values[1]>-0.0001 and values[1]<0.0001 :             #when time = 1 fm 
   A1.append(values[0])                         # put eta
   B1.append(values[4])                         # put epsilon

A2, B2  = [], []
for line in open(sys.argv[2], 'r'):
  values = [float(s) for s in line.split()]
  if values[0]>-0.0001 and values[0]<0.0001 :             #when time = 1 fm 
   A2.append(values[1])                         # put eta
   B2.append(values[4])                         # put epsilon


A3, B3  = [], []
for line in open(sys.argv[3], 'r'):
  values = [float(s) for s in line.split()]
  if values[0]>-0.0001 and values[0]<0.0001 :             #when time = 1 fm 
   A3.append(values[3])                         # put eta
   B3.append(values[4])                         # put epsilon



## xlim ylim
fig = plt.figure()
ax = fig.gca()
plt.xlim(-12.5,12.5)
#plt.ylim(17.8,18.2)



#plotting ... 
plt.plot(X1, Y1, ls = '-', label = r'MUSIC 13.0fm', color = 'orange')
plt.plot(X2, Y2, ls = '-', label = r'MUSIC 13.0fm', color = 'deepskyblue')
plt.plot(X3, Y3, ls = '-', label = r'MUSIC 13.0fm', color = 'lawngreen')
plt.plot(A1, B1, ls = '-', label = r'vischydro 13.0fm', color = 'blue')
plt.plot(A2, B2, ls = '-', label = r'vischydro 13.0fm', color = 'green')
plt.plot(A3, B3, ls = '-', label = r'vischydro 13.0fm', color = 'red')




#plot label ...
plt.xlabel(r'x/ y/ $\eta$')
plt.ylabel(r'$\epsilon$(GeV/fm$^3$)')
plt.yscale('linear')
plt.legend()



#some texts in plot
textstr1 = r' $\eta$/s(T) & $\zeta/s(T)$'
textstr3 = r'3+1D, s95pv1 EoS, $\tau_{0}=1.0 fm$'
textstr2 = 'optical Glauber IC'
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr1, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax.text(0.05, 0.90, textstr2, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
ax.text(0.05, 0.85, textstr3, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)



ax.tick_params(bottom=False, top=True, left=True, right=True)
#ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True)


#plt.show()
plt.savefig('music_vs_vischydro_13_fm.pdf', format = 'pdf')
                
    

    

