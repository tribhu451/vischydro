import matplotlib.pyplot as plt
import sys


###################  MUSIC  #######################
X1, Y1  = [],[]
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  if values[0]>0.999 and values[0]<1.000 : # when y = 0 
   X1.append(values[1])                       # put y
   Y1.append(values[2])                 # put epsilon in GeV/fm^3


## xlim ylim
fig = plt.figure()
ax = fig.gca()
#plt.xlim(-0.5,0.5)
#plt.ylim(17.8,18.2)



#plotting ... 
#plt.plot(X1, Y1, ls = '--', label = r'cornelius T_{f} = 150 MeV', c = 'red')
plt.scatter(X1, Y1, marker = 'o', s =10, label = r'cornelius', c = 'red')




#plot label ...
plt.xlabel(r'$y$(fm) [$\eta$=0, x=0]')
plt.ylabel(r'$\tau$')
plt.yscale('linear')
plt.legend()



#some texts in plot
textstr1 = r' (optical glauber IC, $\eta$/s = 0.0)'
textstr3 = r'2+1D, s95pv1 EoS, $\tau_{0}=0.4 fm$'
textstr2 = r'circles ($\eta$/s = 0.04)'
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr1, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
#ax.text(0.05, 0.90, textstr2, transform=ax.transAxes, fontsize=8,
#        verticalalignment='top', bbox=props)
ax.text(0.05, 0.90, textstr3, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)



ax.tick_params(bottom=False, top=True, left=True, right=True)
#ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True)


#plt.show()
plt.savefig('hypersurface2.pdf', format = 'pdf')
                
    

    

