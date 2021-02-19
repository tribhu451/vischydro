#import matplotlib.pyplot as plt
import sys

for ix in range(8,18,1):
 plt = 'plt'+str(ix)
 import matplotlib.pyplot as plt
 col1 = int(0); # which coloumn for x-axis
 col2 = int(ix); # which coloumn (y-axis)

 if col2 == 8 :
  _ylabel = r'$\pi^{\tau \tau}$'
  
 if col2 == 9 :
  _ylabel = r'$\pi^{\tau x}$'

 if col2 == 10 :
  _ylabel = r'$\pi^{\tau y}$'

 if col2 == 11 :
  _ylabel = r'$\tau \pi^{\tau \eta}$'

 if col2 == 12 :
  _ylabel = r'$\pi^{x x}$'

 if col2 == 13 :
  _ylabel = r'$\pi^{y x}$'

 if col2 == 14 :
  _ylabel = r'$\tau \pi^{\eta x}$'

 if col2 == 15 :
  _ylabel = r'$\pi^{y y}$'

 if col2 == 16 :
  _ylabel = r'$\tau \pi^{\eta y}$'

 if col2 == 17 :
  _ylabel = r'$\tau^2 \pi^{\eta \eta}$'


 x1, y1 = [], []
 for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  x1.append(values[col1]) 
  y1.append(values[col2]) 

 x2, y2 = [], []
 for line in open(sys.argv[2], 'r'):
  values = [float(s) for s in line.split()]
  x2.append(values[col1]) 
  y2.append(values[col2]) 

 x3, y3 = [], []
 for line in open(sys.argv[3], 'r'):
  values = [float(s) for s in line.split()]
  x3.append(values[col1]) 
  y3.append(values[col2]) 

 x4, y4 = [], []
 for line in open(sys.argv[4], 'r'):
  values = [float(s) for s in line.split()]
  x4.append(values[col1]) 
  y4.append(values[col2]) 
  
 x5, y5 = [], []
 for line in open(sys.argv[5], 'r'):
  values = [float(s) for s in line.split()]
  x5.append(values[col1]) 
  y5.append(values[col2]) 


 fig = plt.figure()
 ax = fig.gca()
 #plt.xlim(0.4,9.5)
 #plt.ylim(0.0,0.155)
 plt.plot(x1, y1, ls = '--', label = r'$\tau$ = 0.5')
 plt.plot(x2,y2, ls = '--', label = r'$\tau$ = 1.0')
 plt.plot(x3,y3, ls = '--', label = r'$\tau$ = 2.0')
 plt.plot(x4,y4, ls = '--', label = r'$\tau$ = 3.0')
 plt.plot(x5,y5, ls = '--', label = r'$\tau$ = 4.0')
 plt.legend()
 plt.xlabel(r'$x$')
 plt.ylabel(_ylabel)
 plt.yscale('linear')
 
 ax.tick_params(bottom=False, top=True, left=True, right=True)
 #ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=True)
 
 aname = str(col2) + '.pdf'
 #plt.show()
 plt.savefig(aname, format = 'pdf')
                
 del plt  

    

