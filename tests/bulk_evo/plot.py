import matplotlib.pyplot as plt
import sys


x_axis_var = 0   #this coloumn will be x-axis
y_axis_var = 1   #this coloumn will be y-axis

  
X1, Y1 = [], []
for line in open(sys.argv[1], 'r'):
  values = [float(s) for s in line.split()]
  X1.append(values[x_axis_var])
  Y1.append(values[y_axis_var])
  
X2, Y2 = [], []
for line in open(sys.argv[2], 'r'):
  values = [float(s) for s in line.split()]
  X2.append(values[x_axis_var])
  Y2.append(values[y_axis_var])

plt.plot(X1, Y1, label = 'analytical')
plt.plot(X2, Y2, label = 'numerical', ls ='--')
plt.legend()
plt.yscale('linear')
#plt.show()
plt.savefig('bulk_evo.eps', format = 'eps')
                
    

    

