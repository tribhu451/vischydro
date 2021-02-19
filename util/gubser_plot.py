import matplotlib.pyplot as plt
import sys


x_axis_var = 0
y_axis_var = 3

argc = len(sys.argv) 

if argc == 5 :
    
    
    X1, Y1 = [], []
    for line in open(sys.argv[1], 'r'):
        values = [float(s) for s in line.split()]
        X1.append(values[x_axis_var])
        Y1.append(values[y_axis_var])
        
    X2, Y2 = [], []
    for line2 in open(sys.argv[2], 'r'):
        values = [float(s) for s in line2.split()]
        X2.append(values[x_axis_var])
        Y2.append(values[y_axis_var])
            
            
    X3, Y3 = [], []
    for line in open(sys.argv[3], 'r'):
        values = [float(s) for s in line.split()]
        X3.append(values[x_axis_var])
        Y3.append(values[y_axis_var])
                
                
    X4, Y4 = [], []
    for line2 in open(sys.argv[4], 'r'):
        values = [float(s) for s in line2.split()]
        X4.append(values[x_axis_var])
        Y4.append(values[y_axis_var])
                    
        
    plt.plot(X1, Y1, label = 'analytical (3fm)')
    plt.plot(X2, Y2, ls = '--', label =  'numerical (3fm)')
    plt.plot(X3, Y3, label = 'analytical (5fm)')
    plt.plot(X4, Y4, ls = '--', label =  'numerical (5fm)')
    plt.legend()
    plt.yscale('log')


    plt.xlabel('x')
    plt.ylabel('$\epsilon$')
    #plt.show()
    plt.savefig('myimage.eps', format='eps', dpi=1200)

elif argc == 3 :
   
    X1, Y1 = [], []
    for line in open(sys.argv[1], 'r'):
        values = [float(s) for s in line.split()]
        X1.append(values[x_axis_var])
        Y1.append(values[y_axis_var])
        
    X2, Y2 = [], []
    for line2 in open(sys.argv[2], 'r'):
        values = [float(s) for s in line2.split()]
        X2.append(values[x_axis_var])
        Y2.append(values[y_axis_var])

    plt.plot(X1, Y1, label = 'analytical')
    plt.plot(X2, Y2, ls = '--', label =  'numerical')
    plt.legend()
    plt.yscale('linear')
    plt.show()
                
    
else :
    
    print('plz give 8 file names/2 file names after python3 ./gubser_plot.py ...')
