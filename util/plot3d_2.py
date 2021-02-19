import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys


x, y, z  = [],[],[]
for line in open(sys.argv[1], 'r'):
 values = [float(s) for s in line.split()]
 x.append(values[1])                       
 y.append(values[2])                       
 z.append(values[0])        

fig = plt.figure(figsize =(16, 9))   
ax = plt.axes(projection ='3d')   
  
# Creating color map 
my_cmap = plt.get_cmap('jet') 
    
# Creating plot 
trisurf = ax.plot_trisurf(x, y, z, 
                         cmap = my_cmap, lw = 2.)
                         #linewidth = 0.5,  
                         #antialiased = True, 
                         #edgecolor = 'grey'
                            
fig.colorbar(trisurf, ax = ax, shrink = 0.5, aspect = 5) 
ax.set_title('hypersurface') 
  
# Adding labels 
ax.set_xlabel('X', fontweight ='bold')  
ax.set_ylabel('Y', fontweight ='bold')  
ax.set_zlabel(r'\tau', fontweight ='bold') 
   
#ax.view_init(20.25, -160.0)
  
# Info
'''
Vary azim
The azimuth is the rotation around the z axis e.g.:
0 means "looking from +x"
90 means "looking from +y"

dist seems to be the distance from the center visible point in data coordinates.


From this we understand that elev is the angle between the eye and the xy plane.

'''


 
ax.azim = 0
ax.dist = 10
ax.elev = 0




# show plot 
#plt.show() 
plt.savefig('3d_optical_b0.png', format = 'png')
