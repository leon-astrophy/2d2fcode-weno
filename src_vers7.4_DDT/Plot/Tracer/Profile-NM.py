#import required package#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

#loadeddata#
parameter=np.genfromtxt('../../Outfile/Star_WENO_Parameter.dat')
#load parameter#
factor = 1.01
rhoc1 = parameter[0,1]
rhoc2 = parameter[1,1]
rhoa1 = parameter[2,1]
rhoa2 = parameter[3,1]
fileno = int(parameter[8,1])

for i in range (0,fileno+1):
   essential = np.loadtxt('../../Tracer/Star_WENO_PPT_' + str(i) + '.dat', max_rows=1)
   time = essential[1]
   #load data#
   ppt = np.genfromtxt('../../Tracer/Star_WENO_PPT_' + str(i) + '.dat', skip_header=2)
   # Generate Contour Plot #
   plt.scatter(ppt[:,3], ppt[:,4], s=5, c="r")
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.grid('True')
   plt.title('Tracer Distribution For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('PPT-Tracer'+str(i))
   plt.clf()

   # Generate Contour Plot #
   plt.plot(np.sqrt(ppt[:,3]**2 + ppt[:,4]**2), ppt[:,1], linestyle="--",label="Density")
   plt.xlabel('Distance')
   plt.ylabel('Density')
   plt.grid('True')
   plt.title('Tracer Density Profile For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.savefig('PPT-Tracer-Density'+str(i))
   plt.clf()

   # Generate Contour Plot #
   plt.plot(np.sqrt(ppt[:, 3] ** 2 + ppt[:, 4] ** 2), ppt[:, 2], linestyle="--", label="Temperature")
   plt.xlabel('Distance')
   plt.ylabel('Temperature')
   plt.grid('True')
   plt.title('Tracer Temperature Profile For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.savefig('PPT-Tracer-Temperature'+str(i))
   plt.clf()
