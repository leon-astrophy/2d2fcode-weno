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
fileno = int(parameter[5,1])

for i in range (0,fileno+1):
   essential=np.loadtxt('../../Outfile/Flame/Star_WENO_scaG_'+str(i)+'.dat',max_rows=3)
   #read in essential data#
   jr = int(essential[1,0])
   jz = int(essential[1,1])
   #time#
   time = essential[0,1]
   #assign grid value to x and z axis#
   dx = essential[2,0]
   gridr = np.zeros(jr)
   for j in range (0, jr):
      gridr[j] = dx * (float(j) - 0.5)
   gridz = np.zeros(jz)
   for j in range (0, jz):
      gridz[j] = dx * (float(j) - 0.5)
   # mesh the grid for x and z#
   R,Z=np.meshgrid(gridr,gridz)

   #load data#
   scag = np.loadtxt('../../Outfile/Flame/Star_WENO_scaG_'+str(i)+'.dat', skiprows=3)
   scag2 = np.loadtxt('../../Outfile/Flame/Star_WENO_scaG2_' + str(i) + '.dat', skiprows=3)
   temp = np.loadtxt('../../Outfile/Hydro/Star_WENO_Temperature_NM_' + str(i) + '.dat', skiprows=3)

   # Generate Contour Plot #
   plt.clf()
   level = [0]
   contoursf = plt.contourf(R, Z, temp, 100, cmap='jet')
   contours = plt.contour(R, Z, scag, levels = level, cmap='binary')
   fmt = ticker.LogFormatterMathtext()
   plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('Temperature Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('FlameTemperatureContours'+str(i))

   # Generate Contour Plot #
   plt.clf()
   level = [0]
   contoursf = plt.contourf(R, Z, temp, 100, cmap='jet')
   contours = plt.contour(R, Z, scag2, levels = level, cmap='Greys')
   fmt = ticker.LogFormatterMathtext()
   plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('Temperature Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('DetonTemperatureContours'+str(i))

   #generate desnity profile#
   plt.clf()
   plt.plot(gridr,temp[0,:],label='r-axis')
   plt.plot(gridz,temp[:,0],label='z-axis')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('Temperature')
   plt.title('Temperature Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('TemperatureProfile'+str(i))