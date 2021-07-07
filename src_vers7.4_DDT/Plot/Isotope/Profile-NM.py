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
fileno = int(parameter[7,1])

for i in range (0,fileno+1):
   essential=np.loadtxt('../../Outfile/Isotope/Star_WENO_XC12_'+str(i)+'.dat',max_rows=3)
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
   ye = np.loadtxt('../../Outfile/Hydro/Star_WENO_Ye_NM_'+str(i)+'.dat', skiprows=3)
   c12 = np.loadtxt('../../Outfile/Isotope/Star_WENO_XC12_'+str(i)+'.dat', skiprows=3)
   o16 = np.loadtxt('../../Outfile/Isotope/Star_WENO_XO16_' + str(i) + '.dat', skiprows=3)

   # Generate Contour Plot #
   plt.clf()
   contoursf = plt.contourf(R, Z, ye, 100, cmap='jet')
   #contours = plt.contour(R, Z, ye, 10, cmap='binary')
   fmt = ticker.LogFormatterMathtext()
   #plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('Ye Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('YeContours'+str(i))
   #plt.show()

   # Generate Contour Plot #
   plt.clf()
   contoursf = plt.contourf(R, Z, c12, 100, cmap='jet')
   #contours = plt.contour(R, Z, c12, 10, cmap='binary')
   fmt = ticker.LogFormatterMathtext()
   #plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('C12 Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('C12Contours'+str(i))
   #plt.show()

   # Generate Contour Plot #
   plt.clf()
   contoursf = plt.contourf(R, Z, o16, 100, cmap='jet')
   #contours = plt.contour(R, Z, o16, 10, cmap='binary')
   fmt = ticker.LogFormatterMathtext()
   #plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('O16 Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('O16Contours'+str(i))
   #plt.show()

   #generate desnity profile#
   plt.clf()
   plt.plot(gridr,c12[0,:],label='r-axis')
   plt.plot(gridz,c12[:,0],label='z-axis')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('C12')
   plt.title('C12 Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('C12Profile'+str(i))

   #generate desnity profile#
   plt.clf()
   plt.plot(gridr,o16[0,:],label='r-axis')
   plt.plot(gridz,o16[:,0],label='z-axis')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('O16')
   plt.title('O16 Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('O16Profile'+str(i))

   #generate desnity profile#
   plt.clf()
   plt.plot(gridr,ye[0,:],label='r-axis')
   plt.plot(gridz,ye[:,0],label='z-axis')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('Ye')
   plt.title('Ye Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('YeProfile'+str(i))
