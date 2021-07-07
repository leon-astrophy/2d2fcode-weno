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
fileno = int(parameter[4,1])

for i in range (0,fileno+1):
   essential=np.loadtxt('../../Outfile/Hydro/Star_WENO_Density_NM_'+str(i)+'.dat',max_rows=3)
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
   rhonm = np.loadtxt('../../Outfile/Hydro/Star_WENO_Density_NM_' + str(i) + '.dat', skiprows=3)
   velnmr = np.loadtxt('../../Outfile/Hydro/Star_WENO_Velocity_NM_r_' + str(i) + '.dat', skiprows=3)
   velnmz = np.loadtxt('../../Outfile/Hydro/Star_WENO_Velocity_NM_z_' + str(i) + '.dat', skiprows=3)
   # Generate Contour Plot #
   plt.clf()
   levls2 = np.logspace(np.log10(rhoa2), np.log10(rhoc2), 10)
   levls1 = np.logspace(np.log10(rhoa2), np.log10(rhoc2), 100)
   contoursf = plt.contourf(R, Z, rhonm, 100, cmap='jet')
   #contours = plt.contour(R, Z, rhonm, 10, cmap='binary')
   fmt = ticker.LogFormatterMathtext()
   #plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   quiv = plt.quiver(R[::20, ::20], Z[::20, ::20], velnmr[::20, ::20], velnmz[::20, ::20])
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('Denisty Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('NMDensityContours'+str(i))
   #plt.show()

   #load data#
   phinm = np.loadtxt('../../Outfile/Hydro/Star_WENO_Potential_NM_' + str(i) + '.dat', skiprows=3)
   # Generate Contour Plot #
   plt.clf()
   contoursf = plt.contourf(R, Z, phinm, 100, cmap='jet')
   #contours = plt.contour(R, Z, phinm, 10, cmap='binary')
   fmt = ticker.LogFormatterMathtext()
   #plt.clabel(contours, inline=True, fontsize=10, fmt=fmt)
   plt.colorbar(contoursf)
   plt.xlabel('X-Axis')
   plt.ylabel('Z-Axis')
   plt.title('Potentials Contours For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.axis('square')
   plt.savefig('NMPotentialContours' + str(i))
   #plt.show()

   #generate desnity profile#
   plt.clf()
   plt.plot(gridr,velnmr[0,:],label='r-axis, vr')
   plt.plot(gridz,velnmz[:,0],label='z-axis, vz')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('Velocity')
   plt.title('Velocity Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('NMVelocityProfile'+str(i))

   #generate desnity profile#
   plt.clf()
   plt.plot(gridr,np.log10(rhonm[0,:]),label='r-axis')
   plt.plot(gridz,np.log10(rhonm[:,0]),label='z-axis')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('Log10 Density')
   plt.title('Density Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('NMDensityProfile'+str(i))

   #generate potential profile#
   plt.clf()
   plt.plot(gridr,phinm[0,:],label='r-axis')
   plt.plot(gridz,phinm[:,0],label='z-axis')
   plt.legend(loc="lower right")
   plt.xlabel('Distance')
   plt.ylabel('Potentials')
   plt.title('Potentials Profile Along Different Axis', y=1.08)
   plt.grid(True)
   plt.tight_layout()
   plt.savefig('NMPotentialProfile' + str(i))