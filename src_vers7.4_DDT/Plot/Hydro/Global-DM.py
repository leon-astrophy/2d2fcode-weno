#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
mass1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Mass_DM_0.dat',comments='"')
energy1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Energy_DM_0.dat',comments='"')
centralrho1=np.loadtxt('../../Outfile/Hydro/Star_WENO_CentralDensity_DM_0.dat',comments='"')

#plot#
plt.plot(centralrho1[:,0], centralrho1[:,1], label='Central Density')
plt.plot(centralrho1[:,0], centralrho1[:,2], label='Maximum Density')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Density (CODE)')
plt.title('Density Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('DMCentralDensity')
plt.clf()

#plot#
plt.plot(mass1[:,0], np.abs((mass1[:,1] - mass1[0,1])/(mass1[0,1])))
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Fractional Change')
plt.title('Mass Conservation')
plt.grid(True)
plt.tight_layout()
plt.savefig('DMMass')
plt.clf()

#plot#
plt.plot(energy1[:,0], energy1[:,1], label='Total DM Energy')
plt.plot(energy1[:,0], energy1[:,2], label='Kinetic Energy')
plt.plot(energy1[:,0], energy1[:,3], label='Internal Energy')
plt.plot(energy1[:,0], energy1[:,4], label='Gravitational Energy')
plt.legend(loc="lower right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy (CODE)')
plt.title('Energy Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('DMEnergy')
plt.clf()

#plot#
plt.plot(energy1[:,0], np.abs((energy1[:,1] - energy1[0,1])/(energy1[0,1])), label='Total DM Energy')
plt.legend(loc="lower right")
plt.xlabel('Time (CODE)')
plt.ylabel('Fractional Change')
plt.title('Energy Conservation')
plt.grid(True)
plt.tight_layout()
plt.savefig('DMEnergyChange')
plt.clf()