#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
mass2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Mass_NM_0.dat',comments='"')
energy2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Energy_NM_0.dat',comments='"')
totalenergy=np.loadtxt('../../Outfile/Hydro/Star_WENO_Energy_Total_0.dat',comments='"')
centralrho2=np.loadtxt('../../Outfile/Hydro/Star_WENO_CentralDensity_NM_0.dat',comments='"')

#plot#
plt.plot(centralrho2[:,0], centralrho2[:,1], label='Central Density')
plt.plot(centralrho2[:,0], centralrho2[:,2], label='Maximum Density')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Density (CODE)')
plt.title('Density Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('NMCentralDensity')
plt.clf()

#plot#
plt.plot(mass2[:,0], np.abs((mass2[:,1] - mass2[0,1])/(mass2[0,1])))
plt.legend(loc="lower right")
plt.xlabel('Time (CODE)')
plt.ylabel('Fractional Change')
plt.title('Mass Conservation')
plt.grid(True)
plt.tight_layout()
plt.savefig('NMMass')
plt.clf()

#plot#
plt.plot(energy2[:,0], energy2[:,1], label='Total NM Energy')
plt.plot(energy2[:,0], energy2[:,2], label='Kinetic Energy')
plt.plot(energy2[:,0], energy2[:,3], label='Internal Energy')
plt.plot(energy2[:,0], energy2[:,4], label='Gravitational Energy')
plt.legend(loc="lower right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy (CODE)')
plt.title('Energy Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('NMEnergy')
plt.clf()

#plot#
plt.plot(energy2[:,0], np.abs((energy2[:,1] - energy2[0,1])/(energy2[0,1])), label='Total NM Energy')
plt.plot(totalenergy[:,0], np.abs((totalenergy[:,1] - totalenergy[0,1])/(totalenergy[0,1])), label='Total Energy')
plt.legend(loc="lower right")
plt.xlabel('Time (CODE)')
plt.ylabel('Fractional Change')
plt.title('Energy Conservation')
plt.grid(True)
plt.tight_layout()
plt.savefig('TotalEnergyChange')
plt.clf()

#plot#
plt.plot(totalenergy[:,0], totalenergy[:,2], label='Energy Input')
plt.legend(loc="lower right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy (CODE)')
plt.title('Commutative Energy Generation')
plt.grid(True)
plt.tight_layout()
plt.savefig('EnergyGeneration')
plt.clf()