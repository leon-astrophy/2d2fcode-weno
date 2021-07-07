#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
cye2=np.loadtxt('../../Outfile/Hydro/Star_WENO_CentralYe_NM_0.dat',comments='"')
cabar2=np.loadtxt('../../Outfile/Isotope/Star_WENO_CentralABar_0.dat',comments='"')
czbar2=np.loadtxt('../../Outfile/Isotope/Star_WENO_CentralZBar_0.dat',comments='"')
lum=np.loadtxt('../../Outfile/Isotope/Star_WENO_Luminosity_0.dat',comments='"')
xmass=np.loadtxt('../../Outfile/Isotope/Star_WENO_XMass_Selected_0.dat',comments='"')

#plot#
plt.plot(xmass[:,0], xmass[:,1], label='He4')
plt.plot(xmass[:,0], xmass[:,2], label='C12')
plt.plot(xmass[:,0], xmass[:,3], label='O16')
plt.plot(xmass[:,0], xmass[:,4], label='Ne20')
plt.plot(xmass[:,0], xmass[:,5], label='Mg24')
plt.plot(xmass[:,0], xmass[:,6], label='Si28')
plt.plot(xmass[:,0], xmass[:,7], label='Ni56')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Mass (CODE)')
plt.title('Isotope Mass Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('IsotopeMass')
plt.clf()

#plot#
plt.plot(cye2[:,0], cye2[:,1], label='Central Ye')
plt.plot(cye2[:,0], cye2[:,2], label='Maximum Ye')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Ye (CODE)')
plt.title('Ye Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('Ye')
plt.clf()

#plot#
plt.plot(cabar2[:,0], cabar2[:,1], label='Central Abar')
plt.plot(czbar2[:,0], czbar2[:,1], label='Central Zbar')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('AZBar (CODE)')
plt.title('AZBar Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('AZBar')
plt.clf()

#plot#
plt.plot(lum[:,0], lum[:,1], label='Total')
plt.plot(lum[:,0], lum[:,2], label='Flame')
plt.plot(lum[:,0], lum[:,2], label='Deton')
plt.plot(lum[:,0], lum[:,3], label='Burn')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Luminosity (CODE)')
plt.title('Luminosity Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('lum')
plt.clf()