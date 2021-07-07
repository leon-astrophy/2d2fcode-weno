#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
loss=np.loadtxt('../../Outfile/Neutrino/Star_WENO_Energy_Loss_0.dat',comments='"')
spec=np.loadtxt('../../Outfile/Neutrino/Star_WENO_NuSpec_0.dat',comments='"')

#plot#
plt.plot(loss[:,0], loss[:,1], label='Thermo')
#plt.plot(loss[:,0], loss[:,2], label='Ecap')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy Loss Rate (CODE)')
plt.title('Neutrino Energy Loss Versus Time')
plt.grid(True)
plt.savefig('NeuLoss')
plt.clf()

#plot#
plt.plot(spec[:,0], spec[:,1], label='1 MeV')
plt.plot(spec[:,0], spec[:,2], label='2 MeV')
plt.plot(spec[:,0], spec[:,3], label='3 MeV')
plt.plot(spec[:,0], spec[:,4], label='4 MeV')
plt.plot(spec[:,0], spec[:,5], label='5 MeV')
plt.plot(spec[:,0], spec[:,6], label='6 MeV')
plt.plot(spec[:,0], spec[:,7], label='7 MeV')
plt.plot(spec[:,0], spec[:,8], label='8 MeV')
plt.plot(spec[:,0], spec[:,9], label='9 MeV')
plt.plot(spec[:,0], spec[:,10], label='10 MeV')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Spectra (CODE)')
plt.title('Neutrino Spectra Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('NeuSpec')
plt.clf()

