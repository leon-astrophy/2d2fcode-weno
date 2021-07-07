#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
temp2=np.loadtxt('../../Outfile/Hydro/Star_WENO_CentralTemperature_NM_0.dat',comments='"')

#plot#
plt.plot(temp2[:,0], temp2[:,1], label='Central Temp')
plt.plot(temp2[:,0], temp2[:,2], label='Maximum Temp')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Temperature (CODE)')
plt.title('Temperature Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('NMTemp')
plt.clf()

