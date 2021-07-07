#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
maxq=np.loadtxt('../../Outfile/Turbulence/Star_WENO_MaxQ_0.dat',comments='"')
totalq=np.loadtxt('../../Outfile/Turbulence/Star_WENO_TotalQ_0.dat',comments='"')

#plot#
plt.plot(maxq[:,0], maxq[:,1])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Turbulence (CODE)')
plt.title('Maximum Turbulence Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('MaxTurbulence')
plt.clf()

#plot#
plt.plot(totalq[:,0], totalq[:,1])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Turbulence (CODE)')
plt.title('Total Turbulence Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('TotalTurbulence')
plt.clf()
