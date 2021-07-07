#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
gwave=np.loadtxt('../../Outfile/GWave/Star_WENO_GravQ_0_NM.dat',comments='"')

#plot#
plt.plot(gwave[:,0], gwave[:,1], label='ae220')
plt.plot(gwave[:,0], gwave[:,2], label='hTT')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('GWave (CODE)')
plt.title('GWave Versus Time')
plt.grid(True)
plt.tight_layout()
plt.savefig('GWave')
plt.clf()
