import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate
import sys
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d
import os

current_dir = os.getcwd()
data_file = '../datos/'

#Read the data file and save it in the variables
#Cross sections files
En, cs = np.genfromtxt(data_file+"cs_eE_Zh.dat", dtype=float, unpack=True, skip_header=True) 
En2, cs2 = np.genfromtxt(data_file+"cs_eE_Wh.dat", dtype=float, unpack=True, skip_header=True) 

Ener = np.linspace(200, 1000, 800, endpoint=True)

#Interpolation of the CheckMATE analysis
SIG_Zh = interp1d(En, cs, kind='linear')
SIG_1 = SIG_Zh(Ener)
SIG_Wh = interp1d(En2, cs2, kind='linear')
SIG_2 = SIG_Wh(Ener)

SUM_SIG = SIG_1+SIG_2

fig = plt.subplots()
plt.tick_params(axis='both',labelsize=15)
#Name of axes 
plt.xlabel("$\\sqrt{s}\; (\\mathrm{GeV})$", fontsize=20)
plt.ylabel("$\\sigma(e^+e^-\\rightarrow Xh)\; (\\mathrm{fb})$", fontsize=20)
plt.xlim(200,1000)
#plt.ylim(0.01,10000)
#plt.yscale('log')

plt.plot(Ener, SIG_1, linestyle='-',color='r',linewidth=2, label='$Zh$ higgs-strahlung')
plt.plot(Ener, SIG_2, linestyle='-',color='g',linewidth=2, label='$WW$ fusion')
plt.plot(Ener, SUM_SIG, linestyle='-',color='k',linewidth=2, label='both')

#plt.text(60, 0.0000015, '', verticalalignment='top', horizontalalignment='center', color='black', fontsize=12,zorder=9)

plt.legend(loc="upper right", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 17,labelspacing=0.1,handletextpad=0.1)

plt.grid()
plt.tight_layout()
#plt.savefig('eE_Zh.pdf')
plt.savefig('eE_xx.pdf')
# plt.show()
plt.clf()



