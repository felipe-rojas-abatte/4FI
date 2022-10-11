import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate
import sys
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d
import math

#------ Parametros -------
MW = 80.385
v = 246
g = 2*MW/v
MH = 125.7
GHT = 0.004255
mb = 4.18
mtau = 1.776 
mtop = 172
Nc = 3

#------ Scalar mass ------
M_phi = np.linspace(10**2, 10**5, 10000)

#--- We define the delta Width decay for the scalar and pseudo-scalar particle in the 2HDM. The parameters of the function are: Color factor, Cut-off, new coupling constant, scalar mass, mass of fermion in the final state
 
def delta_WDS(N,Lam,y_new,Ms,Mf):
	v=246
	ys=Mf*y_new/v	 
	WDS = N*(3/128)*(ys)**2*(Ms)/(np.pi*Lam**2)*(Ms**2-2*Mf**2)*(1-(2*Mf/Ms)**2)**(3/2)*np.log(Lam**2/Ms**2)
	return WDS

def delta_WDP(N,Lam,y_new,Mp,Mf):
	v=246
	yp=Mf*y_new/v 
	WDP = N*(3/128)*(yp)**2*(Mp)/(np.pi*Lam**2)*(Mp**2-6*Mf**2)*(1-(2*Mf/Mp)**2)**(1/2)*np.log(Lam**2/Mp**2)
	return WDP

#------ Total width decay of heavy scalars
## H0
(Mh2,WD2)=np.genfromtxt("/home/felipe/Dropbox/4_Fermiones/plots_python/Mh2_WD.dat",dtype=float,unpack=True,skip_header=False) 
func_wdH0=interp1d(Mh2, WD2, kind='linear')

## A0
(Mh3,WD3)=np.genfromtxt("/home/felipe/Dropbox/4_Fermiones/plots_python/Mh3_WD.dat",dtype=float,unpack=True,skip_header=False) 
func_wdA0=interp1d(Mh3, WD3, kind='linear')


#--- Delta Branching Ratio (phi->t tbar) as a function of (pseudo)scalar mass in 2HDM ------------

fig1 = plt.subplots()
plt.tick_params(axis='both',labelsize=12)
plt.title("2 Higgs doublet Model \n Type I (unconstrained by flauvor-physics)")
plt.xlabel("$M_{\\varphi}$ (GeV)", fontsize=15)
plt.ylabel("$\\delta Br^{4F} (\\varphi\\rightarrow t\\overline{t})$ (GeV)", fontsize=15)
#plt.xscale('log')
plt.yscale('log')
plt.xscale('log')
plt.xlim(100,100000)
plt.ylim(0.00001,0.1)
plt.plot(M_phi, delta_WDS(Nc,3000,2.69,M_phi,mtop)/func_wdH0(M_phi), linestyle='-', color='k', linewidth=1.2, label='$\\Lambda = $3 TeV' )
plt.plot(M_phi, delta_WDS(Nc,15000,2.69,M_phi,mtop)/func_wdH0(M_phi), linestyle='-', color='b', linewidth=1.2, label='$\\Lambda = $15 TeV')
plt.plot(M_phi, delta_WDS(Nc,30000,2.69,M_phi,mtop)/func_wdH0(M_phi), linestyle='-', color='r', linewidth=1.2, label='$\\Lambda = $30 TeV')

plt.plot(M_phi, delta_WDP(Nc,3000,2.77,M_phi,mtop)/func_wdA0(M_phi), linestyle='--', color='k', linewidth=2, label='$\\Lambda = $3 TeV' )
plt.plot(M_phi, delta_WDP(Nc,15000,2.77,M_phi,mtop)/func_wdA0(M_phi), linestyle='--', color='b', linewidth=2, label='$\\Lambda = $15 TeV')
plt.plot(M_phi, delta_WDP(Nc,30000,2.77,M_phi,mtop)/func_wdA0(M_phi), linestyle='--', color='r', linewidth=2, label='$\\Lambda = $30 TeV')


#Agregamos la leyenda
plt.legend(loc="upper right", ncol=2, title="scalar H              pseudo-scalar A", handlelength=2.5, borderaxespad=0.1, fancybox=True, shadow=True, fontsize = 10, labelspacing=0.1, handletextpad=0.5) # bbox_to_anchor=(1.5, 1.2))

plt.text(20000, 0.005, '$y_s=2.69; \\quad  y_p=2.77$', verticalalignment='bottom', horizontalalignment='center', color='black', fontsize=12)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('Br_toptop_2HDM.pdf')

#--- Delta Branching Ratio (phi->b bbar) as a function of (pseudo)scalar mass in 2HDM ------------

fig1 = plt.subplots()
plt.tick_params(axis='both',labelsize=15)
plt.title("2 Higgs doublet Model \n Type I (unconstrained by flauvor-physics)")
plt.xlabel("$M_{\\varphi}$ (GeV)", fontsize=15)
plt.ylabel("$\\delta Br^{4F} (\\varphi\\rightarrow b\\overline{b})$ (GeV)", fontsize=15)
#plt.xscale('log')
plt.yscale('log')
plt.xscale('log')
plt.xlim(400,100000)
#plt.ylim(0.00001,0.1)
plt.plot(M_phi, delta_WDS(Nc,3000,2.69,M_phi,mb)/func_wdH0(M_phi), linestyle='-', color='k', linewidth=1.2, label='$\\Lambda = $3 TeV' )
plt.plot(M_phi, delta_WDS(Nc,15000,2.69,M_phi,mb)/func_wdH0(M_phi), linestyle='-', color='b', linewidth=1.2, label='$\\Lambda = $15 TeV')
plt.plot(M_phi, delta_WDS(Nc,30000,2.69,M_phi,mb)/func_wdH0(M_phi), linestyle='-', color='r', linewidth=1.2, label='$\\Lambda = $30 TeV')

plt.plot(M_phi, delta_WDP(Nc,3000,2.77,M_phi,mb)/func_wdA0(M_phi), linestyle='--', color='k', linewidth=2, label='$\\Lambda = $3 TeV' )
plt.plot(M_phi, delta_WDP(Nc,15000,2.77,M_phi,mb)/func_wdA0(M_phi), linestyle='--', color='b', linewidth=2, label='$\\Lambda = $15 TeV')
plt.plot(M_phi, delta_WDP(Nc,30000,2.77,M_phi,mb)/func_wdA0(M_phi), linestyle='--', color='r', linewidth=2, label='$\\Lambda = $30 TeV')


#Agregamos la leyenda
plt.legend(loc="upper left", ncol=2, title="scalar H              pseudo-scalar A", handlelength=2.5, borderaxespad=0.1, fancybox=True, shadow=True, fontsize = 10, labelspacing=0.1, handletextpad=0.1) # bbox_to_anchor=(1.5, 1.2))

plt.text(370, 20, '$y_s=2.69; y_p=2.77$', verticalalignment='bottom', horizontalalignment='center', color='black', fontsize=16)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('Br_bb_2HDM.pdf')


