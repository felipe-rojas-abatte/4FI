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

# ------ Escala de energia ------
L = np.linspace(10**3, 10**5, 10000)

#------ Delta Branching 4 Fermiones Br(H->ff)
#  We define the delta Width decay for the Higgs Boson. The parameters of the function are: Color factor, Cut-off, mass of fermion in the final state

def dBRH4F(N,Lam,Mf): 
	BRH4F = N*(3/512)*(MH)/(np.pi*GHT)*((g*Mf)/(MW*Lam))**2 *(MH**2-2*Mf**2)*(1-(2*Mf/MH)**2)**(3/2)*np.log(Lam**2/MH**2)
	return BRH4F

#---------- Delta Branchign ratio Br(h->bb) and Br(h->tautau)--------------

fig1 = plt.subplots()
plt.tick_params(axis='both',labelsize=12)
#Nombre de los ejes 
plt.xlabel("$\\Lambda$ (GeV)", fontsize=15)
plt.ylabel("$\\delta Br^{4F} \\equiv \\frac{\\delta\\Gamma^{4F}(h\\rightarrow XX)}{\\Gamma(h \\rightarrow all)}$", fontsize=15)
plt.xlim(1000,100000)
plt.ylim(0.0000001,0.1)
#Escala de los ejes
plt.xscale('log')
plt.yscale('log')

plt.plot(L, dBRH4F(Nc,L,mb), linestyle='-',color='k',linewidth=2, label='$h \\rightarrow b\\overline{b}$')
plt.plot(L, dBRH4F(1,L,mtau), linestyle='--',color='k',linewidth=2, label='$h \\rightarrow \\tau^+ \\tau^-$')

#Agregamos la leyenda
plt.legend(loc="upper right", handlelength=1.5,borderaxespad=0.5,fancybox=True, shadow=True,fontsize = 11,labelspacing=0.5,handletextpad=0.7)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('Br4f_AA.pdf')


print("\n Correction 4FI Br(h->bb) a 1 TeV: % ",dBRH4F(Nc,1000,mb)*100)
print(" Correction 4FI Br(h->tautau) A 1 TeV: % ",dBRH4F(1,1000,mtau)*100)
print("\n Correction 4FI Br(h->bb) a 30 TeV: % ",dBRH4F(Nc,30000,mb)*100)
print(" Correction 4FI Br(h->tautau) A 30 TeV: % ",dBRH4F(1,30000,mtau)*100)
print("\n")

######################################################################
#----------------- sigma x Br ---------------------
#--- seccion eficaz en (fb) para ILC en las diferentes
csZH_250 = 318.0; csZH_500 = 95.5; csZH_1T  = 22.3
csWH_250 = 36.6; csWH_500 = 163.0; csWH_1T  = 425.0
cs_sum_250 = csZH_250 + csWH_250
cs_sum_500 = csZH_500 + csWH_500
cs_sum_1T = csZH_1T + csWH_1T

#--- Luminosidad ILC en (fb-1)
L_250 = 250; L_500 = 500; L_1T  = 1000

#--- Ancho decaimiento Higgs segun modelo estandar en GeV
GSM = 0.004255

#--- Br(H->bb) y Br(H->tau+tau-) segun modelo estandar (1310.0763v3.pdf page 19)
BRbbSM = 0.577
BRttSM = 0.0632

#--- numero de eventos: cross section x Br_4F x L
NE_250_ZH = csZH_250*L_250*dBRH4F(Nc,L,mb)
NE_500_ZH = csZH_500*L_500*dBRH4F(Nc,L,mb)
NE_1T_ZH = csZH_1T*L_1T*dBRH4F(Nc,L,mb)

NE_250_WH = csWH_250*L_250*dBRH4F(Nc,L,mb)
NE_500_WH = csWH_500*L_500*dBRH4F(Nc,L,mb)
NE_1T_WH = csWH_1T*L_1T*dBRH4F(Nc,L,mb)

NE_250 = cs_sum_250*L_250*dBRH4F(Nc,L,mb)
NE_500 = cs_sum_500*L_500*dBRH4F(Nc,L,mb)
NE_1T = cs_sum_1T*L_1T*dBRH4F(Nc,L,mb)

#--- raiz numero eventos Modelo Estandar
y_250 = math.sqrt(cs_sum_250*L_250*BRbbSM) 
y_500 = math.sqrt(cs_sum_500*L_500*BRbbSM)
y_1T = math.sqrt(cs_sum_1T*L_1T*BRbbSM)  

#--- Significancia
S_250 = NE_250 / y_250 * 0.3
S_500 = NE_500 / y_500 * 0.3
S_1T = NE_1T / y_1T * 0.3

#-------- cross section x Br a diferentes luminosidades SUMADOS ------------

fig4 = plt.subplots()
plt.tick_params(axis='both',labelsize=12)
#Name of axes 
plt.xlabel("$\\Lambda$ (GeV)", fontsize=15)
plt.ylabel("$S_L \\equiv \\frac{\\sigma \\cdot L \\cdot Br^{4F}}{\\sqrt{\\sigma \\cdot L \\cdot Br^{SM}}}\\epsilon_f$ ", fontsize=15)
#plt.title('Production Cross section $\\sigma(e^+e^- \\rightarrow XH)$ at ILC, \n HZ + WW Fusion', fontsize=15)
plt.xlim(1000,100000)
plt.ylim(0,4)
plt.xscale('log')
#plt.yscale('log')

plt.plot(L, S_250, linestyle='-.',color='g',linewidth=2, label='S$_L$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(L, S_500, linestyle='-',color='r',linewidth=2, label='S$_L$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(L, S_1T, linestyle='--',color='k',linewidth=2, label='S$_L$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')

plt.legend(loc="upper right", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 11,labelspacing=0.1,handletextpad=0.5)

#plt.grid()
plt.tight_layout()
#plt.grid()
plt.savefig('S.pdf')


#------ Scalar mass ------
M_phi = np.linspace(10**2, 10**5, 10000)

#--- We define the delta Width decay for the scalar and pseudo-scalar particle in the 2HDM. The parameters of the function are: Color factor, Cut-off, new coupling constant, scalar mass, mass of fermion in the final state
 
def delta_WDS(N,Lam,y_new,Ms,Mf):
	v=246
	ys=Mf*y_new/v	 
	WDS = N*(3/128)*Ms*ys**2.0/(np.pi*Lam**2.0)*(Ms**2.0-2*Mf**2.0)*(1-(2*Mf/Ms)**2.0)**(3.0/2)*np.log(Lam**2.0/Ms**2.0)
	return WDS

def delta_WDP(N,Lam,y_new,Mp,Mf):
	v=246
	yp=Mf*y_new/v 
	WDP = N*(3/128)*(yp)**2*(Mp)/(np.pi*Lam**2)*(Mp**2-6*Mf**2)*(1-(2*Mf/Mp)**2)**(1/2)*np.log(Lam**2/Mp**2)
	return WDP

#################### Plots #######################

#-------- Delta Width Decay (H->bb) and (H->tautau) as a function of the Lambda in 2HDM ------------
fig1 = plt.subplots()
plt.tick_params(axis='both',labelsize=12)
plt.xlabel("$\\Lambda$ (GeV)", fontsize=15)
plt.ylabel("$\\delta \\Gamma^{4F} (h\\rightarrow f\\overline{f})$ (GeV)", fontsize=15)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1000,100000)

HDM1, = plt.plot(L, delta_WDS(Nc,L,-0.99,MH,mb), linestyle='-', color='k', linewidth=1.2, label='$h \\rightarrow b \\overline{b}$' )
HDM2, = plt.plot(L, delta_WDS(Nc,L,-0.99,MH,mtau), linestyle='-', color='b', linewidth=1.2, label='$h \\rightarrow \\tau^+\\tau^-$')
HDM3, = plt.plot(L, delta_WDS(Nc,L,-0.91,MH,mb), linestyle='--', color='k', linewidth=2, label='$h \\rightarrow b \\overline{b}$' )
HDM4, = plt.plot(L, delta_WDS(Nc,L,-0.91,MH,mtau), linestyle='--', color='b', linewidth=2, label='$h \\rightarrow \\tau^+\\tau^-$')

#Creamos la primera leyenda
legend1 = plt.legend(handles=[HDM1,HDM2,HDM3,HDM4], loc="upper right", ncol=2, title="                  2HDM \n constrained       unconstranied", handlelength=2.5, borderaxespad=0.1, fancybox=True, shadow=True, fontsize = 11, labelspacing=0.1, handletextpad=0.1) # bbox_to_anchor=(1.5, 1.2))

#Agregamos la primera leyenda
plt.gca().add_artist(legend1)

#Creamos la segunda leyenda
SM1, = plt.plot(L, delta_WDS(Nc,L,1,MH,mb), linestyle='-.', color='r', linewidth=2, label='$h \\rightarrow b \\overline{b}$' )
SM2, = plt.plot(L, delta_WDS(Nc,L,1,MH,mtau), linestyle='-.', color='g', linewidth=2, label='$h \\rightarrow \\tau^+\\tau^-$')

#Agregamos segunda leyenda
plt.legend(handles=[SM1,SM2], loc="lower left", ncol=1, title=" SM ", handlelength=3, borderaxespad=0.1, fancybox=True, shadow=True, fontsize = 11, labelspacing=0.1, handletextpad=0.1)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('WDH_2HDM.pdf')

#Comparacion entre las 2 funciones
#print("SM1: ", dBRH4F(Nc,3000,mb)*GHT)
#print("SM2: ", delta_WDS(Nc,3000,1,MH,mb))


#-------- Delta Width decay as a function of (pseudo)scalar mass in 2HDM ------------

fig1 = plt.subplots()
plt.tick_params(axis='both',labelsize=12)
plt.title("2 Higgs doublet Model \n Type I (unconstrained by flauvor-physics)")
plt.xlabel("$M_{\\varphi}$ (GeV)", fontsize=15)
plt.ylabel("$\\delta \\Gamma^{4F} (\\varphi\\rightarrow t\\overline{t})$ (GeV)", fontsize=15)
#plt.xscale('log')
plt.yscale('log')
plt.xscale('log')
plt.xlim(100,100000)
plt.ylim(0.0001,10000)
plt.plot(M_phi, delta_WDS(Nc,3000,2.69,M_phi,mtop), linestyle='-', color='k', linewidth=1.2, label='$\\Lambda = $3 TeV' )
plt.plot(M_phi, delta_WDS(Nc,15000,2.69,M_phi,mtop), linestyle='-', color='b', linewidth=1.2, label='$\\Lambda = $15 TeV')
plt.plot(M_phi, delta_WDS(Nc,30000,2.69,M_phi,mtop), linestyle='-', color='r', linewidth=1.2, label='$\\Lambda = $30 TeV')

plt.plot(M_phi, delta_WDP(Nc,3000,2.77,M_phi,mtop), linestyle='--', color='k', linewidth=2, label='$\\Lambda = $3 TeV' )
plt.plot(M_phi, delta_WDP(Nc,15000,2.77,M_phi,mtop), linestyle='--', color='b', linewidth=2, label='$\\Lambda = $15 TeV')
plt.plot(M_phi, delta_WDP(Nc,30000,2.77,M_phi,mtop), linestyle='--', color='r', linewidth=2, label='$\\Lambda = $30 TeV')

# Pseudo-scalar A^o
## y_b y_tau = 0.36 , y_t = 2.77   unconstrained
## y_b y_tau = 0.50 , y_t = 2     constrained

# Scalar H^o
## y_b y_tau = 0.37 , y_t = 2.69  unconstrained
## y_b y_tau = 0.64 , y_t = 1.84  constrained
 

#Agregamos la leyenda
plt.legend(loc="upper left", ncol=2, title="scalar H              pseudo-scalar A", handlelength=2.5, borderaxespad=0.1, fancybox=True, shadow=True, fontsize = 10, labelspacing=0.1, handletextpad=0.5) # bbox_to_anchor=(1.5, 1.2))

plt.text(400, 40, '$y_s=2.69; \\quad y_p=2.77$', verticalalignment='bottom', horizontalalignment='center', color='black', fontsize=12)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('WD_toptop_2HDM.pdf')


