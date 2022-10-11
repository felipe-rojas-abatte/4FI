import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate
import sys
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d
import math
import os

current_dir = os.getcwd()
data_file = '../datos/'

#Branching Ratios H->bb y H->tt contribucion 4 fermiones
Lam_bb, BR_bb_4F = np.genfromtxt(data_file+'dFbb_ax.dat', dtype=float, unpack=True, skip_header=True) 
Lam_tt, BR_tt_4F = np.genfromtxt(data_file+'dFtt_ax.dat', dtype=float, unpack=True, skip_header=True) 

#New scalar width decay into bb (S->bb), contribucion 4 fermiones
MH_bb3T, W_bb3T = np.genfromtxt(data_file+'dFSAbb_3T.dat', dtype=float, unpack=True, skip_header=True) 
MH_bb15T, W_bb15T = np.genfromtxt(data_file+'dFSAbb_15T.dat', dtype=float, unpack=True, skip_header=True) 
MH_bb30T, W_bb30T = np.genfromtxt(data_file+'dFSAbb_30T.dat', dtype=float, unpack=True, skip_header=True) 

#New scalar width decay into bb (S->tautau), contribucion 4 fermiones
MH_tt3T, W_tt3T = np.genfromtxt(data_file+'dFSAtt_3T.dat', dtype=float, unpack=True, skip_header=True) 
MH_tt15T, W_tt15T = np.genfromtxt(data_file+'dFSAtt_15T.dat', dtype=float, unpack=True, skip_header=True) 
MH_tt30T, W_tt30T = np.genfromtxt(data_file+'dFSAtt_30T.dat', dtype=float, unpack=True, skip_header=True) 


#--------------- Branchign ratio interaccion de 4 fermiones---------------------
fig1 = plt.subplots()
#plt.tick_params(axis='both',labelsize=15)
#Nombre de los ejes y titulo del grafico
plt.xlabel("$\\Lambda$ (GeV)", fontsize=20)
plt.ylabel("$\\delta Br^{4F} \\equiv \\frac{\\delta\\Gamma^{4F}(H\\rightarrow XX)}{\\Gamma(H \\rightarrow all)}$", fontsize=18)
#plt.title('axial--axial interaction', fontsize=15)
#plt.xlim(200,1000)
#plt.ylim(10,0.000001)
#Escala de los ejes
plt.xscale('log')
plt.yscale('log')

#Graficamos los pares ordenados Lam_xx y BR_xx_4F
plt.plot(Lam_bb, BR_bb_4F, linestyle='-',color='k',linewidth=2, label='$H \\rightarrow b\\overline{b}$')
plt.plot(Lam_tt, BR_tt_4F, linestyle='--',color='k',linewidth=2, label='$H \\rightarrow \\tau^+ \\tau^-$')

#Agregamos la leyenda
plt.legend(loc="upper right", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 17,labelspacing=0.1,handletextpad=0.1)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('Br4f_AA.pdf')

#----------------- sigma x Br ---------------------
#--- seccion eficaz en (fb)
csZH_250 = 318.0; csZH_500 = 95.5; csZH_1T  = 22.3
csWH_250 = 36.6; csWH_500 = 163.0; csWH_1T  = 425.0
cs_sum_250 = csZH_250 + csWH_250
cs_sum_500 = csZH_500 + csWH_500
cs_sum_1T = csZH_1T + csWH_1T

#--- Luminosidad ILC en (GeV)
L_250 = 250; L_500 = 500; L_1T  = 1000

#--- Ancho decaimiento Higgs segun modelo estandar en GeV
GSM = 0.004255

#--- Br(H->bb) y Br(H->tt) segun modelo estandar (1310.0763v3.pdf page 19)
BRbbSM = 0.0577
BRttSM = 0.0632

#--- cross section x Br_4F x L
NE_250_ZH = csZH_250*L_250*BR_bb_4F
NE_500_ZH = csZH_500*L_500*BR_bb_4F
NE_1T_ZH = csZH_1T*L_1T*BR_bb_4F

NE_250_WH = csWH_250*BR_bb_4F*L_250
NE_500_WH = csWH_500*BR_bb_4F*L_500
NE_1T_WH = csWH_1T*BR_bb_4F*L_1T

NE_250 = cs_sum_250*BR_bb_4F*L_250
NE_500 = cs_sum_500*BR_bb_4F*L_500
NE_1T = cs_sum_1T*BR_bb_4F*L_1T

#--- backgraound ZH
x_back = [1000, 100000]
y_ZH_250 = math.sqrt(csZH_250*L_250*BRbbSM) 
y_ZH_500 = math.sqrt(csZH_500*L_500*BRbbSM) 
y_ZH_1T = math.sqrt(csZH_1T*L_1T*BRbbSM) 
y_back_250_ZH = [y_ZH_250,y_ZH_250]
y_back_500_ZH = [y_ZH_500,y_ZH_500]
y_back_1T_ZH = [y_ZH_1T,y_ZH_1T]

#--- backgraound WWH
y_WH_250 = math.sqrt(csWH_250*L_250*BRbbSM) 
y_WH_500 = math.sqrt(csWH_500*L_500*BRbbSM)
y_WH_1T = math.sqrt(csWH_1T*L_1T*BRbbSM)  
y_back_250_WH = [y_WH_250,y_WH_250]
y_back_500_WH = [y_WH_500,y_WH_500]
y_back_1T_WH = [y_WH_1T,y_WH_1T]

#--- background sumado
y_250 = math.sqrt(cs_sum_250*L_250*BRbbSM) 
y_500 = math.sqrt(cs_sum_500*L_500*BRbbSM)
y_1T = math.sqrt(cs_sum_1T*L_1T*BRbbSM)  
y_back_250 = [y_250,y_250]
y_back_500 = [y_500,y_500]
y_back_1T = [y_1T,y_1T]


#-------- cross section x Br a diferentes luminosidades ------------

fig2 = plt.subplots()
plt.tick_params(axis='both',labelsize=15)
#Name of axes 
plt.xlabel("$\\Lambda$ (GeV)", fontsize=20)
plt.ylabel("$\\sigma\\cdot L \\cdot Br$", fontsize=20)
plt.title('Production Cross section $\\sigma(e^+e^- \\rightarrow HZ)$ at ILC \n Higgs-strahlung', fontsize=15)
#plt.xlim(1000,10000)
#plt.ylim(1,1000)
plt.xscale('log')
plt.yscale('log')

plt.plot(Lam_bb, NE_250_ZH, linestyle='-',color='r',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(Lam_bb, NE_500_ZH, linestyle='-',color='b',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(Lam_bb, NE_1T_ZH, linestyle='-',color='g',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')
plt.plot(x_back, y_back_250_ZH, linestyle='-.',color='r',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(x_back, y_back_500_ZH, linestyle='-.',color='b',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(x_back, y_back_1T_ZH, linestyle='-.',color='g',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')

plt.legend(loc="lower left", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 12,labelspacing=0.1,handletextpad=0.1)

#plt.grid()
plt.tight_layout()
plt.savefig('ZH_BR.pdf')


#-------- cross section x Br a diferentes luminosidades ------------

fig3 = plt.subplots()
plt.tick_params(axis='both',labelsize=15)
#Name of axes 
plt.xlabel("$\\Lambda$ (GeV)", fontsize=20)
plt.ylabel("$\\sigma\\cdot L \\cdot Br$", fontsize=20)
plt.title('Production Cross section $\\sigma(e^+e^- \\rightarrow \\nu \\overline{\\nu} H)$ at ILC \n WW Fusion', fontsize=15)
#plt.xlim(200,1000)
#plt.ylim(0.01,10000)
plt.xscale('log')
plt.yscale('log')

plt.plot(Lam_bb, NE_250_WH, linestyle='-',color='r',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(Lam_bb, NE_500_WH, linestyle='-',color='b',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(Lam_bb, NE_1T_WH, linestyle='-',color='g',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')
plt.plot(x_back, y_back_250_WH, linestyle='-.',color='r',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(x_back, y_back_500_WH, linestyle='-.',color='b',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(x_back, y_back_1T_WH, linestyle='-.',color='g',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')


plt.legend(loc="lower left", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 11,labelspacing=0.1,handletextpad=0.1)

#plt.grid()
plt.tight_layout()
plt.savefig('WH_BR.pdf')

#-------- cross section x Br a diferentes luminosidades SUMADOS ------------

fig4 = plt.subplots()
plt.tick_params(axis='both',labelsize=15)
#Name of axes 
plt.xlabel("$\\Lambda$ (GeV)", fontsize=20)
plt.ylabel("$\\sigma\\cdot L \\cdot Br$", fontsize=20)
plt.title('Production Cross section $\\sigma(e^+e^- \\rightarrow XH)$ at ILC, \n HZ + WW Fusion', fontsize=15)
#plt.xlim(200,1000)
#plt.ylim(0.01,10000)
plt.xscale('log')
plt.yscale('log')

plt.plot(Lam_bb, NE_250, linestyle='--',color='g',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(Lam_bb, NE_500, linestyle='--',color='r',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(Lam_bb, NE_1T, linestyle='--',color='k',linewidth=2, label='$\\sigma \\cdot L \\cdot Br^{4F}$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')
plt.plot(x_back, y_back_250, linestyle='-',color='g',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(x_back, y_back_500, linestyle='-',color='r',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(x_back, y_back_1T, linestyle='-',color='k',linewidth=2, label='$\\sqrt{\\sigma \\cdot L \\cdot Br^{SM} }$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')


plt.legend(loc="lower left", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 11,labelspacing=0.1,handletextpad=0.1)

#plt.grid()
plt.tight_layout()
plt.savefig('ZH+WH_BR.pdf')


###### Nueva parte

ne_250 = NE_250/y_250
ne_500 = NE_500/y_500
ne_1T = NE_1T/y_1T

#-------- cross section x Br a diferentes luminosidades SUMADOS ------------

fig4 = plt.subplots()
plt.tick_params(axis='both',labelsize=15)
#Name of axes 
plt.xlabel("$\\Lambda$ (GeV)", fontsize=20)
plt.ylabel("$S_L \\equiv \\frac{\\sigma \\cdot L \\cdot Br^{4F}}{\\sqrt{\\sigma \\cdot L \\cdot Br^{SM}}}$ ", fontsize=20)
#plt.title('Production Cross section $\\sigma(e^+e^- \\rightarrow XH)$ at ILC, \n HZ + WW Fusion', fontsize=15)
plt.xlim(1000,50000)
plt.ylim(0,4)
plt.xscale('log')
#plt.yscale('log')

plt.plot(Lam_bb, ne_250, linestyle='-.',color='g',linewidth=2, label='S$_L$ with $L=250 (fb^{-1})$, $\\sqrt{s}=250$ GeV')
plt.plot(Lam_bb, ne_500, linestyle='-',color='r',linewidth=2, label='S$_L$ with $L=500 (fb^{-1})$, $\\sqrt{s}=500$ GeV')
plt.plot(Lam_bb, ne_1T, linestyle='--',color='k',linewidth=2, label='S$_L$ with $L=1000 (fb^{-1})$, $\\sqrt{s}=1000$ GeV')

plt.legend(loc="upper right", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 11,labelspacing=0.1,handletextpad=0.1)

#plt.grid()
plt.tight_layout()
plt.savefig('S.pdf')


#--------------- NEW SCALAR Width decay with 4 fermion interaction ---------------------
fig1 = plt.subplots()
#plt.tick_params(axis='both',labelsize=15)
#Nombre de los ejes y titulo del grafico
plt.xlabel("$M_{\\varphi}$ (GeV)", fontsize=20)
plt.ylabel("$\\delta\\Gamma^{4F}(\\varphi\\rightarrow \overline{b}b)$", fontsize=18)
#plt.title('axial--axial interaction', fontsize=15)
plt.xlim(100,10000)  # en GeV
plt.ylim(0.00000001,0.003)
#Escala de los ejes
plt.xscale('log')
plt.yscale('log')

plt.plot(MH_bb3T, W_bb3T, linestyle='-',color='k',linewidth=2, label='$\\Lambda = 3$ TeV')
plt.plot(MH_bb15T, W_bb15T, linestyle='--',color='b',linewidth=2, label='$\\Lambda = 15$ TeV')
plt.plot(MH_bb30T, W_bb30T, linestyle='-.',color='r',linewidth=2, label='$\\Lambda = 30$ TeV')


#Agregamos la leyenda
plt.legend(loc="upper left", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 17,labelspacing=0.1,handletextpad=0.1)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('W4f_bb_MS.pdf')

#--------------- NEW SCALAR Width decay with 4 fermion interaction ---------------------
fig1 = plt.subplots()
#plt.tick_params(axis='both',labelsize=15)
#Nombre de los ejes y titulo del grafico
plt.xlabel("$M_{\\varphi}$ (GeV)", fontsize=20)
plt.ylabel("$\\delta\\Gamma^{4F}(\\varphi\\rightarrow \\tau^+\\tau^-)$", fontsize=18)
#plt.title('axial--axial interaction', fontsize=15)
plt.xlim(100,10000)
plt.ylim(0.000000001,0.0003)
#Escala de los ejes
plt.xscale('log')
plt.yscale('log')

plt.plot(MH_tt3T, W_tt3T, linestyle='-',color='k',linewidth=2, label='$\\Lambda = 3$ TeV')
plt.plot(MH_tt15T, W_tt15T, linestyle='--',color='b',linewidth=2, label='$\\Lambda = 15$ TeV')
plt.plot(MH_tt30T, W_tt30T, linestyle='-.',color='r',linewidth=2, label='$\\Lambda = 30$ TeV')


#Agregamos la leyenda
plt.legend(loc="upper left", handlelength=2,borderaxespad=0.1,fancybox=True, shadow=True,fontsize = 17,labelspacing=0.1,handletextpad=0.1)

#Guardamos en un formato mas comprimido
plt.tight_layout()
#Exportamos el grafico como extension pdf
plt.savefig('W4f_tt_MS.pdf')



