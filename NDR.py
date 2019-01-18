#%matplotlib inline

import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

import time
start_time = time.time()
matplotlib.rcParams.update({'font.size': 16})

T=4.5e-3 # temperature
Ed=3.270 # charging energy for spin down
Eu=3.149 # charging energy for spin up
muR=3.43 # chemical potential of the right electrode

wL=1e10  # transition rates for left electrode
wR=1e10  # transition rates for right electrode




def f(muL,muR,Vg,spin, electrode,direction):
 if spin == "up":
 	E = Eu
 else: E = Ed
 if electrode == "L":
 	mu = muL
 else: mu = muR
 if direction == "in":
 	sign = 1
 else: sign = -1
 return 1.0/(1.0+np.exp(sign*(E-Vg-mu)/T))


    
def current(muL,muR,Vg):
 Wp0L=wL*(f(muL,muR,Vg, "up", "L", "in"))
 Wp0R=wR*(f(muL,muR,Vg, "up", "R", "in"))
 Wm0L=wL*(f(muL,muR,Vg,"dn","L","in"))
 Wm0R=wR*(f(muL,muR,Vg,"dn","R","in"))
 
 W0pL=wL*(f(muL,muR,Vg,"up","L","out"))
 W0pR=wR*(f(muL,muR,Vg,"up","R","out"))
 W0mL=wL*(f(muL,muR,Vg,"dn","L","out"))
 W0mR=wR*(f(muL,muR,Vg,"dn","R","out"))
 
 W0p=W0pL+W0pR
 W0m=W0mL+W0mR
 Wp0=Wp0L+Wp0R
 Wm0=Wm0L+Wm0R
 
 
 eta=W0p*W0m+Wp0*W0m+W0p*Wm0
 rho0=W0p*W0m/eta
 rhop=Wp0*W0m/eta
 rhom=W0p*Wm0/eta
 
 Iup=W0pR*rhop-Wp0R*rho0
 Idn=W0mR*rhom-Wm0R*rho0
 Itot=Iup+Idn
 
 return [Iup,rho0,rhop,rhom,Idn,Itot]
	    




def diamond():
 def fmt(x, pos):
     return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')
 Ngrid=201  # number of points for plotting
 temp=0.2
 Lim=1.3
 
 xlist = np.linspace(Eu-temp, Eu+temp, Ngrid,endpoint=True)
 ylist = np.linspace(0, temp*2, Ngrid,endpoint=True)
 X, Y = np.meshgrid(xlist, ylist)
 Iup=np.zeros((Ngrid,Ngrid))
 Idn=np.zeros((Ngrid,Ngrid))
 Itot=np.zeros((Ngrid,Ngrid))
#	dIup=np.zeros((Ngrid,Ngrid))
#	dIdn=np.zeros((Ngrid,Ngrid))
#	dItot=np.zeros((Ngrid,Ngrid))
	
 for i in range(Ngrid):
  for j in range(Ngrid):
   temp=current(Y[j][0]/2.0,-Y[j][0]/2.0,X[0][i])
   Iup[j,i]=temp[0]*1.60217662*1e-10
   Idn[j,i]=temp[4]*1.60217662*1e-10
   Itot[j,i]=temp[5]*1.60217662*1e-10

 fig=plt.figure()
 Ncolor=16
 
 p1=fig.add_subplot(221)
 cp = p1.contourf(X,Y , Iup,np.linspace(0,Lim, Ncolor),cmap=plt.cm.Reds)
 
 plt.title('(a) Up')
 plt.xticks(np.arange(3.0,3.4,0.1))
 plt.xlabel('$V_g$ (V)')
 plt.ylabel('$V_b$ (V)')
 p1.yaxis.set_major_locator(MultipleLocator(0.1))

#	plt.gca().set_aspect('equal', adjustable='box')
	
 p2=fig.add_subplot(222)
 cp = p2.contourf(X,Y ,Idn,np.linspace(0,Lim, Ncolor),cmap=plt.cm.Reds)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))

 plt.title('(b) Down')
 plt.xlabel('$V_g$ (V)')
 plt.ylabel('$V_b$ (V)')
 p2.yaxis.set_major_locator(MultipleLocator(0.1))
 plt.xticks(np.arange(3.0,3.4,0.1))
 cbaxes = fig.add_axes([0.85, 0.6, 0.03, 0.33]) 
#plt.colorbar(cp, cax = cbaxes,format=ticker.FuncFormatter(fmt))  
 plt.colorbar(cp, cax = cbaxes)  

#plt.gca().set_aspect('eVgual', adjustable='box')
	
 p3=fig.add_subplot(223)
 cp = p3.contourf( X,Y, Itot,np.linspace(0,Lim, Ncolor),cmap=plt.cm.Reds)
#plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
 plt.title('(c) Total')
 plt.xticks(np.arange(3.0,3.4,0.1))
 plt.xlabel('$V_g$ (V)')
 plt.ylabel('$V_b$ (V)')
 p3.yaxis.set_major_locator(MultipleLocator(0.1))

#plt.gca().set_aspect('equal', adjustable='box')
	
#p4=fig.add_subplot(224)
 Nplot=100
 Bias=np.zeros(Nplot)
 Iup=np.zeros(Nplot)
 rho0=np.zeros(Nplot)
 rhop=np.zeros(Nplot)
 rhom=np.zeros(Nplot)
 #muC=3.149
 Vg=Eu+0.01
 biasrange=0.4
 
 for i in range(Nplot):
  Bias[i]=biasrange/Nplot*i
  Iup[i]=current(Bias[i]/2.0,-Bias[i]/2.0,Vg)[0]*1.60217662*1e-10
  rho0[i]=current(Bias[i]/2.0,-Bias[i]/2.0,Vg)[1]
  rhop[i]=current(Bias[i]/2.0,-Bias[i]/2.0,Vg)[2]
  rhom[i]=current(Bias[i]/2.0,-Bias[i]/2.0,Vg)[3]
 
 
 
 ax1 = fig.add_subplot(224)
 ax2 = ax1.twinx()
 
 line1=ax1.plot(Bias,Iup,label=r'$I_\uparrow$',linewidth=2)
 ax1.set_ylim([0,1])
 ax2.set_ylim([0,1])
 
 line2=ax2.plot(Bias,rhop,label=r'$\rho_\uparrow$',linestyle='dashed',linewidth=2)
 line3=ax2.plot(Bias,rho0,label=r'$\rho_0$',linestyle='dashdot',linewidth=2)
 line4=ax2.plot(Bias,rhom,label=r'$\rho_\downarrow$',linestyle='dotted',linewidth=2)
 ax1.set_xlabel('$V_b$ (V)')
 ax1.set_ylabel('$I$ (nA)')
 ax2.set_ylabel(r'$\rho$')
 ax2.yaxis.set_major_locator(MultipleLocator(0.5))
 
 plt.xticks([0,0.2,0.4])
 
 
 line=line1+line2+line3+line4
 labs = [l.get_label() for l in line]
 
 # Shrink current axis by 20%
 box = ax1.get_position()
 ax1.set_position([box.x0, box.y0, box.width * 1, box.height])
 ax2.set_position([box.x0, box.y0, box.width * 1, box.height])

# Put a legend to the right of the current axis
 ax1.legend(line,labs,ncol=1,loc='center left', bbox_to_anchor=(1.2, 0.5),frameon=False)
#	ax1.legend(line,labs,bbox_to_anchor=(0., 1.0, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)

 plt.title('(d) NDC')
 plt.tight_layout()
 plt.show()
	
	

 return

#centerchange1(0,muR)	
diamond()
print("time:", (time.time() - start_time)/60.0)
