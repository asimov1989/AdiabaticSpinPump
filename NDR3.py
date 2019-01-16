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
from numpy import linalg as LA
from scipy import linalg as LAs
from scipy.linalg import eig as eigscipy
import time
start_time = time.time()
matplotlib.rcParams.update({'font.size': 16})

global T,Eu,Ed,h,hz,wL,wR, a, b ,dchi,Tp,omega
T=4.5e-3
Ed=3.270
#Ed=3.441
Eu=3.149
h=0e-3
hz=0

dchi=1e-6+0j

omega=1e6
Tp=2*math.pi/omega

wL=1e10
wR=1e10


def eigofH(hz,h):
    delta=8e-3
    nom=delta+hz+math.sqrt((hz+delta)**2/4+h**2)
    dem=2*h
    norm=math.sqrt(nom**2+dem**2)
    [a,b]=[nom/norm,dem/norm]
    return [a,b]

a=eigofH(hz,h)[0]
b=eigofH(hz,h)[1]

print("a=",a,"b=",b)





def fLinU(muL,muR,Vg):
    f=1.0/(1.0+np.exp((Eu-Vg-muL)/T))
    return f

def fLoutU(muL,muR,Vg):
    f=1.0/(1.0+np.exp(-(Eu-Vg-muL)/T))
    return f

def fRinU(muL,muR,Vg):
    f=1.0/(1.0+np.exp((Eu-Vg-muR)/T))
    return f

def fRoutU(muL,muR,Vg):
    f=1.0/(1.0+np.exp(-(Eu-Vg-muR)/T))
    return f

def fLinD(muL,muR,Vg):
    f=1.0/(1.0+np.exp((Ed-Vg-muL)/T))
    return f

def fLoutD(muL,muR,Vg):
    f=1.0/(1.0+np.exp(-(Ed-Vg-muL)/T))
    return f

def fRinD(muL,muR,Vg):
    f=1.0/(1.0+np.exp((Ed-Vg-muR)/T))
    return f

def fRoutD(muL,muR,Vg):
    f=1.0/(1.0+np.exp(-(Ed-Vg-muR)/T))
    return f

    
def current(muL,muR,Vg):
    H=np.zeros((3,3),dtype=complex)




    Wp0L=wL*(fLinU(muL,muR,Vg)*a**2+fLinD(muL,muR,Vg)*b**2)
    Wp0R=wR*(fRinU(muL,muR,Vg)*a**2+fRinD(muL,muR,Vg)*b**2)
    Wm0L=wL*(fLinU(muL,muR,Vg)*b**2+fLinD(muL,muR,Vg)*a**2)
    Wm0R=wR*(fRinU(muL,muR,Vg)*b**2+fRinD(muL,muR,Vg)*a**2)
    
    W0pL=wL*(fLoutU(muL,muR,Vg)*a**2+fLoutD(muL,muR,Vg)*b**2)
    W0pR=wR*(fRoutU(muL,muR,Vg)*a**2+fRoutD(muL,muR,Vg)*b**2)
    W0mL=wL*(fLoutU(muL,muR,Vg)*b**2+fLoutD(muL,muR,Vg)*a**2)
    W0mR=wR*(fRoutU(muL,muR,Vg)*b**2+fRoutD(muL,muR,Vg)*a**2)
   
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
	    


def centerchange1(Vg):
#	global wL,wR
	

	Nplot=100
	biasrange=100e-3
	Bias=np.zeros(Nplot)
	Iup=np.zeros(Nplot)
	rho0=np.zeros(Nplot)
	rhop=np.zeros(Nplot)
	rhom=np.zeros(Nplot)
	muR=3.43
	
	for i in range(Nplot):
		Bias[i]=biasrange/Nplot*i
		Iup[i]=current(muR+Bias[i],muR)[0]*1.60217662*1e-7
		rho0[i]=current(muR+Bias[i],muR)[1]
		rhop[i]=current(muR+Bias[i],muR)[2]
		rhom[i]=current(muR+Bias[i],muR)[3]
	
	
	
	fig, ax1 = plt.subplots()
	ax2 = ax1.twinx()

	line1=ax1.plot(Bias*1e3,Iup,label=r'$I_\uparrow$',linewidth='2')
	ax1.set_ylim([0,120])

	line2=ax2.plot(Bias*1e3,rhop,label=r'$\rho_\uparrow$',linestyle='dashed',linewidth='2')
	line3=ax2.plot(Bias*1e3,rho0,label=r'$\rho_0$',linestyle='dashed',linewidth='2')
	line4=ax2.plot(Bias*1e3,rhom,label=r'$\rho_\downarrow$',linestyle='dashed',linewidth='2')
#	plt.axhline(0)
	ax1.set_xlabel('V (mV)')
	ax1.set_ylabel('I (pA)')
	ax2.set_ylabel(r'$\rho$')
	
	line=line1+line2+line3+line4
	labs = [l.get_label() for l in line]
	ax1.legend(line, labs, loc='best')

#	ax1.legend(loc='upper center')
#	ax2.legend(loc='best')
#	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	

	plt.show()
	return


def diamond():
	def fmt(x, pos):
# 	a, b = '{:.2e}'.format(x).split('e')
# 	b = int(b)
# 	return r'${} \times 10^{{{}}}$'.format(a, b)
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')
	Ngrid=201
	temp=0.2
	Lim=1.3

	xlist = np.linspace(Eu-temp, Eu+temp, Ngrid,endpoint=True)
	ylist = np.linspace(0, temp*2, Ngrid,endpoint=True)
	X, Y = np.meshgrid(xlist, ylist)
	Iup=np.zeros((Ngrid,Ngrid))
	Idn=np.zeros((Ngrid,Ngrid))
	Itot=np.zeros((Ngrid,Ngrid))
	dIup=np.zeros((Ngrid,Ngrid))
	dIdn=np.zeros((Ngrid,Ngrid))
	dItot=np.zeros((Ngrid,Ngrid))
#	muC=3.149
	
	for i in range(Ngrid):
		for j in range(Ngrid):
#			temp=current(muR+X[0][i],muR,Y[j][0])
			temp=current(Y[j][0]/2.0,-Y[j][0]/2.0,X[0][i])
#			temp=current(muR+X[0][i]/2,muR-X[0][i]/2,Y[j][0])
			Iup[j,i]=temp[0]*1.60217662*1e-10
			Idn[j,i]=temp[4]*1.60217662*1e-10
			Itot[j,i]=temp[5]*1.60217662*1e-10
	
	
#	print(current(3.48+0.02/2,3.48-0.02/2,0.02))
	
	for i in range(1,Ngrid-1):
		for j in range(Ngrid):
			dIup[j,i]=(Iup[j,i+1]-Iup[j,i-1])/(X[0][i+1]-X[0][i-1])
			dIdn[j,i]=(Idn[j,i+1]-Idn[j,i-1])/(X[0][i+1]-X[0][i-1])
#			dItot[j,i]=(Itot[j,i+1]-Itot[j,i-1])/(X[0][i+1]-X[0][i-1])
			dItot[j,i]=dIup[j,i]+dIdn[j,i]
	

	
	fig=plt.figure()
	Ncolor=16
	
	p1=fig.add_subplot(221)
	cp = p1.contourf(X,Y , Iup,np.linspace(0,Lim, Ncolor),cmap=plt.cm.Reds)
#	cp = p1.contourf(Y,X , Iup,np.linspace(0,120, 256),extend="both")
#	cp.cmap.set_under(color='blue')
#	cp.set_clim(0, 120)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))

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
#	plt.colorbar(cp, cax = cbaxes,format=ticker.FuncFormatter(fmt))  
	plt.colorbar(cp, cax = cbaxes)  

#	plt.gca().set_aspect('eVgual', adjustable='box')
	
	p3=fig.add_subplot(223)
	cp = p3.contourf( X,Y, Itot,np.linspace(0,Lim, Ncolor),cmap=plt.cm.Reds)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	plt.title('(c) Total')
	plt.xticks(np.arange(3.0,3.4,0.1))
	plt.xlabel('$V_g$ (V)')
	plt.ylabel('$V_b$ (V)')
	p3.yaxis.set_major_locator(MultipleLocator(0.1))

#	plt.gca().set_aspect('equal', adjustable='box')
	
#	p4=fig.add_subplot(224)
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
	
	
	
#	fig, ax1 = plt.subplots()
	ax1 = fig.add_subplot(224)
	ax2 = ax1.twinx()

	line1=ax1.plot(Bias,Iup,label=r'$I_\uparrow$',linewidth='2')
	ax1.set_ylim([0,1])
	ax2.set_ylim([0,1])

	line2=ax2.plot(Bias,rhop,label=r'$\rho_\uparrow$',linestyle='dashed',linewidth='2')
	line3=ax2.plot(Bias,rho0,label=r'$\rho_0$',linestyle='dashdot',linewidth='2')
	line4=ax2.plot(Bias,rhom,label=r'$\rho_\downarrow$',linestyle='dotted',linewidth='2')
#	plt.axhline(0)
	ax1.set_xlabel('$V_b$ (V)')
	ax1.set_ylabel('$I$ (nA)')
	ax2.set_ylabel(r'$\rho$')
	ax2.yaxis.set_major_locator(MultipleLocator(0.5))

#	ax1.set_xticklabels([0,0.2,0.4])
	plt.xticks([0,0.2,0.4])

	
	line=line1+line2+line3+line4
	labs = [l.get_label() for l in line]
#	ax1.legend(line, labs, ncol=2,loc='upper right')
	
	# Shrink current axis by 20%
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width * 1, box.height])
	ax2.set_position([box.x0, box.y0, box.width * 1, box.height])

# Put a legend to the right of the current axis
	ax1.legend(line,labs,ncol=1,loc='center left', bbox_to_anchor=(1.2, 0.5),frameon=False)
#	ax1.legend(line,labs,bbox_to_anchor=(0., 1.0, 1., .102), loc=3,
#           ncol=2, mode="expand", borderaxespad=0.)

	plt.title('(d) NDC')

#	plt.tight_layout()
	plt.show()
	
	

	return

#centerchange1(0)	
diamond()
print("time:", (time.time() - start_time)/60.0)
