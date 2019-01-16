#%matplotlib inline

import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from numpy import linalg as LA
from scipy import linalg as LAs
from scipy.linalg import eig as eigscipy
import time
start_time = time.time()


global T,Eu,Ed,h,hz, a, b ,dchi,Tp,omega,muL,muR
#T=5e-3
#Ed=3.498
#Eu=3.441

T=10e-3
#Ed=3.270
Eu=3.149
Ed=Eu
muL=Eu
muR=Eu

h=0e-3
hz=0

dchi=1e-6+0j

omega=1e6
Tp=2*math.pi/omega

#wL=1e9
#wR=1e9


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





def geometric(wLC,wRC,amp):
	def fLinU(muL,muR):
	    f=1.0/(1.0+np.exp((Eu+hz/2.0-muL)/T))
	    return f
	
	def fLoutU(muL,muR):
	    f=1.0-1.0/(1.0+np.exp((Eu+hz/2.0-muL)/T))
	    return f
	
	def fRinU(muL,muR):
	    f=1.0/(1.0+np.exp((Eu+hz/2.0-muR)/T))
	    return f
	
	def fRoutU(muL,muR):
	    f=1.0-1.0/(1.0+np.exp((Eu+hz/2.0-muR)/T))
	    return f
	
	def fLinD(muL,muR):
	    f=1.0/(1.0+np.exp((Ed-hz/2.0-muL)/T))
	    return f
	
	def fLoutD(muL,muR):
	    f=1.0-1.0/(1.0+np.exp((Ed-hz/2.0-muL)/T))
	    return f
	
	def fRinD(muL,muR):
	    f=1.0/(1.0+np.exp((Ed-hz/2.0-muR)/T))
	    return f
	
	def fRoutD(muL,muR):
	    f=1.0-1.0/(1.0+np.exp((Ed-hz/2.0-muR)/T))
	    return f
	
	    
	def Hamilt(wL,wR,chip,chim):
	    H=np.zeros((3,3),dtype=complex)


	
	

	
	    Wp0L=wL*(fLinU(muL,muR)*a**2+fLinD(muL,muR)*b**2)
	    Wp0R=wR*(fRinU(muL,muR)*a**2+fRinD(muL,muR)*b**2)
	    Wm0L=wL*(fLinU(muL,muR)*b**2+fLinD(muL,muR)*a**2)
	    Wm0R=wR*(fRinU(muL,muR)*b**2+fRinD(muL,muR)*a**2)
	    
	    W0pL=wL*(fLoutU(muL,muR)*a**2+fLoutD(muL,muR)*b**2)
	    W0pR=wR*(fRoutU(muL,muR)*a**2+fRoutD(muL,muR)*b**2)
	    W0mL=wL*(fLoutU(muL,muR)*b**2+fLoutD(muL,muR)*a**2)
	    W0mR=wR*(fRoutU(muL,muR)*b**2+fRoutD(muL,muR)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
	    e_chip=np.exp(1j*chip)
	    e_chip2=np.exp(-1j*chip)
	    e_chim=np.exp(1j*chim)
	    e_chim2=np.exp(-1j*chim)
	#    e_chi=1.0+1j*chi
	#    e_chi2=1.0-1j*chi
	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0pL+W0pR*e_chip
	    H[0,2]=W0mL+W0mR*e_chim
	    H[1,0]=Wp0L+Wp0R*e_chip2
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0L+Wm0R*e_chim2
	    H[2,1]=0
	    H[2,2]= -W0m
	    H=-H
	    
	#    print('0pL',W0pL)
	#    print('p0L',Wp0L)
	#    print('0pR',W0pR)
	#    print('p0R',Wp0R)
	#    print(H[0,1])
	#    print(H[1,0])
	    return H
	    

	
	    
	def HdiffL(wL,wR,chip,chim):
	    H=np.zeros((3,3),dtype=complex)
	
	    a=eigofH(hz,h)[0]
	    b=eigofH(hz,h)[1]
	
	    Wp0L=(fLinU(muL,muR)*a**2+fLinD(muL,muR)*b**2)
	    Wp0R=0
	    Wm0L=(fLinU(muL,muR)*b**2+fLinD(muL,muR)*a**2)
	    Wm0R=0
	    
	    W0pL=(fLoutU(muL,muR)*a**2+fLoutD(muL,muR)*b**2)
	    W0pR=0
	    W0mL=(fLoutU(muL,muR)*b**2+fLoutD(muL,muR)*a**2)
	    W0mR=0
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
	    e_chip=np.exp(1j*chip)
	    e_chip2=np.exp(-1j*chip)
	    e_chim=np.exp(1j*chim)
	    e_chim2=np.exp(-1j*chim)
	#    e_chi=1.0+1j*chi
	#    e_chi2=1.0-1j*chi
	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0pL+W0pR*e_chip
	    H[0,2]=W0mL+W0mR*e_chim
	    H[1,0]=Wp0L+Wp0R*e_chip2
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0L+Wm0R*e_chim2
	    H[2,1]=0
	    H[2,2]= -W0m
	    H=-H
	    return H
	
	def HdiffR(wL,wR,chip,chim):
	    H=np.zeros((3,3),dtype=complex)
	
	    a=eigofH(hz,h)[0]
	    b=eigofH(hz,h)[1]
	
	    Wp0L=0
	    Wp0R=(fRinU(muL,muR)*a**2+fRinD(muL,muR)*b**2)
	    Wp0Ru=(fRinU(muL,muR)*a**2)
	    Wp0Rd=(fRinD(muL,muR)*b**2)
	    Wm0L=0
	    Wm0R=(fRinU(muL,muR)*b**2+fRinD(muL,muR)*a**2)
	    Wm0Ru=(fRinU(muL,muR)*b**2)
	    Wm0Rd=(fRinD(muL,muR)*a**2)
	    
	    W0pL=0
	    W0pR=(fRoutU(muL,muR)*a**2+fRoutD(muL,muR)*b**2)
	    W0pRu=(fRoutU(muL,muR)*a**2)
	    W0pRd=(fRoutD(muL,muR)*b**2)
	    W0mL=0
	    W0mR=(fRoutU(muL,muR)*b**2+fRoutD(muL,muR)*a**2)
	    W0mRu=(fRoutU(muL,muR)*b**2)
	    W0mRd=(fRoutD(muL,muR)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
	    e_chip=np.exp(1j*chip)
	    e_chip2=np.exp(-1j*chip)
	    e_chim=np.exp(1j*chim)
	    e_chim2=np.exp(-1j*chim)
	#    e_chi=1.0+1j*chi
	#    e_chi2=1.0-1j*chi
	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0pL+W0pRu*e_chip+W0pRd*e_chim
	    H[0,2]=W0mL+W0mRu*e_chip+W0mRd*e_chim
	    H[1,0]=Wp0L+Wp0Ru*e_chip2+Wp0Rd*e_chim2
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0L+Wm0Ru*e_chip2+Wm0Rd*e_chim2
	    H[2,1]=0
	    H[2,2]= -W0m
	    H=-H
	    return H
	
	def F(wL,wR,chip,chim):
	#    H=Hamilt(t,chi)
	#    Hdagger=np.conj(H.transpose())
	#    print(H)
	#    print(Hdagger)
	
	    H=Hamilt(wL,wR,chip,chim)
	    A=eigscipy(H,left=True)
	    imin=A[0].argsort()[0]
	    imed=A[0].argsort()[1]
	    imax=A[0].argsort()[2]
	    eig0=A[0][imin]
	    eig1=A[0][imed]
	    eig2=A[0][imax]
#	    print(eig2)
	    phi0R=A[2].transpose()[imin]
	    phi0L=A[1].transpose()[imin]
	    phi1R=A[2].transpose()[imed]
	    phi1L=A[1].transpose()[imed]
	    phi2R=A[2].transpose()[imax]
	    phi2L=A[1].transpose()[imax]
	    
	    phi0L=np.conj(phi0L)
	    phi1L=np.conj(phi1L)
	    phi2L=np.conj(phi2L)
	    
	    phi0L=phi0L/np.matmul(phi0L,phi0R)
	    phi1L=phi1L/np.matmul(phi1L,phi1R)
	    phi2L=phi2L/np.matmul(phi2L,phi2R)

	#    phi1L=LA.eig(H1L_T)[1].transpose()[0]
	#    phi1L=np.conj(phi1L)
	
	#    print('eig of H1L',LA.eig(H1L))
	#    print('eig of H1L_T',LA.eig(H1L_T))
	#    print('left',phi1L)
	#    print('left from scipy',eigscipy(H1L,left=True))
	#    print('right',LA.eig(H1L)[1].transpose()[0])
	    
	    a=np.matmul(phi0L,HdiffL(wL,wR,chip,chim))
	    b1L=np.matmul(a,phi1R)
	    a=np.matmul(phi1L,HdiffR(wL,wR,chip,chim))
	    c1L=np.matmul(a,phi0R)
	    
	    a=np.matmul(phi0L,HdiffR(wL,wR,chip,chim))
	    b1R=np.matmul(a,phi1R)
	    a=np.matmul(phi1L,HdiffL(wL,wR,chip,chim))
	    c1R=np.matmul(a,phi0R)
	    
	    element1=(b1L*c1L-b1R*c1R)/((eig0-eig1))**2.0
	    
	    a=np.matmul(phi0L,HdiffL(wL,wR,chip,chim))
	    b2L=np.matmul(a,phi2R)
	    a=np.matmul(phi2L,HdiffR(wL,wR,chip,chim))
	    c2L=np.matmul(a,phi0R)
	    
	    a=np.matmul(phi0L,HdiffR(wL,wR,chip,chim))
	    b2R=np.matmul(a,phi2R)
	    a=np.matmul(phi2L,HdiffL(wL,wR,chip,chim))
	    c2R=np.matmul(a,phi0R)
	    
	    element2=(b2L*c2L-b2R*c2R)/((eig0-eig2))**2.0
	    
	    intg=element1+element2
	    return intg
	
	def BerryCurvature(wL,wR):
	#    BC=(-1j)*(integrand(muL,muR,dchi/2)-integrand(muL,muR,-dchi/2))/(dchi)
	    #BC=(-1j)*2*F(muL,muR,dchi/2)/(dchi)
	    temp1=F(wL,wR,dchi/2,0)
	    temp2=F(wL,wR,0,dchi/2)

	    
	    BC_U=2*(temp1.imag)/(dchi)
	    BC_D=2*(temp2.imag)/(dchi)
	#    print("BC",BC_C)
	    return [BC_U,BC_D]
	

#	print(BerryCurvature(Eu,Ed))	
#	print(BerryCurvature(Ed,Eu))
#	print(BerryCurvature(Eu,Eu))	
#	print(BerryCurvature(Ed,Ed))	

	
	def current():
	    Current_U=0.0
	    Current_D=0.0
	    Ngrid=41
	    #Ngrid=int(2*amp/0.5e-3)
	    xlist = np.linspace(muLC-amp, muLC+amp, Ngrid,endpoint=True)
	    ylist = np.linspace(muRC-amp, muRC+amp, Ngrid,endpoint=True)
	    X, Y = np.meshgrid(xlist, ylist)
	    #Z=np.zeros((Ngrid,Ngrid),dtype=complex)
	    for i in range(Ngrid):
	        for j in range(Ngrid):
	            if np.sqrt((X[0][i]-muLC)**2.0+(Y[j][0]-muRC)**2.0)<amp:
	                temp=BerryCurvature(X[0][i],Y[j][0])
	                Current_U=Current_U+temp[0]
	                Current_D=Current_D+temp[1]
	    Current_U=Current_U*(amp*2.0/Ngrid)**2.0*omega
	    Current_D=Current_D*(amp*2.0/Ngrid)**2.0*omega
	    return [Current_U,Current_D]
	

	def lineInt():
			Ngrid=2000
			
			phiL=np.zeros((Ngrid,3),dtype=complex)
			phiR=np.zeros((Ngrid,3),dtype=complex)
			chip=dchi
			chim=0.0
			for i in range(Ngrid):
				theta=i*2*np.pi/Ngrid
				muL=muLC+amp*np.cos(theta)
				muR=muRC+amp*np.sin(theta)
				H=Hamilt(muL,muR,chip,chim)
				A=eigscipy(H,left=True)
				imin=A[0].argsort()[0]
				#imed=A[0].argsort()[1]
				#imax=A[0].argsort()[2]
				eig0=A[0][imin]
				#eig1=A[0][imed]
				#eig2=A[0][imax]
				#print(eig2)
				phi0R=A[2].transpose()[imin]
				phi0L=A[1].transpose()[imin]
				#phi1R=A[2].transpose()[imed]
				#phi1L=A[1].transpose()[imed]
				#phi2R=A[2].transpose()[imax]
				#phi2L=A[1].transpose()[imax]
				
				phi0L=np.conj(phi0L)
				#phi1L=np.conj(phi1L)
				#phi2L=np.conj(phi2L)
				
				phi0L=phi0L/np.matmul(phi0L,phi0R)
				#phi1L=phi1L/np.matmul(phi1L,phi1R)
				#phi2L=phi2L/np.matmul(phi2L,phi2R)
				
				phiL[i,:]=phi0L
				phiR[i,:]=phi0R		
		
			int1=1.0
			for i in range(Ngrid-1):
				int1=int1*np.matmul(phiL[i],phiR[i+1])	
			int1=int1*np.matmul(phiL[Ngrid-1],phiR[0])
			int1=np.log(int1).imag/dchi*omega
			
			phiL=np.zeros((Ngrid,3),dtype=complex)
			phiR=np.zeros((Ngrid,3),dtype=complex)
			chip=0.0
			chim=dchi
			for i in range(Ngrid):
				theta=i*2*np.pi/Ngrid
				muL=muLC+amp*np.cos(theta)
				muR=muRC+amp*np.sin(theta)
				H=Hamilt(muL,muR,chip,chim)
				A=eigscipy(H,left=True)
				imin=A[0].argsort()[0]
				#imed=A[0].argsort()[1]
				#imax=A[0].argsort()[2]
				eig0=A[0][imin]
				#eig1=A[0][imed]
				#eig2=A[0][imax]
				#print(eig2)
				phi0R=A[2].transpose()[imin]
				phi0L=A[1].transpose()[imin]
				#phi1R=A[2].transpose()[imed]
				#phi1L=A[1].transpose()[imed]
				#phi2R=A[2].transpose()[imax]
				#phi2L=A[1].transpose()[imax]
				
				phi0L=np.conj(phi0L)
				#phi1L=np.conj(phi1L)
				#phi2L=np.conj(phi2L)
				
				phi0L=phi0L/np.matmul(phi0L,phi0R)
				#phi1L=phi1L/np.matmul(phi1L,phi1R)
				#phi2L=phi2L/np.matmul(phi2L,phi2R)
				
				phiL[i,:]=phi0L
				phiR[i,:]=phi0R		
		
			int2=1.0
			for i in range(Ngrid-1):
				#print(np.matmul(phiL[i],phiR[i+1]))
				int2=int2*np.matmul(phiL[i],phiR[i+1])	
			int2=int2*np.matmul(phiL[Ngrid-1],phiR[0])
			int2=np.log(int2).imag/dchi*omega			
			
			return [int1,int2]
			
#	print(lineInt())
		
		

	def fmt(x, pos):
  		a, b = '{:.2e}'.format(x).split('e')
  		b = int(b)
  		return r'${} \times 10^{{{}}}$'.format(a, b)
			return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')
   
	def BC_plot1():
	    Ngrid=11
	    xlist = np.linspace(1, 10, Ngrid,endpoint=True)
	    ylist = np.linspace(1, 10, Ngrid,endpoint=True)
	    X, Y = np.meshgrid(xlist, ylist)
	    BC_C=np.zeros((Ngrid,Ngrid),dtype=complex)
	    BC_S=np.zeros((Ngrid,Ngrid),dtype=complex)
	    BC_U=np.zeros((Ngrid,Ngrid),dtype=complex)
	    BC_D=np.zeros((Ngrid,Ngrid),dtype=complex)
	    for i in range(Ngrid):
	        for j in range(Ngrid):
	            temp=BerryCurvature(X[0][i],Y[j][0])
	            BC_C[i,j]=(temp[0]+temp[1])
	            BC_S[i,j]=(temp[0]-temp[1])
	            BC_U[i,j]=(temp[0])
	            BC_D[i,j]=(temp[1])

#### attention, the order of index

	            BC_C[j,i]=(temp[0]+temp[1])
	            BC_S[j,i]=(temp[0]-temp[1])
	            BC_U[j,i]=(temp[0])
	            BC_D[j,i]=(temp[1])
	            #print("pos=",X[0][i],Y[j][0])
	            #print("polarization=",BC_D[i,j]/BC_U[i,j])
	

	    fig=plt.figure()

	    p1=fig.add_subplot(221)
	    cp = p1.contourf(X, Y, BC_C,128)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(a) Charge')
	    plt.xlabel('$\mu_L$ (eV)')
	    plt.ylabel('$\mu_R$ (eV)')
	    plt.gca().set_aspect('equal', adjustable='box')
	    
	    p2=fig.add_subplot(222)
	    cp = p2.contourf(X, Y, BC_S,128)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(b) Spin')
	    plt.xlabel('$\mu_L$ (eV)')
	    plt.ylabel('$\mu_R$ (eV)')
	    plt.gca().set_aspect('equal', adjustable='box')

	    p3=fig.add_subplot(223)
	    cp = p3.contourf(X, Y, BC_U,128)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(c) Up')
	    plt.xlabel('$\mu_L$ (eV)')
	    plt.ylabel('$\mu_R$ (eV)')
	    plt.gca().set_aspect('equal', adjustable='box')

	    p4=fig.add_subplot(224)
	    cp = p4.contourf(X, Y, BC_D,128)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(d) Down')
	    plt.xlabel('$\mu_L$ (eV)')
	    plt.ylabel('$\mu_R$ (eV)')
	    plt.gca().set_aspect('equal', adjustable='box')

	    plt.tight_layout()
	    plt.show()
	    return

	BC_plot1()

#	return current()
	return lineInt()
	return

geometric(0,0,0)


#=====================================================
#=====================================================
#=====================================================
#=====================================================
#=====================================================
#=====================================================

def dynamic(muLC,muRC,amp):
#	muLC=8.5e-3
#	muRC=8e-3
#	amp=8e-3
	def muL(t):
	    m=muLC+amp*math.cos(omega*t)
	    return m
	
	def muR(t):
	    m=muRC+amp*math.sin(omega*t)
	    return m
	
	def fLinU(t):
	    f=1/(1+np.exp((Eu+hz/2-muL(t))/T))
	    return f
	
	def fLoutU(t):
	    f=1-1/(1+np.exp((Eu+hz/2-muL(t))/T))
	    return f
	
	def fRinU(t):
	    f=1/(1+np.exp((Eu+hz/2-muR(t))/T))
	    return f
	
	def fRoutU(t):
	    f=1-1/(1+np.exp((Eu+hz/2-muR(t))/T))
	    return f
	
	def fLinD(t):
	    f=1/(1+np.exp((Ed-hz/2-muL(t))/T))
	    return f
	
	def fLoutD(t):
	    f=1-1/(1+np.exp((Ed-hz/2-muL(t))/T))
	    return f
	
	def fRinD(t):
	    f=1/(1+np.exp((Ed-hz/2-muR(t))/T))
	    return f
	
	def fRoutD(t):
	    f=1-1/(1+np.exp((Ed-hz/2-muR(t))/T))
	    return f
	
	    
	def Hamilt(t,chip,chim):
	    H=np.zeros((3,3),dtype=complex)
	
	
	    Wp0L=wL*(fLinU(t)*a**2+fLinD(t)*b**2)
	    Wp0R=wR*(fRinU(t)*a**2+fRinD(t)*b**2)
	    Wm0L=wL*(fLinU(t)*b**2+fLinD(t)*a**2)
	    Wm0R=wR*(fRinU(t)*b**2+fRinD(t)*a**2)
	    
	    W0pL=wL*(fLoutU(t)*a**2+fLoutD(t)*b**2)
	    W0pR=wR*(fRoutU(t)*a**2+fRoutD(t)*b**2)
	    W0mL=wL*(fLoutU(t)*b**2+fLoutD(t)*a**2)
	    W0mR=wR*(fRoutU(t)*b**2+fRoutD(t)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
	    e_chip=np.exp(1j*chip)
	    e_chip2=np.exp(-1j*chip)
	    e_chim=np.exp(1j*chim)
	    e_chim2=np.exp(-1j*chim)
	#    e_chi=1.0+1j*chi
	#    e_chi2=1.0-1j*chi
	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0pL+W0pR*e_chip
	    H[0,2]=W0mL+W0mR*e_chim
	    H[1,0]=Wp0L+Wp0R*e_chip2
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0L+Wm0R*e_chim2
	    H[2,1]=0
	    H[2,2]= -W0m
	    H=-H
	    return H
	    
	def eig0(t,chip,chim):
			H=Hamilt(t,chip,chim)
			A=eigscipy(H)
			imin=A[0].argsort()[0]
			eig0=A[0][imin]
			return eig0
	
	def eig0chip(t):
	#		dchi=1e-6
			dertv=1j*(eig0(t,dchi/2,0)-eig0(t,-dchi/2,0))/dchi
			return dertv
			
	def eig0chim(t):
	#		dchi=1e-6
			dertv=1j*(eig0(t,0,dchi/2)-eig0(t,0,-dchi/2))/dchi
			return dertv
	
	def Idyn():
	    N=2000
	    CU=0
	    CD=0
	    for i in range(N):
	        CU=CU+eig0chip(Tp/N*i)/N
	        CD=CD+eig0chim(Tp/N*i)/N    #*Tp is included
	    CU=-CU
	    CD=-CD
	    return [CU,CD]
	

	    
	def steady(t):
	    Wp0L=wL*(fLinU(t)*a**2+fLinD(t)*b**2)
	    Wp0R=wR*(fRinU(t)*a**2+fRinD(t)*b**2)
	    Wm0L=wL*(fLinU(t)*b**2+fLinD(t)*a**2)
	    Wm0R=wR*(fRinU(t)*b**2+fRinD(t)*a**2)
	    
	    W0pL=wL*(fLoutU(t)*a**2+fLoutD(t)*b**2)
	    W0pR=wR*(fRoutU(t)*a**2+fRoutD(t)*b**2)
	    W0mL=wL*(fLoutU(t)*b**2+fLoutD(t)*a**2)
	    W0mR=wR*(fRoutU(t)*b**2+fRoutD(t)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
	    
	    eta=W0p*W0m+Wp0*W0m+W0p*Wm0
	    rho0=W0p*W0m/eta
	    rhop=Wp0*W0m/eta
	    rhom=W0p*Wm0/eta
	    
	    up=W0pR*rhop-Wp0R*rho0
	    dw=W0mR*rhom-Wm0R*rho0
	    return [up,dw]			
	
	#print('%8.2e' % currentp(),'%8.2e' % currentm())
	def check():
		N=10
		tt=np.zeros(N)
		eigU=np.zeros(N)
		eigD=np.zeros(N)
		styU=np.zeros(N)
		styD=np.zeros(N)		
		for i in range(N):
			tt[i]=Tp/N*i
#			print("t=",tt[i],":up:","eigen=",eig0chip(tt[i]),"steady=",steady(tt[i])[0])
#			print("t=",tt[i],":down:","eigen=",eig0chim(tt[i]),"steady=",steady(tt[i])[1])
			eigU[i]=eig0chip(tt[i]).real
			styU[i]=steady(tt[i])[0]
			eigD[i]=eig0chim(tt[i]).real
			styD[i]=steady(tt[i])[1]
		plt.plot(tt,eigU,"-")
		plt.plot(tt,styU,"<")
		plt.plot(tt,eigD,"r-")
		plt.plot(tt,styD,"r<")
		plt.show()
		return
#	check()
		
	return Idyn()
	
######################################################
######################################################
######################################################

def centerchange():
	global wL,wR
	fig=plt.figure()
	
	p1=fig.add_subplot(421)
	wL=1e9
	wR=1e9
	Nplot=1
	centrange=140e-3
	amp=80e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=centrange/Nplot*(i+1)
	    muRC[i]=muLC[i]
	    temp1=geometric(muLC[i],muRC[i],amp)
	    temp2=dynamic(muLC[i],muRC[i],amp)
	    CU[i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
	    CD[i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
	    print(CU[i],CD[i])
	p1.plot(muLC*1e3,CU+CD,'r',label='Charge')
	p1.plot(muLC*1e3,CU-CD,'b',label='Spin')
	p1.plot(muLC*1e3,CU,'r--',label='Up')
	p1.plot(muLC*1e3,CD,'b--',label='Down')
	plt.axhline(0)
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('I (fA)')
	plt.title('(a)')
#	plt.legend(loc='upper left')
#	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
 
	
	p2=fig.add_subplot(422)
	wL=1e9
	wR=1e9
#	Nplot=60
	centrange=140e-3
	amp=80e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=centrange/Nplot*(i+1)+0e-3
	    muRC[i]=muLC[i]+0.025e-3
	    temp1=geometric(muLC[i],muRC[i],amp)
	    temp2=dynamic(muLC[i],muRC[i],amp)
	    CU[i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
	    CD[i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
	    print(CU[i],CD[i])
	p2.plot(muLC*1e3,CU+CD,'r',label='Charge')
	p2.plot(muLC*1e3,CU-CD,'b',label='Spin')
	p2.plot(muLC*1e3,CU,'r--',label='Up')
	p2.plot(muLC*1e3,CD,'b--',label='Down')
	plt.axhline(0)
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('I (fA)')
	plt.title('(b)')

#	plt.legend(loc='upper left')
#	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
	p3=fig.add_subplot(423)
	wL=1.005e9
	wR=1.0e9
#	Nplot=60
	centrange=140e-3
	amp=80e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=centrange/Nplot*(i+1)+20e-3
	    muRC[i]=muLC[i]
	    temp1=geometric(muLC[i],muRC[i],amp)
	    temp2=dynamic(muLC[i],muRC[i],amp)
	    CU[i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
	    CD[i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
	    print(CU[i],CD[i])
	p3.plot(muLC*1e3,CU+CD,'r',label='Charge')
	p3.plot(muLC*1e3,CU-CD,'b',label='Spin')
	p3.plot(muLC*1e3,CU,'r--',label='Up')
	p3.plot(muLC*1e3,CD,'b--',label='Down')
	plt.axhline(0)
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('I (fA)')
	plt.title('(c)')
#	plt.legend(loc='upper left')
#	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
	p4=fig.add_subplot(424)
	wL=10e9
	wR=1e9
#	Nplot=60
	centrange=140e-3
	amp=80e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=centrange/Nplot*(i+1)
	    muRC[i]=muLC[i]-10e-3
	    temp1=geometric(muLC[i],muRC[i],amp)
	    temp2=dynamic(muLC[i],muRC[i],amp)
	    CU[i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
	    CD[i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
	    print(CU[i],CD[i])
#	p4.plot(muLC*1e3,CU+CD,'r',label='Charge')
#	p4.plot(muLC*1e3,CU-CD,'b',label='Spin')
#	p4.plot(muLC*1e3,CU,'r--',label='Up')
#	p4.plot(muLC*1e3,CD,'b--',label='Down')
	l1=p4.plot(muLC*1e3,CU+CD,'r')
	l2=p4.plot(muLC*1e3,CU-CD,'b')
	l3=p4.plot(muLC*1e3,CU,'r--')
	l4=p4.plot(muLC*1e3,CD,'b--')

	plt.axhline(0)
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('I (fA)')
	plt.title('(d)')
#	plt.legend(loc='upper left')
#	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 


	
	plt.tight_layout() 
#	plt.legend((l1, l2,l3,l4), ('Charge', 'Spin','Up','Down'), loc='ower center')

	plt.show()
	return


def centerchange1():
	global wL,wR
	fig=plt.figure()
	
	p1=fig.add_subplot(211)
	wL=1e9
	wR=1e9
	Nplot=51
	centrange=400e-3
	amp=160e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=3.0+centrange/(Nplot-1)*i
	    muRC[i]=muLC[i]
	    temp1=geometric(muLC[i],muRC[i],amp)
	    temp2=dynamic(muLC[i],muRC[i],amp)
	    CU[i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
	    CD[i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
	    print(CU[i],CD[i])
	p1.plot(muLC,CU+CD,'r',label='Charge')
	p1.plot(muLC,CU-CD,'b',label='Spin')
	p1.plot(muLC,CU,'r--',label='Up')
	p1.plot(muLC,CD,'b--',label='Down')
	plt.axhline(0)
	plt.xlabel('$\mu_L^{(C)}$ (eV)')
	plt.ylabel('I (fA)')
	plt.title('(a)')
#	plt.legend(loc='upper left')
#	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
 
	
	p2=fig.add_subplot(212)
	wL=1e9
	wR=1e9
#	Nplot=51
	centrange=200e-3
	amp=160e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=3.10+centrange/(Nplot-1)*(i)
	    muRC[i]=muLC[i]+0.07e-3
	    temp1=geometric(muLC[i],muRC[i],amp)
	    temp2=dynamic(muLC[i],muRC[i],amp)
	    CU[i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
	    CD[i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
	    print(CU[i],CD[i])
	p2.plot(muLC,CU+CD,'r',label='Charge')
	p2.plot(muLC,CU-CD,'b',label='Spin')
	p2.plot(muLC,CU,'r--',label='Up')
	p2.plot(muLC,CD,'b--',label='Down')
	plt.axhline(0)
	plt.xlabel('$\mu_L^{(C)}$ (eV)')
	plt.ylabel('I (fA)')
	plt.title('(b)')




	
	plt.tight_layout() 
#	plt.legend((l1, l2,l3,l4), ('Charge', 'Spin','Up','Down'), loc='ower center')

	plt.show()
	return



def contour():
	global wL,wR
	fig=plt.figure()
	wL=10e9
	wR=1e9
	Ngrid=17

	amp=80e-3
	xlist = np.linspace(Eu-15e-3, Ed+15e-3, Ngrid,endpoint=True)
	ylist = np.linspace(Eu-15e-3, Ed+15e-3, Ngrid,endpoint=True)
	X, Y = np.meshgrid(xlist, ylist)
	C_C=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_S=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_U=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_D=np.zeros((Ngrid,Ngrid),dtype=complex)
	for i in range(Ngrid):
		for j in range(Ngrid):
			temp1=geometric(X[0][i],Y[j][0],amp)
#			temp1=np.array([0,0])
			temp2=dynamic(X[0][i],Y[j][0],amp)
			C_U[j,i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
			C_D[j,i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
			C_C[j,i]=C_U[j,i]+C_D[j,i]
			C_S[j,i]=C_U[j,i]-C_D[j,i]



	amp=0e-3
#	Ngrid=3
	xlist = np.linspace(Eu-15e-3, Ed+15e-3, Ngrid,endpoint=True)
	ylist = np.linspace(Eu-15e-3, Ed+15e-3, Ngrid,endpoint=True)
	X, Y = np.meshgrid(xlist, ylist)
	C_C2=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_S2=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_U2=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_D2=np.zeros((Ngrid,Ngrid),dtype=complex)
	for i in range(Ngrid):
		for j in range(Ngrid):
			temp1=geometric(X[0][i],Y[j][0],amp)
#			temp1=np.array([0,0])
			temp2=dynamic(X[0][i],Y[j][0],amp)
			C_U2[j,i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-4
			C_D2[j,i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-4
			C_C2[j,i]=C_U2[j,i]+C_D2[j,i]
			C_S2[j,i]=C_U2[j,i]-C_D2[j,i]
			
	matplotlib.rcParams['contour.negative_linestyle'] = 'solid'    
	p5=fig.add_subplot(221)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p5.contour(X*1e3, Y*1e3, C_C,colors='k')
#	cp.collections[2].set_color('red')
	cp.collections[2].set_linewidth(3.5)  
	plt.clabel(cp,fmt='%4.0f')
	cp = p5.contour(X*1e3, Y*1e3, C_C2,colors='red')
#	cp.collections[2].set_color('red')
	cp.collections[2].set_linewidth(3.5)  
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('$\mu_R^{(C)}$ (meV)')
	plt.title('(a) Charge')

	p6=fig.add_subplot(222)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p6.contour(X*1e3, Y*1e3, C_S,colors='k')
	plt.clabel(cp,fmt='%4.0f')
	cp = p6.contour(X*1e3, Y*1e3, C_S2,colors='red')
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('$\mu_R^{(C)}$ (meV)')
	plt.title('(b) Spin')

	p7=fig.add_subplot(223)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p7.contour(X*1e3, Y*1e3, C_U,colors='k')
#	cp.collections[2].set_color('red')
	cp.collections[2].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	cp = p7.contour(X*1e3, Y*1e3, C_U2,colors='red')
#	cp.collections[2].set_color('red')
	cp.collections[2].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('$\mu_R^{(C)}$ (meV)')
	plt.title('(c) Up')

	p8=fig.add_subplot(224)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p8.contour(X*1e3, Y*1e3, C_D,colors='k')
#	cp.collections[2].set_color('red')
	cp.collections[2].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	cp = p8.contour(X*1e3, Y*1e3, C_D2,colors='red')
#	cp.collections[2].set_color('red')
	cp.collections[1].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}$ (meV)')
	plt.ylabel('$\mu_R^{(C)}$ (meV)')
	plt.title('(d) Down')
	
	plt.gca().set_aspect('equal', adjustable='box')
	plt.tight_layout() 
#	plt.legend((l1, l2,l3,l4), ('Charge', 'Spin','Up','Down'), loc='ower center')

	plt.show()

	return

def contour1():
	global wL,wR
	fig=plt.figure()
	wL=10e9
	wR=1e9
	Ngrid=33


	xlist = np.linspace(-10e-3, 150e-3, Ngrid,endpoint=True)
	ylist = np.linspace(0, 160e-3, Ngrid,endpoint=True)
	print(xlist,ylist)
	X, Y = np.meshgrid(xlist, ylist)
	C_C=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_S=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_U=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_D=np.zeros((Ngrid,Ngrid),dtype=complex)
	for i in range(Ngrid):
		for j in range(Ngrid):
			temp1=geometric(3.13+X[0][i],3.13,Y[j][0])
#			temp1=np.array([0,0])
			temp2=dynamic(3.13+X[0][i],3.13,Y[j][0])
			C_U[j,i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-7
			C_D[j,i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-7
			C_C[j,i]=C_U[j,i]+C_D[j,i]
			C_S[j,i]=C_U[j,i]-C_D[j,i]




	xlist = np.linspace(-10e-3, 150e-3, Ngrid,endpoint=True)
	ylist = np.linspace(0, 160e-3, Ngrid,endpoint=True)
	X, Y = np.meshgrid(xlist, ylist)
	C_C2=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_S2=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_U2=np.zeros((Ngrid,Ngrid),dtype=complex)
	C_D2=np.zeros((Ngrid,Ngrid),dtype=complex)
	for i in range(Ngrid):
		for j in range(Ngrid):
			temp1=geometric(3.19+X[0][i],3.19,Y[j][0])
#			temp1=np.array([0,0])
			temp2=dynamic(3.19+X[0][i],3.19,Y[j][0])
			C_U2[j,i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-7
			C_D2[j,i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-7
			C_C2[j,i]=C_U2[j,i]+C_D2[j,i]
			C_S2[j,i]=C_U2[j,i]-C_D2[j,i]
			
	matplotlib.rcParams['contour.negative_linestyle'] = 'solid'    
	p5=fig.add_subplot(221)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p5.contour(X*1e3, Y*1e3, C_C,colors='k')
#	cp.collections[2].set_color('red')
	cp.collections[6].set_linewidth(3.5)  
	plt.clabel(cp,fmt='%4.0f')
	cp = p5.contour(X*1e3, Y*1e3, C_C2,colors='red')
#	cp.collections[2].set_color('red')
	cp.collections[4].set_linewidth(3.5)  
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
	plt.ylabel('$A$ (meV)')
	plt.title('(a) Charge')

	p6=fig.add_subplot(222)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p6.contour(X*1e3, Y*1e3, C_S,colors='k')
	plt.clabel(cp,fmt='%4.0f')
	cp = p6.contour(X*1e3, Y*1e3, C_S2,colors='red')
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
	plt.ylabel('$A$ (meV)')
	plt.title('(b) Spin')

	p7=fig.add_subplot(223)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p7.contour(X*1e3, Y*1e3, C_U,colors='k')
#	cp.collections[2].set_color('red')
	cp.collections[6].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	cp = p7.contour(X*1e3, Y*1e3, C_U2,colors='red')
#	cp.collections[2].set_color('red')
	cp.collections[5].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
	plt.ylabel('$A$ (meV)')
	plt.title('(c) Up')

	p8=fig.add_subplot(224)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = plt.contourf(X*1e3, Y*1e3, C_C,128)
#	plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	cp = p8.contour(X*1e3, Y*1e3, C_D,colors='k')
#	cp.collections[2].set_color('red')
	cp.collections[4].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	cp = p8.contour(X*1e3, Y*1e3, C_D2,colors='red')
#	cp.collections[2].set_color('red')
	cp.collections[4].set_linewidth(3.5) 
	plt.clabel(cp,fmt='%4.0f')
	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
	plt.ylabel('$A$ (meV)')
	plt.title('(d) Down')
	
#	plt.gca().set_aspect('equal', adjustable='box')
	plt.tight_layout() 
#	plt.legend((l1, l2,l3,l4), ('Charge', 'Spin','Up','Down'), loc='ower center')

	plt.show()

	return



def geometricplot():
    Nplot=2
    centrange=30e-3
    amp=8e-3
    muLC=np.zeros(Nplot)
    muRC=np.zeros(Nplot)
    CU=np.zeros(Nplot)
    CD=np.zeros(Nplot)
    for i in range(Nplot):
        muLC[i]=centrange/Nplot*(i+1)+0e-3
        muRC[i]=muLC[i]+0.0e-3
        temp1=geometric(muLC[i],muRC[i],amp)
#        temp2=[0.0,0,0]
        CU[i]=(temp1[0])/2.0/np.pi
        CD[i]=(temp1[1])/2.0/np.pi
        print("up",CU[i])
        print("dn",CD[i])

    plt.plot(muLC*1e3,CU+CD,'r',label='Charge')
    plt.plot(muLC*1e3,CU-CD,'b',label='Spin')
    plt.plot(muLC*1e3,CU,'r--',label='Up')
    plt.plot(muLC*1e3,CD,'b--',label='Down')
    plt.axhline(0)
    plt.xlabel('Center of chemical pot (meV)')
    plt.ylabel('Current')
    plt.legend(loc='upper left')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.show()
    


    return

#geometricplot()

#centerchange()	

centerchange1()	

wL=1e9
wR=1e9
muL=3.43
muR=3.43
#print(dynamic(muL,muR,0e-3)[0]*1.60217662*1e-4,dynamic(muL,muR,0e-3)[1]*1.60217662*1e-4)
print(geometric(muL,muR,0e-3))


#contour1()




print("time:", (time.time() - start_time)/60.0)



#x=[0,12,12,20,20,30]
#y=[0,0, 1e5*0.94,1e5*0.94,0,0]
#y2=[0,0, -1e5*0.94,-1e5*0.94,0,0]
#plt.plot(x,y)
#plt.plot(x,y2,'b--')
#
#plt.axhline(0)
#plt.xlabel('Center of chemical pot (mV)')
#plt.ylabel('Current')
#plt.legend(loc='upper left')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.show()


######################################################
######################################################
######################################################

# delta/temperature
