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

label_size = 17
#matplotlib.rcParams['xtick.labelsize'] = label_size 
#matplotlib.rcParams['ytick.labelsize'] = label_size 
matplotlib.rcParams.update({'font.size': 16})
#matplotlib.rcParams['xtick.major.pad']='0'

global T,Eu,Ed,wL,wR, dchi,Tp,omega

shift=0.0
scale=1e0
T=10e-3 # temperature
Ed=scale*(3.270-3.149+shift)
Eu=scale*(0.0+shift)


dchi=8e-5+0j
#print(dchi)

omega=1e9  # frequency
Tp=2*math.pi/omega   # period

wL=1e10
wR=1e10


def geometric(muLC,muRC,amp):
	def fLinU(muL,muR):
	    f=1.0/(1.0+np.exp((Eu-muL)/T))
	    return f
	
	def fLoutU(muL,muR):
	    f=1.0/(1.0+np.exp(-(Eu-muL)/T))
	    return f
	
	def fRinU(muL,muR):
	    f=1.0/(1.0+np.exp((Eu-muR)/T))
	    return f
	
	def fRoutU(muL,muR):
	    f=1.0/(1.0+np.exp(-(Eu-muR)/T))
	    return f
	
	def fLinD(muL,muR):
	    f=1.0/(1.0+np.exp((Ed-muL)/T))
	    return f
	
	def fLoutD(muL,muR):
	    f=1.0/(1.0+np.exp(-(Ed-muL)/T))
	    return f
	
	def fRinD(muL,muR):
	    f=1.0/(1.0+np.exp((Ed-muR)/T))
	    return f
	
	def fRoutD(muL,muR):
	    f=1.0/(1.0+np.exp(-(Ed-muR)/T))
	    return f
	
	    
	def Hamilt(muL,muR,chip,chim):
	    H=np.zeros((3,3),dtype=np.dtype(np.complex128))
	
	

	
	    Wp0L=wL*(fLinU(muL,muR))
	    Wp0R=wR*(fRinU(muL,muR))
	    Wm0L=wL*(fLinD(muL,muR))
	    Wm0R=wR*(fRinD(muL,muR))
	    
	    W0pL=wL*(fLoutU(muL,muR))
	    W0pR=wR*(fRoutU(muL,muR))
	    W0mL=wL*(fLoutD(muL,muR))
	    W0mR=wR*(fRoutD(muL,muR))
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
	    e_chip=np.exp(1j*chip)
	    e_chip2=np.exp(-1j*chip)
	    e_chim=np.exp(1j*chim)
	    e_chim2=np.exp(-1j*chim)



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
	    
	def fLinUdiff(muL,muR):
	    a=np.exp((Eu-muL)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=1.0/(1.0/a+2+a)/T
	    return f
	
	def fLoutUdiff(muL,muR):
	    a=np.exp((Eu-muL)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=-1.0/(1.0/a+2+a)/T
	    return f
	
	def fRinUdiff(muL,muR):
	    a=np.exp((Eu-muR)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=1.0/(1.0/a+2+a)/T
	    return f
	
	def fRoutUdiff(muL,muR):
	    a=np.exp((Eu-muR)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=-1.0/(1.0/a+2+a)/T
	    return f
	
	def fLinDdiff(muL,muR):
	    a=np.exp((Ed-muL)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=1.0/(1.0/a+2+a)/T
	    return f
	
	def fLoutDdiff(muL,muR):
	    a=np.exp((Ed-muL)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=-1.0/(1.0/a+2+a)/T
	    return f
	
	def fRinDdiff(muL,muR):
	    a=np.exp((Ed-muR)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=1.0/(1.0/a+2+a)/T
	    return f
	
	def fRoutDdiff(muL,muR):
	    a=np.exp((Ed-muR)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=-1.0/(1.0/a+2+a)/T
	    return f
	
	    
	def HdiffL(muL,muR,chip,chim):
	    H=np.zeros((3,3),dtype=np.dtype(np.complex128))
	

	
	    Wp0L=wL*(fLinUdiff(muL,muR))
	    Wp0R=0
	    Wm0L=wL*(fLinDdiff(muL,muR))
	    Wm0R=0
	    
	    W0pL=wL*(fLoutUdiff(muL,muR))
	    W0pR=0
	    W0mL=wL*(fLoutDdiff(muL,muR))
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
	
	def HdiffR(muL,muR,chip,chim):
	    H=np.zeros((3,3),dtype=np.dtype(np.complex128))
	
	
	    Wp0L=0
	    Wp0R=wR*(fRinUdiff(muL,muR))
	    Wp0Ru=wR*(fRinUdiff(muL,muR))
	    Wp0Rd=0
	    Wm0L=0
	    Wm0R=wR*(fRinDdiff(muL,muR))
	    Wm0Ru=0
	    Wm0Rd=wR*(fRinDdiff(muL,muR))
	    
	    W0pL=0
	    W0pR=wR*(fRoutUdiff(muL,muR))
	    W0pRu=wR*(fRoutUdiff(muL,muR))
	    W0pRd=0
	    W0mL=0
	    W0mR=wR*(fRoutDdiff(muL,muR))
	    W0mRu=0
	    W0mRd=wR*(fRoutDdiff(muL,muR))
	   
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
	
	def F(muL,muR,chip,chim):
	#    H=Hamilt(t,chi)
	#    Hdagger=np.conj(H.transpose())
	#    print(H)
	#    print(Hdagger)
	
	    H=Hamilt(muL,muR,chip,chim)
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
	    
	    a=np.matmul(phi0L,HdiffL(muL,muR,chip,chim))
	    b1L=np.matmul(a,phi1R)
	    a=np.matmul(phi1L,HdiffR(muL,muR,chip,chim))
	    c1L=np.matmul(a,phi0R)
	    
	    a=np.matmul(phi0L,HdiffR(muL,muR,chip,chim))
	    b1R=np.matmul(a,phi1R)
	    a=np.matmul(phi1L,HdiffL(muL,muR,chip,chim))
	    c1R=np.matmul(a,phi0R)
	    
	    element1=(b1L*c1L-b1R*c1R)/((eig0-eig1))**2.0
	    
	    a=np.matmul(phi0L,HdiffL(muL,muR,chip,chim))
	    b2L=np.matmul(a,phi2R)
	    a=np.matmul(phi2L,HdiffR(muL,muR,chip,chim))
	    c2L=np.matmul(a,phi0R)
	    
	    a=np.matmul(phi0L,HdiffR(muL,muR,chip,chim))
	    b2R=np.matmul(a,phi2R)
	    a=np.matmul(phi2L,HdiffL(muL,muR,chip,chim))
	    c2R=np.matmul(a,phi0R)
	    
	    element2=(b2L*c2L-b2R*c2R)/((eig0-eig2))**2.0
	    
	    intg=element1+element2
	    return intg
	
	def BerryCurvature(muL,muR):

#	    temp1=F(muL,muR,dchi/2,0)
#	    temp2=F(muL,muR,0,dchi/2)
#
#	    
#	    BC_U=2*(temp1.imag)/(dchi)
#	    BC_D=2*(temp2.imag)/(dchi)
	    
	    temp1=8.0*F(muL,muR,dchi/2.0,0)-F(muL,muR,dchi,0)
	    temp2=8.0*F(muL,muR,0,dchi/2.0)-F(muL,muR,0,dchi)
	    BC_U=(temp1.imag)/(3.0*dchi)
	    BC_D=(temp2.imag)/(3.0*dchi)

#	    temp1=8*F(muL,muR,dchi/2,0)-8*F(muL,muR,-dchi/2,0)-F(muL,muR,dchi,0)+F(muL,muR,-dchi,0)
#	    temp2=8*F(muL,muR,0,dchi/2)-8*F(muL,muR,0,-dchi/2)-F(muL,muR,0,dchi)+F(muL,muR,0,-dchi)

	    
#	    BC_U=(temp1.imag)/(6*dchi)
#	    BC_D=(temp2.imag)/(6*dchi)


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
	    #Z=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
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
			
			phiL=np.zeros((Ngrid,3),dtype=np.dtype(np.complex128))
			phiR=np.zeros((Ngrid,3),dtype=np.dtype(np.complex128))
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
			
			phiL=np.zeros((Ngrid,3),dtype=np.dtype(np.complex128))
			phiR=np.zeros((Ngrid,3),dtype=np.dtype(np.complex128))
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
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')
    
	def BC_plot():
	    Ngrid=41 #41
	    xlist = np.linspace(Eu-50e-3, Ed+50e-3, Ngrid,endpoint=True)
	    ylist = np.linspace(Eu-50e-3, Ed+50e-3, Ngrid,endpoint=True)
	    X, Y = np.meshgrid(xlist, ylist)
	    BC_C=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	    BC_S=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	    BC_U=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	    BC_D=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	    for i in range(Ngrid):
	        for j in range(Ngrid):
	            temp=BerryCurvature(X[0][i],Y[j][0])
#	            BC_C[i,j]=(temp[0]+temp[1])
#	            BC_S[i,j]=(temp[0]-temp[1])
#	            BC_U[i,j]=(temp[0])
#	            BC_D[i,j]=(temp[1])

##### attention, the order of index

	            BC_C[j,i]=(temp[0]+temp[1])
	            BC_S[j,i]=(temp[0]-temp[1])
	            BC_U[j,i]=(temp[0])

	            if i>j:
	            	BC_U[j,i]=(BC_U[j,i]+BC_U[i,j])/2.0
	            	BC_U[i,j]=BC_U[j,i]
	            BC_D[j,i]=(temp[1])

	    fig=plt.figure()

#	    p1=fig.add_subplot(221)
#	    cp = p1.contourf(X, Y, BC_C,128)
#	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
#	    plt.title('(a) Charge')
#	    plt.xlabel('$\mu_L$ (eV)')
#	    plt.ylabel('$\mu_R$ (eV)')
#	    plt.gca().set_aspect('equal', adjustable='box')
#	    
#	    p2=fig.add_subplot(222)
#	    cp = p2.contourf(X, Y, BC_S,128)
#	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
#	    plt.title('(b) Spin')
#	    plt.xlabel('$\mu_L$ (eV)')
#	    plt.ylabel('$\mu_R$ (eV)')
#	    plt.gca().set_aspect('equal', adjustable='box')

	    p3=fig.add_subplot(121)
	    cp = p3.contourf(X*1e3, Y*1e3, BC_U,16)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(a) Up')
	    plt.xlabel('$\mu_L$ (meV)',fontsize=19)
	    plt.ylabel('$\mu_R$ (meV)',fontsize=19)
#	    plt.gca().set_aspect('equal', adjustable='box')

	    p4=fig.add_subplot(122)
	    cp = p4.contourf(X*1e3, Y*1e3, BC_D,16)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(b) Down')
	    plt.xlabel('$\mu_L$ (meV)',fontsize=19)
	    plt.ylabel('$\mu_R$ (meV)',fontsize=19)
#	    plt.gca().set_aspect('equal', adjustable='box')

	    plt.tight_layout()
	    plt.show()
	    return

	BC_plot()

#	return current()
	return lineInt()


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
	    f=1/(1+np.exp((Eu-muL(t))/T))
	    return f
	
	def fLoutU(t):
	    f=1-1/(1+np.exp((Eu-muL(t))/T))
	    return f
	
	def fRinU(t):
	    f=1/(1+np.exp((Eu-muR(t))/T))
	    return f
	
	def fRoutU(t):
	    f=1-1/(1+np.exp((Eu-muR(t))/T))
	    return f
	
	def fLinD(t):
	    f=1/(1+np.exp((Ed-muL(t))/T))
	    return f
	
	def fLoutD(t):
	    f=1-1/(1+np.exp((Ed-muL(t))/T))
	    return f
	
	def fRinD(t):
	    f=1/(1+np.exp((Ed-muR(t))/T))
	    return f
	
	def fRoutD(t):
	    f=1-1/(1+np.exp((Ed-muR(t))/T))
	    return f
	
	    
	def Hamilt(t,chip,chim):
	    H=np.zeros((3,3),dtype=np.dtype(np.complex128))
	
	
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
#	Nplot=51
	Nplot=2
	centrange=350e-3
	amp=160e-3
	muLC=np.zeros(Nplot)
	muRC=np.zeros(Nplot)
	CU=np.zeros(Nplot)
	CD=np.zeros(Nplot)
	for i in range(Nplot):
	    muLC[i]=-0.1+centrange/(Nplot-1)*i
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
	plt.legend(loc='upper left')
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
	    muLC[i]=-0.05+centrange/(Nplot-1)*(i)
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
	C_C=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_S=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_U=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_D=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
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
	C_C2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_S2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_U2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_D2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
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
	wL=1e10
	wR=1e10
#	Ngrid=3
	Ngrid=21


	xlist = np.linspace(-20e-3, 10e-3, Ngrid,endpoint=True)
	ylist = np.linspace(50e-3, 200e-3, Ngrid,endpoint=True)
	print(xlist,ylist)
	X, Y = np.meshgrid(xlist, ylist)
	C_C=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_S=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_U=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
	C_D=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))

	Z=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))

	for i in range(Ngrid):
		for j in range(Ngrid):
#			muRC=3.13-3.149
			muRC=0
			temp1=geometric(muRC+X[0][i],muRC,Y[j][0])
#			temp1=np.array([0,0])
			temp2=dynamic(muRC+X[0][i],muRC,Y[j][0])
			C_U[j,i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-7
			C_D[j,i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-7
			C_C[j,i]=C_U[j,i]+C_D[j,i]
			C_S[j,i]=C_U[j,i]-C_D[j,i]
			Z[j,i]=0
			print('i=',i)




#	xlist = np.linspace(-10e-3, 150e-3, Ngrid,endpoint=True)
#	ylist = np.linspace(0, 160e-3, Ngrid,endpoint=True)
#	X, Y = np.meshgrid(xlist, ylist)
#	C_C2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
#	C_S2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
#	C_U2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
#	C_D2=np.zeros((Ngrid,Ngrid),dtype=np.dtype(np.complex128))
#	for i in range(Ngrid):
#		for j in range(Ngrid):
#			muRC=3.19-3.149
#			temp1=geometric(muRC+X[0][i],muRC,Y[j][0])
#			temp2=dynamic(muRC+X[0][i],muRC,Y[j][0])
#			C_U2[j,i]=(temp1[0]/2.0/np.pi+temp2[0])*1.60217662*1e-7
#			C_D2[j,i]=(temp1[1]/2.0/np.pi+temp2[1])*1.60217662*1e-7
#			C_C2[j,i]=C_U2[j,i]+C_D2[j,i]
#			C_S2[j,i]=C_U2[j,i]-C_D2[j,i]
#			print('i=',i)
			
	matplotlib.rcParams['contour.negative_linestyle'] = 'solid'    
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	

#	p5=fig.add_subplot(221)
#	cp = p5.contour(X*1e3, Y*1e3, C_C,colors='k')
#	cp.collections[6].set_linewidth(3.5)  
#	plt.clabel(cp,fmt='%4.0f')
#	cp = p5.contour(X*1e3, Y*1e3, C_C2,colors='red')
#	cp.collections[4].set_linewidth(3.5)  
#	plt.clabel(cp,fmt='%4.0f')
#	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
#	plt.ylabel('$A$ (meV)')
#	plt.title('(a) Charge')

	p5=fig.add_subplot(221,projection='3d')
	p5.plot_surface(X*1e3, Y*1e3, C_C,rstride=1, cstride=1,linewidth=0.5,  alpha=0.3,cmap=cm.coolwarm)
#	p5.plot_surface(X*1e3, Y*1e3, Z,rstride=1, cstride=1,linewidth=0.5,  alpha=0.3,cmap="winter")
	cset = p5.contour(X*1e3, Y*1e3, C_C, [0],colors='k', zdir='z', offset=0)
	cset.collections[0].set_linewidth(2)  
	plt.clabel(cset,fmt='%4.0f')
	p5.set_xlabel('\n'+r'$e \bar{V}_b$ (meV)',linespacing=1)
	p5.set_ylabel('\n $A$ (meV)',linespacing=1)
	p5.zaxis.set_rotate_label(False)
	p5.set_zlabel('$I_C$ (pA)',rotation=90)
	p5.set_zlim(-80,120)
	p5.xaxis.set_major_locator(MultipleLocator(10))
	p5.yaxis.set_major_locator(MultipleLocator(50))
	p5.zaxis.set_major_locator(MultipleLocator(60))
	p5.view_init(elev=40, azim = 50)
	plt.title('(a) Charge')


#	p6=fig.add_subplot(222)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = p6.contour(X*1e3, Y*1e3, C_S,colors='k')
#	plt.clabel(cp,fmt='%4.0f')
#	cp = p6.contour(X*1e3, Y*1e3, C_S2,colors='red')
#	plt.clabel(cp,fmt='%4.0f')
#	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
#	plt.ylabel('$A$ (meV)')
#	plt.title('(b) Spin')
	p6=fig.add_subplot(222,projection='3d')
	p6.plot_surface(X*1e3, Y*1e3, C_S, rstride=1, cstride=1,linewidth=0.5,  alpha=0.3,cmap=cm.coolwarm)
	cset = p6.contour(X*1e3, Y*1e3, C_S, [0],colors='k', zdir='z', offset=0)
	cset.collections[0].set_linewidth(2)  
	p6.set_xlabel('\n'+r'$e \bar{V}_b$ (meV)',linespacing=1)
	p6.set_ylabel('\n $A$ (meV)',linespacing=1)
	p6.zaxis.set_rotate_label(False)
	p6.set_zlabel('$I_S$ (pA)',rotation=90)
	p6.set_zlim(-80,120)
	p6.xaxis.set_major_locator(MultipleLocator(10))
	p6.yaxis.set_major_locator(MultipleLocator(50))
	p6.zaxis.set_major_locator(MultipleLocator(60))
	p6.view_init(elev=40, azim = 50)
	plt.title('(b) Spin')

#	p7=fig.add_subplot(223)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = p7.contour(X*1e3, Y*1e3, C_U,colors='k')
#	plt.clabel(cp,fmt='%4.0f')
#	cp = p7.contour(X*1e3, Y*1e3, C_U2,colors='red')
#	plt.clabel(cp,fmt='%4.0f')
#	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
#	plt.ylabel('$A$ (meV)')
#	plt.title('(c) Up')
	p7=fig.add_subplot(223,projection='3d')
	p7.plot_surface(X*1e3, Y*1e3, C_U, rstride=1, cstride=1,linewidth=0.5, alpha=0.3,cmap=cm.coolwarm)
	cset = p7.contour(X*1e3, Y*1e3, C_U,[0], colors='k', zdir='z', offset=0)
	cset.collections[0].set_linewidth(2)  
	p7.set_xlabel('\n'+r'$e \bar{V}_b$ (meV)',linespacing=1)
	p7.set_ylabel('\n $A$ (meV)',linespacing=1)
	p7.zaxis.set_rotate_label(False)
	p7.set_zlabel(r'$I_\uparrow$ (pA)',rotation=90)
	p7.set_zlim(-80,120)
	p7.xaxis.set_major_locator(MultipleLocator(10))
	p7.yaxis.set_major_locator(MultipleLocator(50))
	p7.zaxis.set_major_locator(MultipleLocator(60))
	p7.view_init(elev=40, azim = 50)
	plt.title('\n (c) Up',linespacing=1)
	
	def fmt(x, pos):                                                            
		return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')	
#	cp = p8.contour(X*1e3, Y*1e3, C_D,colors='k')
#	plt.clabel(cp,fmt='%4.0f')
#	cp = p8.contour(X*1e3, Y*1e3, C_D2,colors='red')
#	plt.clabel(cp,fmt='%4.0f')
#	plt.xlabel('$\mu_L^{(C)}-\mu_R^{(C)}$ (meV)')
#	plt.ylabel('$A$ (meV)')
#	plt.title('(d) Down')
	p8=fig.add_subplot(224,projection='3d')
	p8.plot_surface(X*1e3, Y*1e3, C_D, rstride=1, cstride=1,linewidth=0.5,  alpha=0.3,cmap=cm.coolwarm)
	cset = p8.contour(X*1e3, Y*1e3, C_D,[0], colors='k', zdir='z', offset=0)
	cset.collections[0].set_linewidth(2)  
	p8.set_xlabel('\n'+r'$e \bar{V}_b$ (meV)',linespacing=1)
	p8.set_ylabel('\n $A$ (meV)',linespacing=1)
	p8.zaxis.set_rotate_label(False)
	p8.set_zlabel(r'$I_\downarrow$ (pA)',rotation=90)
	p8.set_zlim(-80,120)
	p8.xaxis.set_major_locator(MultipleLocator(10))
	p8.yaxis.set_major_locator(MultipleLocator(50))
	p8.zaxis.set_major_locator(MultipleLocator(60))
	p8.view_init(elev=40, azim = 50)
	plt.title('\n (d) Down',linespacing=1)        
	
	plt.tight_layout() 

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

#centerchange1()	

#wL=1e0
#wR=1e0
#muL=0
#muR=0
#print(dynamic(muL,muR,0e-3)[0]*1.60217662*1e-4,dynamic(muL,muR,0e-3)[1]*1.60217662*1e-4)
#print(geometric(muL,muR,0e-3))


geometric(0,0,0)

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
