<<<<<<< HEAD
# degenerate perturbation theory (small chi)

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

global T,Eu,Ed,h,hz,wL,wR, a, b ,dchi,Tp,omega
#T=5e-3
#Ed=3.498
#Eu=3.441


Ed=3.270
Eu=3.149
shift=0.0
scale=1e0
T=10e-3
Ed=scale*(Ed-Eu+shift)
Eu=scale*(0.0+shift)
h=0e-3
hz=0

dchi=0e-4+0j
print(dchi)

omega=1e9
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






def geometric(muLC,muRC,amp):
	def fLinU(muL,muR):
	    f=1.0/(1.0+np.exp((Eu+hz/2.0-muL)/T))
	    return f
	
	def fLoutU(muL,muR):
	    f=1.0/(1.0+np.exp(-(Eu+hz/2.0-muL)/T))
	    return f
	
	def fRinU(muL,muR):
	    f=1.0/(1.0+np.exp((Eu+hz/2.0-muR)/T))
	    return f
	
	def fRoutU(muL,muR):
	    f=1.0/(1.0+np.exp(-(Eu+hz/2.0-muR)/T))
	    return f
	
	def fLinD(muL,muR):
	    f=0.0
	    return f
	
	def fLoutD(muL,muR):
	    f=1.0
	    return f
	
	def fRinD(muL,muR):
	    f=0.0
	    return f
	
	def fRoutD(muL,muR):
	    f=1.0
	    return f
	
	    
	def H0(muL,muR):
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
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)



	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0p
	    H[0,2]=W0m
	    H[1,0]=Wp0
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0
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

	def H1D(muL,muR):
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
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)



	    H[0,0]= 0
	    H[0,1]=0
	    H[0,2]=W0mR
	    H[1,0]=0
	    H[1,1]=0
	    H[1,2]=0
	    H[2,0]=-Wm0R
	    H[2,1]=0
	    H[2,2]=0
	    H=-H
	    
	#    print('0pL',W0pL)
	#    print('p0L',Wp0L)
	#    print('0pR',W0pR)
	#    print('p0R',Wp0R)
	#    print(H[0,1])
	#    print(H[1,0])
	    return H

	def H1U(muL,muR):
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
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)



	    H[0,0]= 0
	    H[0,1]=W0pR
	    H[0,2]=0
	    H[1,0]=-Wp0R
	    H[1,1]=0
	    H[1,2]=0
	    H[2,0]=0
	    H[2,1]=0
	    H[2,2]=0
	    H=-H
	    
	#    print('0pL',W0pL)
	#    print('p0L',Wp0L)
	#    print('0pR',W0pR)
	#    print('p0R',Wp0R)
	#    print(H[0,1])
	#    print(H[1,0])
	    return H
	    
	def fLinUdiff(muL,muR):
	    a=np.exp((Eu+hz/2.0-muL)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=1.0/(1.0/a+2+a)/T
	    return f
	
	def fLoutUdiff(muL,muR):
	    a=np.exp((Eu+hz/2.0-muL)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=-1.0/(1.0/a+2+a)/T
	    return f
	
	def fRinUdiff(muL,muR):
	    a=np.exp((Eu+hz/2.0-muR)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=1.0/(1.0/a+2+a)/T
	    return f
	
	def fRoutUdiff(muL,muR):
	    a=np.exp((Eu+hz/2.0-muR)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=-1.0/(1.0/a+2+a)/T
	    return f
	
	def fLinDdiff(muL,muR):
	    a=np.exp((Ed-hz/2.0-muL)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=0.0
	    return f
	
	def fLoutDdiff(muL,muR):
	    a=np.exp((Ed-hz/2.0-muL)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=0.0
	    return f
	
	def fRinDdiff(muL,muR):
	    a=np.exp((Ed-hz/2.0-muR)/T)
#	    f=a/(1.0+a)**2.0/T
	    f=0.0
	    return f
	
	def fRoutDdiff(muL,muR):
	    a=np.exp((Ed-hz/2.0-muR)/T)
#	    f=-a/(1.0+a)**2.0/T
	    f=0.0
	    return f
	
	    
	def H0diffL(muL,muR):
	    H=np.zeros((3,3),dtype=complex)
	
	    a=eigofH(hz,h)[0]
	    b=eigofH(hz,h)[1]
	
	    Wp0L=wL*(fLinUdiff(muL,muR)*a**2+fLinDdiff(muL,muR)*b**2)
	    Wp0R=0
	    Wm0L=wL*(fLinUdiff(muL,muR)*b**2+fLinDdiff(muL,muR)*a**2)
	    Wm0R=0
	    
	    W0pL=wL*(fLoutUdiff(muL,muR)*a**2+fLoutDdiff(muL,muR)*b**2)
	    W0pR=0
	    W0mL=wL*(fLoutUdiff(muL,muR)*b**2+fLoutDdiff(muL,muR)*a**2)
	    W0mR=0
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)


	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0pL+W0pR
	    H[0,2]=W0mL+W0mR
	    H[1,0]=Wp0L+Wp0R
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0L+Wm0R
	    H[2,1]=0
	    H[2,2]= -W0m
	    H=-H
	    return H



	
	def H0diffR(muL,muR):
	    H=np.zeros((3,3),dtype=complex)
	
	    a=eigofH(hz,h)[0]
	    b=eigofH(hz,h)[1]
	
	    Wp0L=0
	    Wp0R=wR*(fRinUdiff(muL,muR)*a**2+fRinDdiff(muL,muR)*b**2)
	    Wm0L=0
	    Wm0R=wR*(fRinUdiff(muL,muR)*b**2+fRinDdiff(muL,muR)*a**2)
	    
	    W0pL=0
	    W0pR=wR*(fRoutUdiff(muL,muR)*a**2+fRoutDdiff(muL,muR)*b**2)
	    W0mL=0
	    W0mR=wR*(fRoutUdiff(muL,muR)*b**2+fRoutDdiff(muL,muR)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)

	    H[0,0]= -Wp0-Wm0
	    H[0,1]=W0p
	    H[0,2]=W0m
	    H[1,0]=Wp0
	    H[1,1]= -W0p
	    H[1,2]=0
	    H[2,0]=Wm0
	    H[2,1]=0
	    H[2,2]= -W0m
	    H=-H
	    return H

	def H1DdiffR(muL,muR):
	    H=np.zeros((3,3),dtype=complex)
	
	    a=eigofH(hz,h)[0]
	    b=eigofH(hz,h)[1]
	
	    Wp0L=0
	    Wp0R=wR*(fRinUdiff(muL,muR)*a**2+fRinDdiff(muL,muR)*b**2)
	    Wm0L=0
	    Wm0R=wR*(fRinUdiff(muL,muR)*b**2+fRinDdiff(muL,muR)*a**2)
	    
	    W0pL=0
	    W0pR=wR*(fRoutUdiff(muL,muR)*a**2+fRoutDdiff(muL,muR)*b**2)
	    W0mL=0
	    W0mR=wR*(fRoutUdiff(muL,muR)*b**2+fRoutDdiff(muL,muR)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)

	    H[0,0]=0
	    H[0,1]=0
	    H[0,2]=W0mR
	    H[1,0]=0
	    H[1,1]=0
	    H[1,2]=0
	    H[2,0]=-Wm0R
	    H[2,1]=0
	    H[2,2]=0
	    H=-H
	    return H

	def H1UdiffR(muL,muR):
	    H=np.zeros((3,3),dtype=complex)
	
	    a=eigofH(hz,h)[0]
	    b=eigofH(hz,h)[1]
	
	    Wp0L=0
	    Wp0R=wR*(fRinUdiff(muL,muR)*a**2+fRinDdiff(muL,muR)*b**2)
	    Wm0L=0
	    Wm0R=wR*(fRinUdiff(muL,muR)*b**2+fRinDdiff(muL,muR)*a**2)
	    
	    W0pL=0
	    W0pR=wR*(fRoutUdiff(muL,muR)*a**2+fRoutDdiff(muL,muR)*b**2)
	    W0mL=0
	    W0mR=wR*(fRoutUdiff(muL,muR)*b**2+fRoutDdiff(muL,muR)*a**2)
	   
	    W0p=W0pL+W0pR
	    W0m=W0mL+W0mR
	    Wp0=Wp0L+Wp0R
	    Wm0=Wm0L+Wm0R
#	    e_chip=np.exp(1j*chip)
#	    e_chip2=np.exp(-1j*chip)
#	    e_chim=np.exp(1j*chim)
#	    e_chim2=np.exp(-1j*chim)

	    H[0,0]= 0
	    H[0,1]=W0pR
	    H[0,2]=0
	    H[1,0]=-Wp0R
	    H[1,1]=0
	    H[1,2]=0
	    H[2,0]=0
	    H[2,1]=0
	    H[2,2]=0
	    H=-H
	    return H
	
	def BerryCurvature_D(muL,muR):
	#    H=Hamilt(t,chi)
	#    Hdagger=np.conj(H.transpose())
	#    print(H)
	#    print(Hdagger)
	
	    H=H0(muL,muR)
#	    print("H=",H)
#	    H=np.array([[0,0,0],[0,1,0],[0,0,2]])
	    A=eigscipy(H,left=True)

	    imin=A[0].argsort()[0]
	    imed=A[0].argsort()[1]
	    imax=A[0].argsort()[2]
	    eig0=A[0][imin]
	    eig1=A[0][imed]
	    eig2=A[0][imax]
	    
#	    print(eig0,eig1,eig2)
	    ev0R=A[2].transpose()[imin]
	    ev0L=A[1].transpose()[imin]
	    ev1R=A[2].transpose()[imed]
	    ev1L=A[1].transpose()[imed]
	    ev2R=A[2].transpose()[imax]
	    ev2L=A[1].transpose()[imax]
	    ev0L=np.conj(ev0L)
	    ev1L=np.conj(ev1L)
	    ev2L=np.conj(ev2L)
	
	
		  
	    ev0L=ev0L/np.matmul(ev0L,ev0R)
	    ev1L=ev1L/np.matmul(ev1L,ev1R)
	    ev2L=ev2L/np.matmul(ev2L,ev2R)

	    if (abs(eig2-eig1)>1e-10):
	    	#print(eig0,eig1,eig2)
	    	ev0R=A[2].transpose()[imin]
	    	ev0L=A[1].transpose()[imin]
	    	ev1R=A[2].transpose()[imed]
	    	ev1L=A[1].transpose()[imed]
	    	ev2R=A[2].transpose()[imax]
	    	ev2L=A[1].transpose()[imax]
	    	ev0L=np.conj(ev0L)
	    	ev1L=np.conj(ev1L)
	    	ev2L=np.conj(ev2L)
	
	
		    
	    	ev0L=ev0L/np.matmul(ev0L,ev0R)
	    	ev1L=ev1L/np.matmul(ev1L,ev1R)
	    	ev2L=ev2L/np.matmul(ev2L,ev2R)
	    	
	    	
	    	
	    	Ha1=H1D(muL,muR)
	    	#	    Ha1=np.array([[0,0,1],[0,0,0],[-1,0,0]])
	    	
	    	# i=0
	    	temp0R=np.matmul(Ha1,ev0R)
	    	# j=1
	    	temp10=np.matmul(ev1L,temp0R)
	    	# j=2
	    	temp20=np.matmul(ev2L,temp0R)
	    	ev0Rp=temp10/(eig0-eig1)*ev1R+temp20/(eig0-eig2)*ev2R
	    	
	    	# i=1
	    	temp1R=np.matmul(Ha1,ev1R)
	    	# j=0
	    	temp01=np.matmul(ev0L,temp1R)
	    	# j=2
	    	temp21=np.matmul(ev2L,temp1R)
	    	ev1Rp=temp01/(eig1-eig0)*ev0R+temp21/(eig1-eig2)*ev2R
	    	
	    	# i=2
	    	temp2R=np.matmul(Ha1,ev2R)
	    	# j=0
	    	temp02=np.matmul(ev0L,temp2R)
	    	# j=1
	    	temp12=np.matmul(ev1L,temp2R)
	    	ev2Rp=temp02/(eig2-eig0)*ev0R+temp12/(eig2-eig1)*ev1R	    
	    	
	    	
	    	# i=0
	    	#	    temp0L=np.matmul(ev0L,Ha1)
	    	# j=1
	    	#	    temp01=np.matmul(temp0L,ev1R)
	    	# j=2
	    	#	    temp02=np.matmul(temp0L,ev2R)
	    	ev0Lp=temp01/(eig0-eig1)*ev1L+temp02/(eig0-eig2)*ev2L
	    	
	    	
	    	#	    print("ev0L=",ev0L)
	    	
	    	
	    	
	    	# i=1
	    	#	    temp1L=np.matmul(ev1L,Ha1)
	    	# j=0
	    	#	    temp10=np.matmul(temp1L,ev0R)
	    	# j=2
	    	#	    temp12=np.matmul(temp1L,ev2R)
	    	ev1Lp=temp10/(eig1-eig0)*ev0L+temp12/(eig1-eig2)*ev2L
	    	
	    		# i=2
	    	#	    temp2L=np.matmul(ev2L,Ha1)
	    	# j=0
	    	#	    temp20=np.matmul(temp2L,ev0R)
	    	# j=1
	    	#	    temp21=np.matmul(temp2L,ev1R)
	    	ev2Lp=temp20/(eig2-eig0)*ev0L+temp21/(eig2-eig1)*ev1L
	    	
	    	ev0Lchi=ev0L+ev0Lp
	    	ev1Rchi=ev1R+ev1Rp

		  #########################3
		    
	    if (abs(eig2-eig1)<1e-10):
#	    	print(muL,muR)
	    	print(eig0,eig1,eig2)
	    	
	    	
	    	
	    	Ha1=H1D(muL,muR)
	    	H11=np.matmul(np.matmul(ev1L,Ha1),ev1R)
	    	H12=np.matmul(np.matmul(ev1L,Ha1),ev2R)
	    	H21=np.matmul(np.matmul(ev2L,Ha1),ev1R)
	    	H22=np.matmul(np.matmul(ev2L,Ha1),ev2R)
	    	H1=np.array([[H11,H12],[H21,H22]])
#	    	print("H1=",H1)
#	    	print("det=",np.linalg.det(H1))
	    	# diagnolize H1
	    	B=eigscipy(H1)
#	    	print("eig=",B[0])
	    	eig1p=B[0][0]
	    	eig2p=B[0][1]
	    	transf=B[1]
	    	
	    	inverse=LA.inv(transf)
#	    	print("inverse",inverse)
#	    	print("D=",np.matmul(np.matmul(inverse,H1),transf))


	    	#print("ident=",np.matmul(transf,np.conj(transf.transpose())))
#	    	print("Delta=",(H11-H22)*(H11-H22)+4*H12*H21)
#	    	print("eig analytical=",0.5*(H11+H22+np.sqrt((H11-H22)*(H11-H22)+4*H12*H21)))
#	    	print("det=",B[0][0]*B[0][1])

	    	evaL=(inverse[0,0])*ev1L+(inverse[0,1])*ev2L
	    	evbL=(inverse[1,0])*ev1L+(inverse[1,1])*ev2L
	    	ev1L=evaL
	    	ev2L=evbL
	    	evaR=B[1][0,0]*ev1R+B[1][1,0]*ev2R
	    	evbR=B[1][0,1]*ev1R+B[1][1,1]*ev2R
	    	ev1R=evaR
	    	ev2R=evbR
	    	H11=np.matmul(np.matmul(ev1L,Ha1),ev1R)
	    	H12=np.matmul(np.matmul(ev1L,Ha1),ev2R)
	    	H21=np.matmul(np.matmul(ev2L,Ha1),ev1R)
	    	H22=np.matmul(np.matmul(ev2L,Ha1),ev2R)
	    	H1=np.array([[H11,H12],[H21,H22]])
	    	#print("D=",H1)
	    	
	    	# i=0 (same with non-degenerate)
	    	temp0R=np.matmul(Ha1,ev0R)
	    	# j=1
	    	temp10=np.matmul(ev1L,temp0R)
	    	# j=2
	    	temp20=np.matmul(ev2L,temp0R)
	    	
	    	# i=1
	    	temp1R=np.matmul(Ha1,ev1R)
	    	# j=0
	    	temp01=np.matmul(ev0L,temp1R)
	    	# j=2
	    	temp21=np.matmul(ev2L,temp1R)
	    	
	    	# i=2
	    	temp2R=np.matmul(Ha1,ev2R)
	    	# j=0
	    	temp02=np.matmul(ev0L,temp2R)
	    	# j=1
	    	temp12=np.matmul(ev1L,temp2R)

	    	ev0Rp=temp10/(eig0-eig1)*ev1R+temp20/(eig0-eig2)*ev2R
	    	ev1Rp=temp01/(eig1-eig0)*ev0R+temp20*temp01/(eig1p-eig2p)/(eig1-eig0)*ev2R
	    	ev2Rp=temp02/(eig2-eig0)*ev0R+temp10*temp02/(eig2p-eig1p)/(eig2-eig0)*ev1R    
	    	
	    	

	    	ev0Lp=temp01/(eig0-eig1)*ev1L+temp02/(eig0-eig2)*ev2L
	    	ev1Lp=temp10/(eig1-eig0)*ev0L+temp10*temp02/(eig1p-eig2p)/(eig1-eig0)*ev2L
	    	ev2Lp=temp20/(eig2-eig0)*ev0L+temp20*temp01/(eig2p-eig1p)/(eig2-eig0)*ev1L
	    	
	    	ev0Lchi=ev0L+ev0Lp
	    	ev1Rchi=ev1R+ev1Rp	


			###################################

	    H0L=H0diffL(muL,muR)
	    H1DL=np.zeros((3,3),dtype=complex)
	    H0R=H0diffR(muL,muR)
	    H1DR=H1DdiffR(muL,muR)

#	    print("ev0L=",ev0L)
#	    print("a11=",np.matmul(np.matmul(ev0L,H0L),ev1R))
	    

# i=1
	    b11=np.matmul(np.matmul(ev0Lp,H0L),ev1R)+\
	        np.matmul(np.matmul(ev0L,H1DL),ev1R)+\
	        np.matmul(np.matmul(ev0L,H0L),ev1Rp)

#	    print("b11=",b11*1j*dchi)
	    
	    a12=np.matmul(np.matmul(ev1L,H0R),ev0R)
	    
#	    print("a12=",a12)
	    b12=np.matmul(np.matmul(ev1Lp,H0R),ev0R)+\
	        np.matmul(np.matmul(ev1L,H1DR),ev0R)+\
	        np.matmul(np.matmul(ev1L,H0R),ev0Rp)
#	    print("b12=",b12)
#	    print(a12+1j*dchi*b12)
	    
	    
	    b112=b11*a12
	    
	    b13=np.matmul(np.matmul(ev0Lp,H0R),ev1R)+\
	        np.matmul(np.matmul(ev0L,H1DR),ev1R)+\
	        np.matmul(np.matmul(ev0L,H0R),ev1Rp)
	    a14=np.matmul(np.matmul(ev1L,H0L),ev0R)
	    b134=b13*a14
	    
	    f1=eig0-eig1
	    BC_D1=(b112-b134)/(f1*f1)
	    
# i=2
	    b11=np.matmul(np.matmul(ev0Lp,H0L),ev2R)+\
	        np.matmul(np.matmul(ev0L,H1DL),ev2R)+\
	        np.matmul(np.matmul(ev0L,H0L),ev2Rp)
	    a12=np.matmul(np.matmul(ev2L,H0R),ev0R)
	    b112=b11*a12
	    
	    b13=np.matmul(np.matmul(ev0Lp,H0R),ev2R)+\
	        np.matmul(np.matmul(ev0L,H1DR),ev2R)+\
	        np.matmul(np.matmul(ev0L,H0R),ev2Rp)
	    a14=np.matmul(np.matmul(ev2L,H0L),ev0R)
	    b134=b13*a14
	    
	    f1=eig0-eig2
	    BC_D2=(b112-b134)/(f1*f1)
	    
	    BC_D=BC_D1+BC_D2


	    return BC_D
	
#	BerryCurvature_D(-0.05,-0.03895)
		
	def BerryCurvature_U(muL,muR):
	#    H=Hamilt(t,chi)
	#    Hdagger=np.conj(H.transpose())
	#    print(H)
	#    print(Hdagger)
	
	    H=H0(muL,muR)
#	    print("H=",H)
#	    H=np.array([[0,0,0],[0,1,0],[0,0,2]])
	    A=eigscipy(H,left=True)

	    imin=A[0].argsort()[0]
	    imed=A[0].argsort()[1]
	    imax=A[0].argsort()[2]
	    eig0=A[0][imin]
	    eig1=A[0][imed]
	    eig2=A[0][imax]
	    
#	    print(eig0,eig1,eig2)
	    ev0R=A[2].transpose()[imin]
	    ev0L=A[1].transpose()[imin]
	    ev1R=A[2].transpose()[imed]
	    ev1L=A[1].transpose()[imed]
	    ev2R=A[2].transpose()[imax]
	    ev2L=A[1].transpose()[imax]
	    ev0L=np.conj(ev0L)
	    ev1L=np.conj(ev1L)
	    ev2L=np.conj(ev2L)
	
	
		  
	    ev0L=ev0L/np.matmul(ev0L,ev0R)
	    ev1L=ev1L/np.matmul(ev1L,ev1R)
	    ev2L=ev2L/np.matmul(ev2L,ev2R)

	    if (abs(eig2-eig1)>1e-10):
	    	#print(eig0,eig1,eig2)
	    	ev0R=A[2].transpose()[imin]
	    	ev0L=A[1].transpose()[imin]
	    	ev1R=A[2].transpose()[imed]
	    	ev1L=A[1].transpose()[imed]
	    	ev2R=A[2].transpose()[imax]
	    	ev2L=A[1].transpose()[imax]
	    	ev0L=np.conj(ev0L)
	    	ev1L=np.conj(ev1L)
	    	ev2L=np.conj(ev2L)
	
	
		    
	    	ev0L=ev0L/np.matmul(ev0L,ev0R)
	    	ev1L=ev1L/np.matmul(ev1L,ev1R)
	    	ev2L=ev2L/np.matmul(ev2L,ev2R)
	    	
	    	
	    	
	    	Ha1=H1U(muL,muR)
	    	#	    Ha1=np.array([[0,0,1],[0,0,0],[-1,0,0]])
	    	
	    	# i=0
	    	temp0R=np.matmul(Ha1,ev0R)
	    	# j=1
	    	temp10=np.matmul(ev1L,temp0R)
	    	# j=2
	    	temp20=np.matmul(ev2L,temp0R)
	    	ev0Rp=temp10/(eig0-eig1)*ev1R+temp20/(eig0-eig2)*ev2R
	    	
	    	# i=1
	    	temp1R=np.matmul(Ha1,ev1R)
	    	# j=0
	    	temp01=np.matmul(ev0L,temp1R)
	    	# j=2
	    	temp21=np.matmul(ev2L,temp1R)
	    	ev1Rp=temp01/(eig1-eig0)*ev0R+temp21/(eig1-eig2)*ev2R
	    	
	    	# i=2
	    	temp2R=np.matmul(Ha1,ev2R)
	    	# j=0
	    	temp02=np.matmul(ev0L,temp2R)
	    	# j=1
	    	temp12=np.matmul(ev1L,temp2R)
	    	ev2Rp=temp02/(eig2-eig0)*ev0R+temp12/(eig2-eig1)*ev1R	    
	    	
	    	
	    	# i=0
	    	#	    temp0L=np.matmul(ev0L,Ha1)
	    	# j=1
	    	#	    temp01=np.matmul(temp0L,ev1R)
	    	# j=2
	    	#	    temp02=np.matmul(temp0L,ev2R)
	    	ev0Lp=temp01/(eig0-eig1)*ev1L+temp02/(eig0-eig2)*ev2L
	    	
	    	
	    	#	    print("ev0L=",ev0L)
	    	
	    	
	    	
	    	# i=1
	    	#	    temp1L=np.matmul(ev1L,Ha1)
	    	# j=0
	    	#	    temp10=np.matmul(temp1L,ev0R)
	    	# j=2
	    	#	    temp12=np.matmul(temp1L,ev2R)
	    	ev1Lp=temp10/(eig1-eig0)*ev0L+temp12/(eig1-eig2)*ev2L
	    	
	    		# i=2
	    	#	    temp2L=np.matmul(ev2L,Ha1)
	    	# j=0
	    	#	    temp20=np.matmul(temp2L,ev0R)
	    	# j=1
	    	#	    temp21=np.matmul(temp2L,ev1R)
	    	ev2Lp=temp20/(eig2-eig0)*ev0L+temp21/(eig2-eig1)*ev1L
	    	
	    	ev0Lchi=ev0L+ev0Lp
	    	ev1Rchi=ev1R+ev1Rp

		  #########################3
		    
	    if (abs(eig2-eig1)<1e-10):
#	    	print(muL,muR)
	    	print(eig0,eig1,eig2)
	    	
	    	
	    	
	    	Ha1=H1U(muL,muR)
	    	H11=np.matmul(np.matmul(ev1L,Ha1),ev1R)
	    	H12=np.matmul(np.matmul(ev1L,Ha1),ev2R)
	    	H21=np.matmul(np.matmul(ev2L,Ha1),ev1R)
	    	H22=np.matmul(np.matmul(ev2L,Ha1),ev2R)
	    	H1=np.array([[H11,H12],[H21,H22]])
#	    	print("H1=",H1)
#	    	print("det=",np.linalg.det(H1))
	    	# diagnolize H1
	    	B=eigscipy(H1)
#	    	print("eig=",B[0])
	    	eig1p=B[0][0]
	    	eig2p=B[0][1]
	    	transf=B[1]
	    	
	    	inverse=LA.inv(transf)
#	    	print("inverse",inverse)
#	    	print("D=",np.matmul(np.matmul(inverse,H1),transf))


	    	#print("ident=",np.matmul(transf,np.conj(transf.transpose())))
#	    	print("Delta=",(H11-H22)*(H11-H22)+4*H12*H21)
#	    	print("eig analytical=",0.5*(H11+H22+np.sqrt((H11-H22)*(H11-H22)+4*H12*H21)))
#	    	print("det=",B[0][0]*B[0][1])

	    	evaL=(inverse[0,0])*ev1L+(inverse[0,1])*ev2L
	    	evbL=(inverse[1,0])*ev1L+(inverse[1,1])*ev2L
	    	ev1L=evaL
	    	ev2L=evbL
	    	evaR=B[1][0,0]*ev1R+B[1][1,0]*ev2R
	    	evbR=B[1][0,1]*ev1R+B[1][1,1]*ev2R
	    	ev1R=evaR
	    	ev2R=evbR
	    	H11=np.matmul(np.matmul(ev1L,Ha1),ev1R)
	    	H12=np.matmul(np.matmul(ev1L,Ha1),ev2R)
	    	H21=np.matmul(np.matmul(ev2L,Ha1),ev1R)
	    	H22=np.matmul(np.matmul(ev2L,Ha1),ev2R)
	    	H1=np.array([[H11,H12],[H21,H22]])
	    	#print("D=",H1)
	    	
	    	# i=0 (same with non-degenerate)
	    	temp0R=np.matmul(Ha1,ev0R)
	    	# j=1
	    	temp10=np.matmul(ev1L,temp0R)
	    	# j=2
	    	temp20=np.matmul(ev2L,temp0R)
	    	
	    	# i=1
	    	temp1R=np.matmul(Ha1,ev1R)
	    	# j=0
	    	temp01=np.matmul(ev0L,temp1R)
	    	# j=2
	    	temp21=np.matmul(ev2L,temp1R)
	    	
	    	# i=2
	    	temp2R=np.matmul(Ha1,ev2R)
	    	# j=0
	    	temp02=np.matmul(ev0L,temp2R)
	    	# j=1
	    	temp12=np.matmul(ev1L,temp2R)

	    	ev0Rp=temp10/(eig0-eig1)*ev1R+temp20/(eig0-eig2)*ev2R
	    	ev1Rp=temp01/(eig1-eig0)*ev0R+temp20*temp01/(eig1p-eig2p)/(eig1-eig0)*ev2R
	    	ev2Rp=temp02/(eig2-eig0)*ev0R+temp10*temp02/(eig2p-eig1p)/(eig2-eig0)*ev1R    
	    	
	    	

	    	ev0Lp=temp01/(eig0-eig1)*ev1L+temp02/(eig0-eig2)*ev2L
	    	ev1Lp=temp10/(eig1-eig0)*ev0L+temp10*temp02/(eig1p-eig2p)/(eig1-eig0)*ev2L
	    	ev2Lp=temp20/(eig2-eig0)*ev0L+temp20*temp01/(eig2p-eig1p)/(eig2-eig0)*ev1L
	    	
	    	ev0Lchi=ev0L+ev0Lp
	    	ev1Rchi=ev1R+ev1Rp	


			###################################

	    H0L=H0diffL(muL,muR)
	    H1UL=np.zeros((3,3),dtype=complex)
	    H0R=H0diffR(muL,muR)
	    H1UR=H1UdiffR(muL,muR)

#	    print("ev0L=",ev0L)
#	    print("a11=",np.matmul(np.matmul(ev0L,H0L),ev1R))
	    

# i=1
	    b11=np.matmul(np.matmul(ev0Lp,H0L),ev1R)+\
	        np.matmul(np.matmul(ev0L,H1UL),ev1R)+\
	        np.matmul(np.matmul(ev0L,H0L),ev1Rp)

#	    print("b11=",b11*1j*dchi)
	    
	    a12=np.matmul(np.matmul(ev1L,H0R),ev0R)
	    
#	    print("a12=",a12)
	    b12=np.matmul(np.matmul(ev1Lp,H0R),ev0R)+\
	        np.matmul(np.matmul(ev1L,H1UR),ev0R)+\
	        np.matmul(np.matmul(ev1L,H0R),ev0Rp)
#	    print("b12=",b12)
#	    print(a12+1j*dchi*b12)
	    
	    
	    b112=b11*a12
	    
	    b13=np.matmul(np.matmul(ev0Lp,H0R),ev1R)+\
	        np.matmul(np.matmul(ev0L,H1UR),ev1R)+\
	        np.matmul(np.matmul(ev0L,H0R),ev1Rp)
	    a14=np.matmul(np.matmul(ev1L,H0L),ev0R)
	    b134=b13*a14
	    
	    f1=eig0-eig1
	    BC_U1=(b112-b134)/(f1*f1)
	    
# i=2
	    b11=np.matmul(np.matmul(ev0Lp,H0L),ev2R)+\
	        np.matmul(np.matmul(ev0L,H1UL),ev2R)+\
	        np.matmul(np.matmul(ev0L,H0L),ev2Rp)
	    a12=np.matmul(np.matmul(ev2L,H0R),ev0R)
	    b112=b11*a12
	    
	    b13=np.matmul(np.matmul(ev0Lp,H0R),ev2R)+\
	        np.matmul(np.matmul(ev0L,H1UR),ev2R)+\
	        np.matmul(np.matmul(ev0L,H0R),ev2Rp)
	    a14=np.matmul(np.matmul(ev2L,H0L),ev0R)
	    b134=b13*a14
	    
	    f1=eig0-eig2
	    BC_U2=(b112-b134)/(f1*f1)
	    
	    BC_U=BC_U1+BC_U2


	    return BC_U		

	def fmt(x, pos):
#  		a, b = '{:.2e}'.format(x).split('e')
#  		b = int(b)
#  		return r'${} \times 10^{{{}}}$'.format(a, b)
			return '{0:.1E}'.format(x).replace('+0', '').replace('-0', '-')
    
	def BC_plot1():
	    Ngrid=21 #41
	    xlist = np.linspace(Eu-50e-3, Ed+50e-3, Ngrid,endpoint=True)
	    ylist = np.linspace(Eu-50e-3, Ed+50e-3, Ngrid,endpoint=True)
	    X, Y = np.meshgrid(xlist, ylist)
	    BC_C=np.zeros((Ngrid,Ngrid),dtype=complex)
	    BC_S=np.zeros((Ngrid,Ngrid),dtype=complex)
	    BC_U=np.zeros((Ngrid,Ngrid),dtype=complex)
	    BC_D=np.zeros((Ngrid,Ngrid),dtype=complex)
	    for i in range(Ngrid):
	        for j in range(Ngrid):
	            BC_D[j,i]=BerryCurvature_D(X[0][i],Y[j][0])
#	            BC_C[i,j]=(temp[0]+temp[1])
#	            BC_S[i,j]=(temp[0]-temp[1])
	            BC_U[i,j]=BerryCurvature_U(X[0][i],Y[j][0])
#	            BC_D[i,j]=(temp[1])

##### attention, the order of index


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
	    plt.xticks([0,50,100])
	    plt.yticks([0,50,100])

#	    plt.gca().set_aspect('equal', adjustable='box')

	    p4=fig.add_subplot(122)
	    cp = p4.contourf(X*1e3, Y*1e3, BC_D,16)
	    plt.colorbar(cp,format=ticker.FuncFormatter(fmt))
	    plt.title('(b) Down')
	    plt.xlabel('$\mu_L$ (meV)',fontsize=19)
	    plt.ylabel('$\mu_R$ (meV)',fontsize=19)
	    plt.xticks([0,50,100])
	    plt.yticks([0,50,100])
#	    plt.gca().set_aspect('equal', adjustable='box')

	    plt.tight_layout()
	    plt.show()
	    return

	BC_plot1()

#	return current()
	return



wL=1e0
wR=1e0
muL=0
muR=0
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



