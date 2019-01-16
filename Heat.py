#%matplotlib inline

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from numpy import linalg as LA
from scipy import linalg as LAs
from scipy.linalg import eig as eigscipy
import time
start_time = time.time()


global hw,k,dchi
hw=25e-3
k=8.617e-5







dchi=1e-6+0j





def NL(TL):
    f=1.0/(np.exp(hw/k/TL)-1.0)
    return f

def NR(TR):
    f=1.0/(np.exp(hw/k/TR)-1.0)
    return f

    
def Hamilt(TL,TR,chi):
    H=np.zeros((2,2),dtype=complex)
    k01L=NL(TL)
    k01R=NR(TR)
    k10L=NL(TL)+1.0
    k10R=NR(TR)+1.0
    
    e_chi=np.exp(1j*chi)
    e_chi2=np.exp(-1j*chi)
#    e_chi=1.0+1j*chi
#    e_chi2=1.0-1j*chi
    H[0,0]= k01L+k01R
    H[0,1]=-k10L-k10R*e_chi
    H[1,0]=-k01L-k01R*e_chi2
    H[1,1]= k10L+k10R
    
#    print('0pL',W0pL)
#    print('p0L',Wp0L)
#    print('0pR',W0pR)
#    print('p0R',Wp0R)
#    print(H[0,1])
#    print(H[1,0])
    return H
    
def NLdiff(TL,TR):
    a=np.exp(hw/k/TL)
    f=a*hw/k/TL**2.0/(a-1.0)**2.0
    return f

def NRdiff(TL,TR):
    a=np.exp(hw/k/TR)
    f=a*hw/k/TR**2.0/(a-1.0)**2.0
    return f

    
def HdiffL(TL,TR,chi):
    H=np.zeros((2,2),dtype=complex)
    k01L=NLdiff(TL,TR)
    k01R=0
    k10L=NLdiff(TL,TR)
    k10R=0
    
    e_chi=np.exp(1j*chi)
    e_chi2=np.exp(-1j*chi)
#    e_chi=1.0+1j*chi
#    e_chi2=1.0-1j*chi
    H[0,0]= k01L+k01R
    H[0,1]=-k10L-k10R*e_chi
    H[1,0]=-k01L-k01R*e_chi2
    H[1,1]= k10L+k10R
    return H

def HdiffR(TL,TR,chi):
    H=np.zeros((2,2),dtype=complex)
    k01L=0
    k01R=NRdiff(TL,TR)
    k10L=0
    k10R=NRdiff(TL,TR)
    
    e_chi=np.exp(1j*chi)
    e_chi2=np.exp(-1j*chi)
#    e_chi=1.0+1j*chi
#    e_chi2=1.0-1j*chi
    H[0,0]= k01L+k01R
    H[0,1]=-k10L-k10R*e_chi
    H[1,0]=-k01L-k01R*e_chi2
    H[1,1]= k10L+k10R
    return H

def F(TL,TR,chi):
#    H=Hamilt(t,chi)
#    Hdagger=np.conj(H.transpose())
#    print(H)
#    print(Hdagger)

    H=Hamilt(TL,TR,chi)
    A=eigscipy(H,left=True)
    imin=A[0].argsort()[0]
    imed=A[0].argsort()[1]
#    imax=A[0].argsort()[2]
    eig0=A[0][imin]
    eig1=A[0][imed]
#    eig2=A[0][imax]
    phi0R=A[2].transpose()[imin]
    phi0L=A[1].transpose()[imin]
    phi1R=A[2].transpose()[imed]
    phi1L=A[1].transpose()[imed]
#    phi2R=A[2].transpose()[imax]
#    phi2L=A[1].transpose()[imax]
    
    phi0L=np.conj(phi0L)
    phi1L=np.conj(phi1L)
#    phi2L=np.conj(phi2L)
    
    phi0L=phi0L/np.matmul(phi0L,phi0R)
    phi1L=phi1L/np.matmul(phi1L,phi1R)
#    phi2L=phi2L/np.matmul(phi2L,phi2R)

#    phi1L=LA.eig(H1L_T)[1].transpose()[0]
#    phi1L=np.conj(phi1L)

#    print('eig of H1L',LA.eig(H1L))
#    print('eig of H1L_T',LA.eig(H1L_T))
#    print('left',phi1L)
#    print('left from scipy',eigscipy(H1L,left=True))
#    print('right',LA.eig(H1L)[1].transpose()[0])
    
    a=np.matmul(phi0L,HdiffL(TL,TR,chi))
    b1L=np.matmul(a,phi1R)
    a=np.matmul(phi1L,HdiffR(TL,TR,chi))
    c1L=np.matmul(a,phi0R)
    
    a=np.matmul(phi0L,HdiffR(TL,TR,chi))
    b1R=np.matmul(a,phi1R)
    a=np.matmul(phi1L,HdiffL(TL,TR,chi))
    c1R=np.matmul(a,phi0R)
    
    element1=(b1L*c1L-b1R*c1R)/(np.abs(eig0-eig1))**2.0
    
#    a=np.matmul(phi0L,HdiffL(TL,TR,chi))
#    b2L=np.matmul(a,phi2R)
#    a=np.matmul(phi2L,HdiffR(TL,TR,chi))
#    c2L=np.matmul(a,phi0R)
    
#    a=np.matmul(phi0L,HdiffR(TL,TR,chi))
#    b2R=np.matmul(a,phi2R)
#    a=np.matmul(phi2L,HdiffL(TL,TR,chi))
#    c2R=np.matmul(a,phi0R)
    
#    element2=(b2L*c2L-b2R*c2R)/(np.abs(eig0-eig2))**2.0
    
#    intg=element1+element2
    intg=element1

    return intg

def BerryCurvature(TL,TR):
#    BC=(-1j)*(integrand(muL,muR,dchi/2)-integrand(muL,muR,-dchi/2))/(dchi)
    #BC=(-1j)*2*F(muL,muR,dchi/2)/(dchi)
    temp1=F(TL,TR,dchi/2)
#    temp2=F(TL,TR,0,dchi/2)

    
    BC=2*(temp1.imag)/(dchi)
#    print("BC",BC_C)
    return BC



def BC_plot():
    Ngrid=21
    xlist = np.linspace(10, 500, Ngrid,endpoint=True)
    ylist = np.linspace(10, 500, Ngrid,endpoint=True)
    X, Y = np.meshgrid(xlist, ylist)
    BC=np.zeros((Ngrid,Ngrid),dtype=complex)

    for i in range(Ngrid):
        for j in range(Ngrid):
            BC[i,j]=BerryCurvature(X[0][i],Y[j][0])
            
    plt.figure()
    cp = plt.contourf(X, Y, BC,128)
    plt.colorbar(cp)
    #plt.title('Filled Contours Plot')
    plt.xlabel('$\mu_L$ (meV)')
    plt.ylabel('$\mu_R$ (meV)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
 
    
    return

BC_plot()






#print(dynamic(12e-3,8e-3,8e-3))

print("time:", (time.time() - start_time)/60.0)
