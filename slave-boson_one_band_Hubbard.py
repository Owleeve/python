#/usr/bin
import sympy
from sympy import diff
from sympy import symbols
from sympy import solve
#from sympy import exp

import numpy as np
from numpy import arange
from numpy import linspace
from numpy import pi

import math 
from math import cos
from math import exp

from scipy.optimize import fsolve


#define the constants:
t=1 # hopping factor
a=1 # lattice vector length
#beta=100 # beta=1/KT
Nx=45
Ny=45
N=Nx*Ny # number of electron or k in first BZ.
  
#get energy of hopping
Ek=np.zeros((2*N))
for x in range(0,Nx):
    for y in range(0,Ny):
        i = Ny*x + y
        Ek[2*i]= -2*t*(cos(2*pi*x/Nx)+cos(2*pi*y/Ny))
        Ek[2*i+1] =Ek[2*i]
Ek=np.sort(Ek)

#average energy per site
def AveE(num):
    energy=0
    for i in range(0,num):
        energy=energy+1/N*Ek[i]
    return energy

def q(e,p,d):
    return (e+d)**2*p**2/((1-d**2-p**2)*(1-e**2-p**2))

# iteration cri
def norm(x_m,x_n):
    err=0
    for i in range(0,3):
        z=(x_m[i]-x_n[i])/x_m[i]
        err=err+abs(z)
    return err
#
accurate=1e-3

u_list=[0.5,1.0001,5]

for u in u_list:
    for delta in linspace(0.21,0.97,20):
        Ne = int(N*(1-delta))
        U=8*u*abs(AveE(Ne))
       
        def functions(unsolved):
            e=unsolved[0]
            p=unsolved[1]
            d=unsolved[2]
            lam=unsolved[3]
            lam_=unsolved[4]
            f1=2*p**2+e**2+d**2-1
            f2=Ne/N-2*(p**2+d**2)
            #f3=AveE(Ne)*q_e(e,p,d)+2*e*lam
            f3=AveE(Ne)*(2*e*p**2*(d + e)**2/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)**2) + p**2*(2*d + 2*e)/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)))+2*e*lam
            #f3=-2*(2*e*p**2*(d + e)**2/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)**2) + p**2*(2*d + 2*e)/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)))+2*e*lam
            #f4=AveE(Ne)*q_d(e,p,d)+2*U*d-4*lam_*d+2*lam*d
            f4=AveE(Ne)*(2*d*p**2*(d + e)**2/((-d**2 - p**2 + 1)**2*(-e**2 - p**2 + 1)) + p**2*(2*d + 2*e)/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)))+2*U*d-4*lam_*d+2*lam*d
            #f4=-2*(2*d*p**2*(d + e)**2/((-d**2 - p**2 + 1)**2*(-e**2 - p**2 + 1)) + p**2*(2*d + 2*e)/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)))+2*U*d-4*lam_*d+2*lam*d
            #f5=AveE(Ne)*q_p(e,p,d)-2*lam_*p+2*lam*p
            f5=AveE(Ne)*(2*p**3*(d + e)**2/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)**2) + 2*p**3*(d + e)**2/((-d**2 - p**2 + 1)**2*(-e**2 - p**2 + 1)) + 2*p*(d + e)**2/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)))-4*lam_*p+4*lam*p
            #f5=-2*(2*p**3*(d + e)**2/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)**2) + 2*p**3*(d + e)**2/((-d**2 - p**2 + 1)**2*(-e**2 - p**2 + 1)) + 2*p*(d + e)**2/((-d**2 - p**2 + 1)*(-e**2 - p**2 + 1)))-2*lam_*p+2*lam*p
            #f6=1-2*(d**2+p**2)-delta
            #f=np.array([f1,f2,f3,f4,f5,f6])
            #f=np.array([f1,f2,f3,f4,f5])
            return [f1,f2,f3,f4,f5]

        xx=[0.5,0.5,0.5,0.5,0.5]
        solved=xx

        item=0
        s=[solved]
        solved=fsolve(functions,solved)
        s.append(solved)
        item +=1
        while norm(s[item],s[item-1]) >= accurate:
            solved=fsolve(functions,solved)
            s.append(solved)
            item +=1

        e=solved[0]
        p=solved[1]
        d=solved[2]
        q_=q(e,p,d)
        outputfile_delta_q='./u_{}_delta-q.dat'.format(u)  # - cannot be used in python
        outputf_delta_q=open(outputfile_delta_q,"a")
        outputf_delta_q.write(str(delta)+" "+str(q_)+"\n")
        


