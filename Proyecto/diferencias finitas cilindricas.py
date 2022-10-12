# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:15:12 2022

@author: rober
"""

#caso explicito

import numpy as np
import matplotlib.pyplot as plt



#condiciones iniciales
L=0.5 #Largo cilindro
A=0.5 #Ancho cilindro
H=0.5 #Alto cilindro
Tcorp=37.0
pnir=1
k=200.0
rho=2000.0
cp=800-0


#parametros de la malla y tiempo
n=10 #numero nodos
dz=H/(n-1)
dx=L/(n-1)
dy=A/(n-1)
dt=1e-5
niter=800

solution=np.zeros(shape=(n,n,n,niter))
tiempo=np.linspace(0,niter,niter)
nodo=np.linspace(0,n,n)
# nodes=np.zeros(shape=dz)
solution[:,:,:,0]=Tcorp
solution[0,:,:,:]=Tcorp
solution[:,0,:,:]=Tcorp
solution[:,:,0,:]=Tcorp
solution[:,:,:,0]=Tcorp

solution[n-1,:,:,:]=Tcorp
solution[:,n-1,:,:]=Tcorp
solution[:,:,n-1,:]=Tcorp
solution[:,:,:,n-1]=Tcorp

##DISCRETIZACIONES

def D2Tcent(i,j,k,d): ##Segunda derivada centrada
    d2tx=(k-2*j+i)/dx**2
    return d2tx

def D2Tdcenti(i,j,k,d): ##Segunda derivada descentrad izquierda
    d2tx=(k-2*j+i)/dx
    return d2tx

def laser(x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    p=np.exp(-r)
    return p

for t in range(1,len(solution[0,0,0,:])):
    for i in range(0,len(solution[:,0,0,0])):
        for j in range(0,len(solution[0,:,0,0])):
            for k in range(0,len(solution[0,0,:,0])):
                if (0<i<n-1) and (0<j<n-1) and (0<k<n-1):
                    a=D2Tcent(solution[i-1,j,k,t-1],solution[i,j,k,t-1],solution[i+1,j,k,t-1],dx)
                    b=D2Tcent(solution[i,j-1,k,t-1],solution[i,j-1,k,t-1],solution[i,j+1,k,t-1],dy)
                    c=D2Tcent(solution[i,j,k-1,t-1],solution[i,j,k,t-1],solution[i,j,k+1,t-1],dz)
                    print(t)
                    
                    # solution[i,j,k,t]=(k*dt*((a+b+c)+laser(i-n/2,j-n/2,k)))+solution[i,j,k,t-1]
                    solution[i,j,k,t]=(k*dt*((a+b+c)+laser(i-n/2,j-n/2,k-n/2)))+solution[i,j,k,t-1]
                    



for i in range(0,n):
    plt.figure(i)
    X=np.arange(n)
    Y=np.arange(n)
    coordX,coordY = np.meshgrid(X,Y)
    T=solution[:,:,n//2,i]
    plt.pcolormesh(coordX,coordY,T,cmap="plasma", shading=("gouraud"))
    plt.colorbar()

plt.figure(n+1)
plt.plot(nodo,solution[:,n//2,n//2,n-1])
plt.figure(n+2)
plt.plot(tiempo,solution[n//2,n//2,n//2,:])