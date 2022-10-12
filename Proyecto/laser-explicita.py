# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:15:12 2022

@author: rober
"""

#caso explicito

import numpy as np
import matplotlib.pyplot as plt



#condiciones iniciales
L=10e-1 #Largo malla
A=10e-1 #Ancho malla
H=10e-1 #Alto malla
Tcorp=0
pnir=1500.0
K=0.19
rho=1000
cp=3900

#propiedades opticas
g=0.9
ua=2.2
us=1220
usp=us*(1-g)
ueff=np.sqrt(3*ua*(ua+usp))
D=1/(3*(usp+ua))


#parametros de la malla y tiempo
n=10 #numero nodos
dz=H/(n-1)
dx=L/(n-1)
dy=A/(n-1)
dt=1e-5
niter=1000

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

# def D2Tdcenti(i,j,k,d): ##Segunda derivada descentrad izquierda
#     d2tx=(k-2*j+i)/dx
#     return d2tx

def laser(x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    p=ua*((pnir*np.exp(-ueff*r))/4*np.pi*D*r)
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
                    
                    solution[i,j,k,t]=dt*(K*(a+b+c))/(rho*cp)+laser((dx*(i-n/2)),dy*(j-n/2),dz*k)+solution[i,j,k,t-1]
                    # solution[i,j,k,t]=((1/(rho*cp))*dt*(K*(a+b+c)+(1/rho*cp)*laser((i-n/2),(j-n/2),k)))+solution[i,j,k,t-1]
                    



for i in range(0,niter,niter//10):
    plt.figure(i)
    X=np.arange(n)
    Y=np.arange(n)
    coordX,coordY = np.meshgrid(X,Y)
    T=solution[:,n//2,:,i]
    T[:]=T[::-1,::-1]
    plt.pcolormesh(coordY,coordX,T,cmap="plasma", shading=("gouraud"))
    plt.colorbar()

plt.figure(n+1)
plt.plot(nodo,solution[:,n//2,n//2,n-1])
plt.figure(n+2)
plt.plot(tiempo,solution[n//2,n//2,n//2,:])