# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:15:12 2022

@author: rober
"""

#caso explicito

import numpy as np
import matplotlib.pyplot as plt

tolerancia=1e-10
error=1.0


#condiciones iniciales
largo=5e-1 #Largo malla
ancho=5e-1 #Ancho malla
alto=5e-1 #Alto malla
Tcorp=0
pnir=1.5
K=0.19
rho=1000
cp=3900

#propiedades opticas
g=0.9
ua=2.2
us=1220
usp=us*(1-g)
ueff=np.sqrt(3*ua*(ua+usp))
Doptic=1/(3*(usp+ua))


#parametros de la malla y tiempo
n=8 #numero nodos
dz=alto/(n-1)
dx=largo/(n-1)
dy=ancho/(n-1)
dt=1
niter=100

solution=np.zeros(shape=(n,n,n,niter))
tiempo=np.linspace(0,niter,niter)
nodo=np.linspace(0,n,n)
A=np.zeros((n**3,n**3))
B=np.zeros(n**3)
x1=np.ones(n**3)
x2=np.zeros(n**3)


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
    p=ua*((pnir*np.exp(-ueff*r))/4*np.pi*Doptic*r)
    return p

def pos(i,j,k):
    pos=i+n*j+n**2*k
    return pos

def pos1(N):
    k=N//n**2
    j=N%n**2//n
    i=N%n**2%n
    return i,j,k
er=0
iter=0
for t in range(1,len(solution[0,0,0,:])):
    for k in range(0,len(solution[:,0,0,0])):
        for j in range(0,len(solution[0,:,0,0])):
            for i in range(0,len(solution[0,0,:,0])):
            
                    if k==0 or i==0 or j==0:
                        A[pos(i,j,k),pos(i,j,k)]=1
                        B[pos(i,j,k)]=Tcorp

                    elif k==n-1 or j==n-1 or i==n-1:
                        A[pos(i,j,k),pos(i,j,k)]=1
                        B[pos(i,j,k)]=Tcorp
                        # print('debug')
                    
                    elif (0<i<n-1) and (0<j<n-1) and (0<k<n-1):
                        
                        A[pos(i,j,k),pos(i+1,j,k)]=K/dx**2
                        A[pos(i,j,k),pos(i-1,j,k)]=K/dx**2
                        
                        A[pos(i,j,k),pos(i,j+1,k)]=K/dy**2
                        A[pos(i,j,k),pos(i,j-1,k)]=K/dy**2
                        
                        A[pos(i,j,k),pos(i,j,k+1)]=K/dz**2
                        A[pos(i,j,k),pos(i,j,k-1)]=K/dz**2
                        
                        A[pos(i,j,k),pos(i,j,k)]=-K*(2/dx**2+2/dy**2+2/dz**2)-rho*cp/dt
                        B[pos(i,j,k)]=-rho*cp/dt*solution[i,j,k,t-1]-laser(dx*i-largo/2,dy*j-ancho/2,dz*k)
                        # B[pos(i,j,k)]=-rho*cp/dt*solution[i,j,k,t-1]
                    
    
    # x2=np.linalg.solve(A,B)
        
    #GAUSS SEDIEL
    np.diag(A)
    D=np.diag(np.diag(A))
    L=np.tril(A,k=-1)
    U=np.triu(A,k=1)
    while error>=tolerancia:
        cj=np.matmul(-np.linalg.inv(L+D),(U))
        
        dj=np.matmul(np.linalg.inv(D+L),B)
        
        x2=np.matmul(cj,x1)+dj
        
        error=np.linalg.norm(x2-x1,2)  #este sera el error que tenemos con la euclidea
        iter+=1
        
        x1=x2
        
    for z in range(n**3):
        solution[pos1(z)[0],pos1(z)[1],pos1(z)[2],t]=x2[z]
    print(t)
    




for i in range(0,niter,niter//10):
    plt.figure(i)
    X=np.arange(n)
    Y=np.arange(n)
    coordX,coordY = np.meshgrid(X,Y)
    T=solution[:,n//2,:,i]
    T[:]=T[::-1,::-1]
    plt.pcolormesh(coordY,coordX,T,cmap="plasma", shading=("gouraud"))
    plt.colorbar()

plt.figure(niter)
plt.plot(nodo,solution[:,n//3,n//5,niter-1])
plt.figure(niter+1)
plt.plot(tiempo,solution[n//3,n//3,n//5,:])