# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:15:12 2022

@author: rober
"""

#caso explicito

import numpy as np
import matplotlib.pyplot as plt

def laser(x,y,z):
    r=np.sqrt(x**2+y**2+(z+dermis)**2)
    p=ua*(pnir*np.exp(-ueff*r))/(4*np.pi*Doptic*r)
    return p

def pos(i,j,k):
    pos=i+n*j+n**2*k
    return pos

def pos1(N):
    k=N//n**2
    j=N%n**2//n
    i=N%n**2%n
    return i,j,k

tolerancia=1e-10


#condiciones iniciales
largo=3e-2 #Largo malla
ancho=3e-2 #Ancho malla
alto=2e-2#Alto malla

Tcorp=0
Tambiente=-5
h=10 #coeficiente convectivo

pnir=1.5

K=0.19
rho=1000
cp=3800
dermis=1e-4

#propiedades opticas
g=0.9
ua=0.8+2.65
# us=1220
usp=1220
ueff=np.sqrt(3*ua*(ua+usp))
Doptic=1/(3*(usp+ua))


#parametros de la malla y tiempo
n=7#numero nodos
dz=alto/(n-1)
dx=largo/(n-1)
dy=ancho/(n-1)
dt=5
niter=90

solution=np.zeros(shape=(n,n,n,niter))
tiempo=np.linspace(0,dt*niter,niter)
nodo=np.linspace(0,n,n)
A=np.zeros((n**3,n**3))
B=np.zeros(n**3)
x1=np.ones(n**3)
x2=np.zeros(n**3)


# nodes=np.zeros(shape=dz)
solution[:,:,:,0]=Tcorp
# solution[1:4,1:4,1:4,0]=Tcorp+20
Tp=np.zeros((n,n,n))
Tp=np.copy(solution[:,:,:,0])
Tf=np.zeros((n,n,n))



iter=0

for t in range(1,niter):
    for k in range(0,n):
        for j in range(0,n):
            for i in range(0,n):
                
                    if i==0:
                        A[pos(i,j,k),pos(i,j,k)]=-1/dx
                        A[pos(i,j,k),pos(i+1,j,k)]=1/dx
                        B[pos(i,j,k)]=0
                    
                    elif i==n-1:
                        A[pos(i,j,k),pos(i,j,k)]=1/dx
                        A[pos(i,j,k),pos(i-1,j,k)]=-1/dx
                        B[pos(i,j,k)]=0
                    elif j==0:
                        A[pos(i,j,k),pos(i,j,k)]=-1/dy
                        A[pos(i,j,k),pos(i,j+1,k)]=1/dy
                        B[pos(i,j,k)]=0
                    
                    elif j==n-1:
                        A[pos(i,j,k),pos(i,j,k)]=1/dy
                        A[pos(i,j,k),pos(i,j-1,k)]=-1/dy
                        B[pos(i,j,k)]=0
                    elif k==0:
                        A[pos(i,j,k),pos(i,j,k)]=K/dz-h
                        A[pos(i,j,k),pos(i,j,k+1)]=-K/dz
                        B[pos(i,j,k)]=h*Tambiente
                    
                    elif k==n-1:
                        A[pos(i,j,k),pos(i,j,k)]=1/dz
                        A[pos(i,j,k),pos(i,j,k-1)]=-1/dz
                        B[pos(i,j,k)]=0
                        
                    
                    elif (0<i<n-1) and (0<j<n-1) and (0<k<n-1):
                        
                        A[pos(i,j,k),pos(i+1,j,k)]=K/dx**2
                        A[pos(i,j,k),pos(i-1,j,k)]=K/dx**2
                        
                        A[pos(i,j,k),pos(i,j+1,k)]=K/dy**2
                        A[pos(i,j,k),pos(i,j-1,k)]=K/dy**2
                        
                        A[pos(i,j,k),pos(i,j,k+1)]=K/dz**2
                        A[pos(i,j,k),pos(i,j,k-1)]=K/dz**2
                        
                        A[pos(i,j,k),pos(i,j,k)]=-2*K/dx**2-2*K/dy**2-2*K/dz**2-rho*cp/dt
                        
                        B[pos(i,j,k)]=-rho*cp/dt*Tp[i,j,k]-laser(dx*i-largo/2,dy*j-ancho/2,dz*k)
                        # B[pos(i,j,k)]=-rho*cp/dt*Tp[i,j,k]
                        # print(B[pos(i,j,k)])
                    
    
    # x2=np.linalg.solve(A,B)
        
    # GAUSS SEDIEL
    D=np.diag(np.diag(A))
    L=np.tril(A,k=-1)
    U=np.triu(A,k=1)
    error = 1.0
    
    while error>=tolerancia:
        cj=np.matmul(-np.linalg.inv(L+D),(U))        
        dj=np.matmul(np.linalg.inv(D+L),B)        
        x2=np.matmul(cj,x1)+dj        
        error=np.linalg.norm(x2-x1,2)
        iter+=1        
        x1=x2
        
    for z in range(n**3):
        solution[pos1(z)[0],pos1(z)[1],pos1(z)[2],t]=x2[z]
        Tf[pos1(z)[0],pos1(z)[1],pos1(z)[2]]=x2[z]
    print(f'{t*100/niter}%')

    # solution[:,:,:,j+1]=Tf
    Tp=Tf
    




for z in range(0,niter,niter//10):
    plt.figure(z)
    X=np.arange(n)
    Y=np.arange(n)
    coordX,coordY = np.meshgrid(X,Y)
    T=solution[:,n//2,:,z]
    T[:]=T[::-1,::-1]
    plt.pcolormesh(coordY,coordX,T,cmap="plasma", shading=("gouraud"))
    plt.colorbar()
    # plt.clim(0,42)

plt.figure(niter)
plt.plot(nodo,solution[:,n//3,n//5,niter-1])
plt.figure(niter+1)
plt.plot(tiempo,solution[n//3,n//3,n//5,:])