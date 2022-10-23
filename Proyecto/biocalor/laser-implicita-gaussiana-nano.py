# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:15:12 2022

@author: rober
"""

#caso implicito

import numpy as np
import matplotlib.pyplot as plt

def pos(i,j,k):
    pos=i+n*j+n**2*k
    return pos

def pos1(N):
    k=N//n**2
    j=N%n**2//n
    i=N%n**2%n
    return i,j,k

def laser(x,y,z,ua):
    r=np.sqrt(x**2+y**2+(z+dermis)**2)
    p=ua*(pnir*np.exp(-ueff*r))/(4*np.pi*Doptic*r)
    return p

def arrhenius(T,te):
    omega = Ar*np.exp(-Ea/(R*T))*te
    return omega
    
    
    
    
#dimensiones de la malla
largo=10e-3 #Largo malla
ancho=10e-3 #Ancho malla
alto=10e-3#Alto malla

#datos condiciones de frontera e iniciales
Tcorp=37
Tambiente=24
h=10 #coeficiente convectivo

#potencia del laser
pnir=1 #W
#parametros del tumor
K=0.5 #conductividad
rho=1052 #densidad
cp=3800 #calor especifico
dermis=1e-3 #distancia inicial desde el laser
Qb=40000 #calor metabolico
#para arrhenius
Ar=5.94e91 #energia de activación
Ea=0.5867e6 #factor de frecuencia
R=8.314 #constante de gas
#parametros de sangre
rhob=1052
cpb=3800
wb=0.01 #perfusion de sangre
 #calor por perfusion sanguinea

#propiedades opticas del tumor
g=0.9 #factor anisotropico
ua=80 #indice de absorcion
us=10000 #indice de dispersion
usp=us*(1-g) #indice de dispersion reducido
uanano=100 #indice de absorcion de nanoparticulas
uspnano=0.1 #indice de dispersion de nanoparticulas
uat=ua+uanano #indice de absorcion total
uspt=usp+uspnano #indice de dispersion total
ueff=np.sqrt(3*ua*(ua+usp))
Doptic=1/(3*(usp+ua)) #densidad optica


#parametros de la malla y tiempo
n=5#numero nodos
dz=alto/(n-1)
dx=largo/(n-1)
dy=ancho/(n-1)
dt=5
niter=120


## M A I N  ##
solution=np.zeros(shape=(n,n,n,niter))
solution_nano=np.copy(solution)
tiempo=np.linspace(0,dt*niter/60,niter) #minutos
nodoz=np.linspace(0,n*dz*1e3,n)

#para gauss sediel
A=np.zeros((n**3,n**3))

B=np.zeros(n**3)
Bnano=np.zeros(n**3)

x1=np.ones(n**3)
x2=np.zeros(n**3)

x1nano=np.ones(n**3)
x2nano=np.zeros(n**3)

tolerancia=1e-10


solution[:,:,:,0]=Tcorp
solution_nano[:,:,:,0]=Tcorp

Tp=np.copy(solution[:,:,:,0])
Tf=np.zeros((n,n,n))

Tpnano=np.copy(Tp)
Tfnano=np.zeros((n,n,n))


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
                        A[pos(i,j,k),pos(i,j,k)]=-1/dz-h
                        A[pos(i,j,k),pos(i,j,k+1)]=1/dz
                        B[pos(i,j,k)]=h*(-Tambiente)
                        
                    # elif k==0:
                    #     A[pos(i,j,k),pos(i,j,k)]=-1/dz
                    #     A[pos(i,j,k),pos(i,j,k+1)]=1/dz
                    #     B[pos(i,j,k)]=-K*h*(Tcorp-Tambiente)
                                        
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
                        
                        A[pos(i,j,k),pos(i,j,k)]=-2*K/dx**2-2*K/dy**2-2*K/dz**2-rho*cp/dt-rhob*cpb*wb
                        # A[pos(i,j,k),pos(i,j,k)]=-2*K/dx**2-2*K/dy**2-2*K/dz**2-rho*cp/dt
                        
                        B[pos(i,j,k)]=-rho*cp/dt*Tp[i,j,k]-Qb-rhob*cpb*wb*Tcorp-laser(dx*i-largo/2,dy*j-ancho/2,dz*k,ua)
                        Bnano[pos(i,j,k)]=-rho*cp/dt*Tpnano[i,j,k]-Qb-rhob*cpb*wb*Tcorp-laser(dx*i-largo/2,dy*j-ancho/2,dz*k,uat)
                        # B[pos(i,j,k)]=-rho*cp/dt*Tp[i,j,k]-Qb-rhob*cpb*wb*Tcorp
                        # Bnano[pos(i,j,k)]=-rho*cp/dt*Tpnano[i,j,k]-Qb-rhob*cpb*wb*Tcorp
        
    # GAUSS SEIDEL 1
    D=np.diag(np.diag(A))
    L=np.tril(A,k=-1)
    U=np.triu(A,k=1)
    error = 1.0
    
    iter1=0
    while error>=tolerancia:
        cj=np.matmul(-np.linalg.inv(L+D),(U))        
        dj=np.matmul(np.linalg.inv(D+L),B)        
        x2=np.matmul(cj,x1)+dj        
        error=np.linalg.norm(x2-x1,2)
        iter1+=1        
        x1=x2
    
    # GAUSS SEDIEL 2
    iter2=0
    error = 1.0
    
    while error>=tolerancia:
        cj=np.matmul(-np.linalg.inv(L+D),(U))        
        dj=np.matmul(np.linalg.inv(D+L),Bnano)        
        x2nano=np.matmul(cj,x1nano)+dj     
        error=np.linalg.norm(x2nano-x1nano,2)
        iter2+=1        
        x1nano=x2nano
        
    for z in range(n**3):
        solution[pos1(z)[0],pos1(z)[1],pos1(z)[2],t]=x2[z]
        Tf[pos1(z)[0],pos1(z)[1],pos1(z)[2]]=x2[z]
        
        solution_nano[pos1(z)[0],pos1(z)[1],pos1(z)[2],t]=x2nano[z]
        Tfnano[pos1(z)[0],pos1(z)[1],pos1(z)[2]]=x2nano[z]
    
    print(f'{round(t*100/niter,2)}% - Gauss Seidel: {iter1} iteraciones.')

    # solution[:,:,:,j+1]=Tf
    Tp=Tf
    Tpnano=Tfnano
    
#ARRHENIUS
omega=np.zeros(shape=(n,n,n,niter))
omega_nano=np.copy(omega)

for t in range(1,niter):
    for k in range(0,n):
        for j in range(0,n):
            for i in range(0,n):
                omega[i,j,k,t]=omega[i,j,k,t-1]+arrhenius(solution[i,j,k,t]+273.15,dt)
                omega_nano[i,j,k,t]=omega_nano[i,j,k,t-1]+arrhenius(solution_nano[i,j,k,t]+273.15,dt)



#plots
X=np.linspace(0,largo*1e3,n)
Y=np.linspace(-ancho*1e3/2,ancho*1e3/2,n)
coordX,coordY = np.meshgrid(X,Y)

plt.figure(1)
T=np.copy(solution[:,n//2+1,:,niter-1])
T[:]=T[::-1,::-1]
plt.pcolormesh(coordY,coordX,T,cmap="plasma", shading=("gouraud"))
plt.colorbar()
# plt.clim(max(x1))
plt.xlabel('mm')
plt.ylabel('mm')
plt.title(f'Laser - t= {dt*niter/60} min')

plt.figure(2)
T[:]=np.copy(solution_nano[:,n//2+1,:,niter-1])
T[:]=T[::-1,::-1]
plt.pcolormesh(coordY,coordX,T,cmap="plasma", shading=("gouraud"))
plt.colorbar()
# plt.clim(max(x2)) 
plt.xlabel('mm')
plt.ylabel('mm')
plt.title(f'Laser + NPs - t= {dt*niter/60} min')

plt.figure(3) #Se define una figura
T=np.copy(omega[:,n//2,:,niter-1])
T[:]=T[::-1,::-1]
plt.title(f'Muerte celular: Laser t= {dt*niter/60} min')
CS=plt.contour(coordY,coordX,T,cmap='plasma')
plt.clabel(CS,inline=1,fontsize=10)
CS.cmap.set_over('red')
CS.cmap.set_under('blue')
plt.clim(0,1)
plt.colorbar()
plt.show()

plt.figure(4) #Se define una figura
T=np.copy(omega_nano[:,n//2,:,niter-1])
T[:]=T[::-1,::-1]
plt.title(f'Muerte celular: Laser + NPs t = {dt*niter/60} min')
CS=plt.contour(coordY,coordX,T,cmap='plasma')
plt.clabel(CS,inline=1,fontsize=10)
CS.cmap.set_over('red')
CS.cmap.set_under('blue')
plt.clim(0,1)
plt.colorbar()
plt.show()


plt.figure(5)
plt.plot(nodoz,solution[n//3,n//3,:,niter-1], label="Laser")
plt.plot(nodoz,solution_nano[n//3,n//3,:,niter-1], label=f"Laser + NPs - t = {dt*niter/60} min")
plt.xlabel('mm')
plt.ylabel('Temperatura (Celsius)')
plt.legend(loc='best') 


plt.figure(6)
plt.plot(tiempo,solution[3,3,3,:], label="Laser")
plt.plot(tiempo,solution_nano[3,3,3,:], label="Laser + Nanoparticulas")
plt.ylabel('Temperatura (Celsius)')
plt.xlabel('Tiempo (min)')
plt.title('Calentamiento en un nodo')
plt.legend(loc='best') 