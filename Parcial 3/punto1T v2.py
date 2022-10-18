# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 22:59:35 2022

@author: Roberto Gomez
"""

import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator,FormatStrFormatter
# import mpl_toolkits.mplot3d.axes3d

#parametros fisicos
k=1e-3
q=1
c=10.0 # condicion de neumann

#parametros geometricos
H=1
W=1
dx=W/6
dy=H/13
Ta=35
Tb=30


n=58 #incognitas
A=np.zeros((n,n))
u=np.zeros(n)
B=np.zeros(n)

##Familias de nodos
F1=[4,7,10,13,16,19,22,25,28] #nodos internos 1
F2=[32,33,34] #nodos internos 2
F3=[38,39,40,41,42,45,46,47,48,49] #nodos internos 3
F4=[51,52,53,54,55,56,57] #nodos arriba
F5=[3,6,9,12,15,18,21,24,27,31,30,37,44] #nodos Ta
F6=[5,8,11,14,17,20,23,26,29,35,36,43,50] #nodos Tb
F7=[0,1,2] #nodos abajo

for i in F1:
    A[i,i+1]=1/dx**2
    A[i,i-1]=1/dx**2
    A[i,i+3]=1/dy**2
    A[i,i-3]=1/dy**2
    A[i,i]=-2/dx**2-2/dy**2
    B[i]=-q/k
    
for i in F2:
    A[i,i+1]=1/dx**2
    A[i,i-1]=1/dx**2
    A[i,i+7]=1/dy**2
    A[i,i-5]=1/dy**2
    A[i,i]=-2/dx**2-2/dy**2
    B[i]=-q/k

for i in F3:
    A[i,i+1]=1/dx**2
    A[i,i-1]=1/dx**2
    A[i,i+7]=1/dy**2
    A[i,i-7]=1/dy**2
    A[i,i]=-2/dx**2-2/dy**2
    B[i]=-q/k

for i in F4:
    A[i,i-7]=-1/dy
    A[i,i]=1/dy
    B[i]=c

for i in F5:
    A[i,i]=1
    B[i]=Ta

for i in F6:
    A[i,i]=1
    B[i]=Tb
    
for i in F7:
    A[i,i+3]=1/dy
    A[i,i]=-1/dy
    B[i]=c

u=np.linalg.solve(A,B)
T=np.zeros((14,7))

T[3:,0:2]=Ta
T[3:,5:7]=Tb
T[13,2:5]=u[0:3]
T[12,2:5]=u[3:6]
T[11,2:5]=u[6:9]
T[10,2:5]=u[9:12]
T[9,2:5]=u[12:15]
T[8,2:5]=u[15:18]
T[7,2:5]=u[18:21]
T[6,2:5]=u[21:24]
T[5,2:5]=u[24:27]
T[4,2:5]=u[27:30]
T[3,:]=u[30:37]
T[2,:]=u[37:44]
T[1,:]=u[44:51]
T[0,:]=u[51:58]


T1=np.copy(T)
plt.figure(1)
X=np.arange(7)
Y=np.arange(14)
coordX,coordY = np.meshgrid(X,Y)
T1[:]=T[::-1,:]
plt.pcolormesh(coordX,coordY,T1,cmap="plasma", shading=("gouraud"))
plt.title('Propagacion de calor en 2D')
plt.colorbar()

plt.figure(2) #Se define una figura
plt.title('Propagacion de calor en 2D - Lineas de contorno')
CS=plt.contour(coordX,coordY,T1,cmap='plasma')
plt.clabel(CS,inline=1,fontsize=10)
plt.colorbar()
plt.show()

plt.figure(3)
fig=plt.figure()
ax=fig.gca(projection='3d')
surf=ax.plot_surface(coordX,coordY,T1,cmap='plasma')
fig.colorbar(surf)
plt.show()
