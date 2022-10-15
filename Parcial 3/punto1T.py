# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 18:49:03 2022

@author: rober
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator,FormatStrFormatter
import mpl_toolkits.mplot3d.axes3d

#parametros fisicos
k=0.49
q=1
c=10.0 # condicion de neumann

#parametros geometricos
H=1
W=1
dx=W/4
dy=H/4
Ta=37
Tb=36


n=10 #incognitas
A=np.zeros((n,n))
u=np.zeros(n)
B=np.zeros(n)

##Familias de nodos
F1=[5] #nodos internos
F2=[1,2,3] #nodos soporte
F3=[7,8,9] #nodos arriba
F4=[4] #nodo izquierda
F5=[6] #nodo derecha
F6=[0] #nodo abajo

for i in F1:
    A[i,i+1]=1/dx**2
    A[i,i-1]=1/dx**2
    A[i,i+3]=1/dy**2
    A[i,i-2]=1/dy**2
    A[i,i]=-2/dx**2-2/dy**2
    B[i]=-q/k
    
for i in F2:
    A[i,i]=1
    B[i]=(Ta+Tb)/2

for i in F3:
    A[i,i-3]=1/dy
    A[i,i]=-1/dy
    B[i]=c

for i in F4:
    A[i,i]=1
    B[i]=Ta

for i in F5:
    A[i,i]=1
    B[i]=Tb
    
for i in F6:
    A[i,i+1]=1/dy
    A[i,i]=-1/dy
    B[i]=c

u=np.linalg.solve(A,B)
T=np.zeros((6,3))

T[2:,0]=Ta
T[2:,2]=Tb
T[0,:]=u[7:10]
T[5,1]=u[0]
T[1,0]=u[4]
T[1,1]=u[5]
T[1,2]=u[6]
T[2:5,1]=u[1:4]

T1=np.copy(T)
plt.figure(1)
X=np.arange(3)
Y=np.arange(6)
coordX,coordY = np.meshgrid(X,Y)
T1[:]=T[::-1,:]
plt.pcolormesh(coordX,coordY,T1,cmap="plasma", shading=("gouraud"))
plt.colorbar()

plt.figure(2)
CS=plt.contour(coordX,coordY,T1,cmap="plasma", shading=("gouraud"))

plt.figure(3)
fig=plt.figure()
ax=fig.gca(projection="3d")
surf=ax.plot_surface(coordX,coordY,T1,cmap="plasma")
