# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 12:22:37 2022

@author: B09S114est
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
Ta=75.0
Tb=100.0

dim=5 #datos por fila
n=dim*dim #incognitas
A=np.zeros((n,n))
u=np.zeros(n)
B=np.zeros(n)

##Familias de nodos
F1=[6,7,8,11,12,13,16,17,18] #nodos internos
F2=[0,4,5,9,10,14,15,19] #nodos laterales
F3=[20,21,22,23,24] #nodos arriba
F4=[1,2,3] #nodos abajo

for i in F1:
    A[i,i+1]=1/dx**2
    A[i,i-1]=1/dx**2
    A[i,i+5]=1/dy**2
    A[i,i-5]=1/dy**2
    A[i,i]=-2/dx**2-2/dy**2
    B[i]=-q/k

for i in F2:
    A[i,i]=1
    B[i]=Ta

for i in F3:
    A[i,i]=1
    B[i]=Tb

for i in F4:
    A[i,i+5]=1/dy
    A[i,i]=-1/dy
    B[i]=c
    
u=np.linalg.solve(A,B)
T=np.zeros((dim,dim))

T[0,:]=u[20:25]
T[1,:]=u[15:20]
T[2,:]=u[10:15]
T[3,:]=u[5:10]
T[4,:]=u[0:5]

X=np.arange(dim)
Y=np.arange(dim)
coordX,coordY = np.meshgrid(X,Y)
T[:]=T[::-1,:]
plt.pcolormesh(coordX,coordY,T,cmap="plasma", shading=("gouraud"))
plt.colorbar()