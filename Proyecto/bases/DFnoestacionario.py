# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:26:44 2022

@author: B09S114est
"""

#caso explicito

import numpy as np
import matplotlib.pyplot as plt

#parametros
k=200.0
rho=2000.0
cp=800-0
L=0.4

Tw=20.0
Tini=100.0

#parametros de la malla y tiempo
n=20 #numero nodos
dx=L/(n-1)
dt=1.0


niter=500 #numero de iteraciones
solution=np.zeros(shape=(n,niter))


#definicion de vectores
xnodes=np.linspace(0,L,n)
tnodes=np.linspace(0,niter*dt,niter+1)

Tpresente=np.zeros(n)
Tfuturo=np.zeros(n)

Tpresente[0]=Tw
Tpresente[-1]=Tw

Tpresente[1:-1]=Tini

for j in range(1,niter+1): #ciclo del tiempo
    for i in range(0,n):
        if i==0:
            Tfuturo[i]=Tpresente[i]
        elif i==n-1:
            Tfuturo[i]=Tpresente[i]
        else:
            Tfuturo[i]=k*dt/(rho*cp)*(Tpresente[i+1]-2*Tpresente[i]+Tpresente[i-1])/dx**2+Tpresente[i]
    
    solution[:,j-1]=Tfuturo
    # plt.figure(j) 
    # plt.plot(xnodes,Tfuturo)        
    # plt.ylim(0,100)
    Tpresente=Tfuturo
plt.plot(xnodes,solution[:,0])