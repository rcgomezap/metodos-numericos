#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 12:47:40 2021

@author: usuario
"""

#caso esquema implicito del tiempo 
#ejemplo pared en una dimension

import numpy as np
import matplotlib.pyplot as plt


#parametros fisicos
k=200. #conductividad del material
rho=2000.  #densidad del material
cp=800. #capacidad calorifica del material
L=0.4 #espesor de pared
Tw=20. #temperatura de la pared
Tini=100. #temperatura inicial de los nodos internos


#parametros de la malla
n=6
dx=L/(n-1)
dt=1.
niter=500

#creacion de arreglos
xnodes=np.linspace(0,L,n)  #vector de las coordenadas x de cada nodo
tnodes=np.linspace(0,niter*dt,niter+1)


A=np.zeros((n,n))
B=np.zeros(n)

Tp=np.zeros(n)
Tf=np.zeros(n)
solution=np.zeros((n,niter+1))

#condiciones iniciales
Tp[1:n-1:1]=Tini
Tp[0]=Tw
Tp[n-1]=Tw
solution[:,0]=Tp

iter=0


for j in range(0,niter): #ciclo del tiempo
    #ensamble matriz A
    for i in range(0,n): 
        if i==0:
            A[i,i]=1
            B[i]=Tw
        
        if i>0 and i<n-1:
            A[i,i-1]=k/dx**2
            A[i,i+1]=k/dx**2
            A[i,i]=-2*k/dx**2-rho*cp/dt
            B[i]=-rho*cp/dt*Tp[i]
        
        if i==n-1:
            A[i,i]=1
            B[i]=Tw    
            
    Tf=np.linalg.solve(A,B)  #solucion del sistema AX=B
    solution[:,j+1]=Tf[:] #agrega a columna de matriz Historial
    Tp=Tf # para la siguiente iteracion se actualiza Tp
    

plt.figure(1)
plt.plot(xnodes,solution[:,0],label="j=0")
plt.plot(xnodes,solution[:,30],label="j=30")
plt.plot(xnodes,solution[:,100],label="j=100")
plt.plot(xnodes,solution[:,300],label="j=300")
plt.plot(xnodes,solution[:,500],label="j=500")
plt.title('solucion implicita')
plt.legend(loc='upper right')

plt.figure(2)
plt.plot(tnodes,solution[1,:],label="nodo1")
#plt.plot(tnodes,solution[2,:],label="nodo2")   
plt.legend(loc='upper right') 
    
    
    
    
    
    
    
        
    
































