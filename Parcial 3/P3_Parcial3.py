# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#caso explicito

import numpy as np
import matplotlib.pyplot as plt

#parametros
B1=5e-1
B2=10e-1
sigma=0.01

#parametros de la malla y tiempo
L=2
n=30 #numero nodos
dx=L/(n-1)
dt=1e-3
niter=500 #numero de iteraciones

solution1=np.zeros(shape=(n,niter))
solution2=np.zeros(shape=(n,niter))

def ini(x):
    nod=(1/(np.sqrt(2*np.pi)*sigma))*np.exp(-x**2/2*sigma**2)
    return nod

#definicion de vectores
xnodes=np.linspace(-L/2,L/2,n)
tnodes=np.linspace(0,niter*dt,niter)

upresente=np.zeros(n)
ufuturo=np.zeros(n)
upresente2=np.zeros(n)
ufuturo2=np.zeros(n)

for i in range(n): #condicion inicial u(x,0)
    upresente[i]=ini(dx*i-L/2)
    upresente2[i]=ini(dx*i-L/2)
    


for j in range(1,niter+1): #ciclo del tiempo
    for i in range(0,n):
        if 0<i<n-1:
            ufuturo[i]=B1*dt*(upresente[i+1]-2*upresente[i]+upresente[i-1])/dx**2+upresente[i]
            ufuturo2[i]=B2*dt*(upresente2[i+1]-2*upresente2[i]+upresente2[i-1])/dx**2+upresente2[i]
        elif i==0:
            ufuturo[i]=ufuturo[i+1]
            ufuturo2[i]=ufuturo2[i+1]
            
        elif i==n-1:
            ufuturo[i]=ufuturo[i-1]
            ufuturo2[i]=ufuturo2[i-1]
    
    solution1[:,j-1]=upresente
    solution2[:,j-1]=upresente2
    upresente=ufuturo
    upresente2=ufuturo2


plt.figure(1)
plt.plot(xnodes,solution1[:,0],color='r',label=str(0*dt)+" Segundos")
plt.plot(xnodes,solution1[:,30],color='b',label=str(30*dt)+" Segundos")
plt.plot(xnodes,solution1[:,100],color='g',label=str(100*dt)+" Segundos")
plt.plot(xnodes,solution1[:,300],color='y',label=str(300*dt)+" Segundos")
plt.plot(xnodes,solution1[:,niter-1],color='c',label=str(500*dt)+" Segundos")
plt.title('Solucion explicita. B=5e-1')
plt.legend(loc='best')

plt.figure(2)
plt.plot(xnodes,solution2[:,0],color='r',label=str(0*dt)+" Segundos")
plt.plot(xnodes,solution2[:,30],color='b',label=str(30*dt)+" Segundos")
plt.plot(xnodes,solution2[:,100],color='g',label=str(100*dt)+" Segundos")
plt.plot(xnodes,solution2[:,300],color='y',label=str(300*dt)+" Segundos")
plt.plot(xnodes,solution2[:,niter-1],color='c',label=str(500*dt)+" Segundos")
plt.title('Solucion explicita. B=10e-1')
plt.legend(loc='best')
