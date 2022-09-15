
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:56:47 2021

@author: usuario
"""

#librerias
import numpy as np
import matplotlib.pyplot as plt
from math import cos


#funcion a resolver

def fode(x,y):
    dYdX=np.array([y[1],cos(2*x)-y[0]])
    return dYdX

#metodo euler 
def meuler(xi,Yi,h,F):
    x=xi+h
    Y=Yi+h*F(xi,Yi)
    
    return x,Y

def rk4(xi,Yi,h,f):
    k1=h*f(xi,Yi)
    k2=h*f(xi+0.5*h,Yi+0.5*k1)
    k3=h*f(xi+0.5*h,Yi+0.5*k2)
    k4=h*f(xi+h,Yi+k3)
    Y=Yi+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4
    x=xi+h
    return x,Y

#*****************************************************
#main del programa


#datos iniciales
xi=0.0
xf=3.0
n=20# numero de datos a analizar
h=(xf-xi)/(n-1) #paso




#definiendo el tamaño de los arreglos
X=np.zeros(n)

YnumEuler=np.zeros(shape=(n,2))

YnumRK=np.zeros(shape=(n,2))



#condicion inicial
Yi=np.array([1.0,0.0])

YnumEuler[0,:]=Yi  #Euler

YnumRK[0,:]=Yi   #Runge Kutta4

# para Euler
for i in range(0,n-1,1):
    X[i+1],YnumEuler[i+1,:]=meuler(X[i],YnumEuler[i,:],h,fode)
   

# para Runge Kutta 4    

for i in range(0,n-1,1):
    X[i+1],YnumRK[i+1,:]=rk4(X[i],YnumRK[i,:],h,fode)



#Graficos de la solución

plt.figure(1)

plt.plot(X,YnumEuler[:,0],color="b",label="Euler:Y vs X") #grafico para Euler

plt.plot(X,YnumRK[:,0],color="c",label="RK4:Y vs X") #grafico para Runge Kutta 4



plt.legend()

#############################

plt.figure(2)

plt.plot(X,YnumEuler[:,1],color="r",label="Euler:dY/dX vs X") #grafico para Euler

plt.plot(X,YnumRK[:,1],color="m",label="RK4:dY/dX vs X")  #grafico para Runge Kutta 4

plt.legend()





















    