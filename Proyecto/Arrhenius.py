# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 20:19:04 2022

@author: rober
"""
#librerias
import numpy as np
import matplotlib.pyplot as plt


#funciones

#EDO a resolver
def arrhenius(x,y):
    dydx=2*np.exp(1/x)
    return dydx

def rk4(xi,Yi,h,f):
    k1=h*f(xi,Yi)
    k2=h*f(xi+0.5*h,Yi+0.5*k1)
    k3=h*f(xi+0.5*h,Yi+0.5*k2)
    k4=h*f(xi+h,Yi+k3)
    Y=Yi+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4
    x=xi+h
    return x,Y

xi=1.0
yi=0
xf=2
n=2
h=(xf-xi)/(n-1) #paso
x=np.zeros(n)
ynum1=np.zeros(n)

#valores iniciales
x[0]=xi
ynum1[0]=yi

#solucion runge kutta
for i in range(1,n):
    x[i],ynum1[i]=rk4(x[i-1],ynum1[i-1],h,arrhenius)
    
#grafico
plt.figure(1)
plt.plot(x,ynum1,'*-',color='c',label='RK4')
plt.xlabel('x')
plt.ylabel('y')
plt.title('paso h:'+str(h))
plt.legend(loc='best')

