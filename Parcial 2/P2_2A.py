# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:19:57 2022

@authors:   Roberto Carlos Gomez Araque: 000423542
            Juan Camilo Fernandez Angarita: 000424705
"""
import numpy as np
import matplotlib.pyplot as plt

#Datos iniciales
b=0.3
c=0.1
m=1
g=9.8
u=0.6

N=m*g


def fode(t,x):
    k=1600
    n=0.1
    v=v0-x[1]
    F=u*N*(1-b*v+c*(v**3))
    F1=k*x[0]
    F2=n*x[1]
    dxdt=np.array([x[1],(F-F1-F2)*m])
    return dxdt


def rk4(xi,Yi,h,f):
    k1=h*f(xi,Yi)
    k2=h*f(xi+0.5*h,Yi+0.5*k1)
    k3=h*f(xi+0.5*h,Yi+0.5*k2)
    k4=h*f(xi+h,Yi+k3)
    Y=Yi+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4
    x=xi+h
    return x,Y

#datos iniciales
ti=0
tf=10
vi=np.array([0.5,1,1.5])
xi=np.array([0,0.5])
n=2000# numero de datos a analizar
h=(tf-ti)/(n-1) #paso


t=np.zeros(n)
x1=np.zeros(shape=(n,2))
x2=np.zeros(shape=(n,2))
x3=np.zeros(shape=(n,2))



#condicion inicial

x1[0,:]=xi
x2[0,:]=xi 
x3[0,:]=xi 

for i in range(0,n-1,1):
    v0=vi[0]
    t[i+1],x1[i+1,:]=rk4(t[i],x1[i,:],h,fode)
    v0=vi[1]
    t[i+1],x2[i+1,:]=rk4(t[i],x2[i,:],h,fode)
    v0=vi[2]
    t[i+1],x3[i+1,:]=rk4(t[i],x3[i,:],h,fode)
    
##Determinacion de la estabilidad del sistema
    

plt.figure(1)
plt.figure(figsize=(15,5))

plt.subplot(2,2,1)
plt.plot(t,x1[:,0],color='r',label='v0=0.5')
plt.legend()
plt.grid()

plt.subplot(2,2,3)
plt.plot(t,x2[:,0],color='r',label='v0=1.0')
plt.legend()
plt.grid()

plt.subplot(2,2,2)
plt.title('0.5')
plt.plot(t,x3[:,0],color='r',label='v0=1.5')
plt.legend()
plt.grid()