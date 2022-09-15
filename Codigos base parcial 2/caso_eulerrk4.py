  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:51:05 2021

@author: usuario
"""

#librerias
import numpy as np
import matplotlib.pyplot as plt

from IPython import get_ipython
ipython = get_ipython()


#funciones

#EDO a resolver
def fode(x,y):
    dydx=x/y
    return dydx



#funcion con el metodo Euler

def meuler(xi,yi,h,f):
    x=xi+h
    y=yi+h*f(xi,yi)
    
    return x,y

def rk4(xi,Yi,h,f):
    k1=h*f(xi,Yi)
    k2=h*f(xi+0.5*h,Yi+0.5*k1)
    k3=h*f(xi+0.5*h,Yi+0.5*k2)
    k4=h*f(xi+h,Yi+k3)
    Y=Yi+(1./6.)*k1+(1./3.)*k2+(1./3.)*k3+(1./6.)*k4
    x=xi+h
    return x,Y


#######################################3
#main

#parametros de entrada

xi=1.0
yi=2.0
xf=3.0
n=5 #numero de datos
h=(xf-xi)/(n-1) #paso


#definiendo el tamaño de los vectores
x=np.zeros(n)
ynum=np.zeros(n)

ynum1=np.zeros(n)


#valores iniciales
x[0]=xi
ynum[0]=yi

#opcion 2
for i in range(1,n):
    x[i],ynum[i]=meuler(x[i-1],ynum[i-1],h,fode)
    

#valores iniciales
x[0]=xi
ynum1[0]=yi

#solucion runge kutta
for i in range(1,n):
    x[i],ynum1[i]=rk4(x[i-1],ynum1[i-1],h,fode)
    
    
#solucion analitica
yanalitica=np.sqrt(x**2+3)  #ecuacion vectorial

#error
error=np.zeros(n) #definimos tamaño vector error
error=np.abs(yanalitica-ynum)


#grafico
plt.figure(1)
plt.plot(x,ynum,'*-',color='r',label='Euler')
plt.plot(x,ynum1,'*-',color='c',label='RK4')
plt.plot(x,yanalitica,'o',color='b',label='analitica')
plt.xlabel('x')
plt.ylabel('y')
plt.title('paso h:'+str(h))
plt.legend(loc='best')

#reporte resultados
# ipython.magic('%clear') #para borrar consola

#opcion 2
solution=np.zeros(shape=(n,4))

solution[:,0]=x
solution[:,1]=ynum
solution[:,2]=yanalitica
solution[:,3]=error

print(solution)







