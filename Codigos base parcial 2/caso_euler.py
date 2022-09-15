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


#valores iniciales
x[0]=xi
ynum[0]=yi


for i in range(1,n):
    x[i],ynum[i]=meuler(x[i-1],ynum[i-1],h,fode)
    


#solucion analitica
yanalitica=np.sqrt(x**2+3)  #ecuacion vectorial

#error
error=np.zeros(n) #definimos tamaño vector error
error=np.abs(yanalitica-ynum)


#grafico
plt.figure(1)
plt.plot(x,ynum,'*-',color='r',label='Euler')
plt.plot(x,yanalitica,'o-',color='b',label='analitica')
plt.xlabel('x')
plt.ylabel('y')
plt.title('paso h:'+str(h))
plt.legend(loc='best')


#opcion 2
solution=np.zeros(shape=(n,4))

solution[:,0]=x
solution[:,1]=ynum
solution[:,2]=yanalitica
solution[:,3]=error

print(solution)







