# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:48:50 2022

@author: Roberto Gomez y Camilo Angarita
"""

import numpy as np
from numpy import linalg as la

def F(x):
    n=np.size(x)
    y=np.zeros(n,float)
    p=x[0]
    s=x[1]
    y[0]=-(3.182e9/((p**2)*np.sqrt(s)))+810
    y[1]=-(1.591e9/(p*s*np.sqrt(s)))+4750/np.sqrt(s)
    return y

def jacobiano_newton(x):
    n=np.size(x)
    j=np.zeros([n,n],float)
    p=x[0]
    s=x[1]
    d1dp=6364000000/(p**3*np.sqrt(s))
    d1ds=1591000000/(p**2*s**(3/2))
    d2dp=1591000000/(p**2*s*np.sqrt(s))
    d2ds=2386500000/(p*s**(5/2))-2375/s**(3/2)
    j=np.array([[d1dp,d2dp],[d1ds,d2ds]])
    return j

def costo(x):
    p=x[0]
    s=x[1]
    c=3.182e9/(p*np.sqrt(s))+9500*np.sqrt(s)+810*p
    return c

#valores iniciales
x0=np.array([1,1])
tol=1e-10
er=1
maxiter=100
i=0



print('Metodo de Newton Multivariable')
print('+++++++++++++++++++++++++++++++')

while (er>tol)and(i<=maxiter):
    j=jacobiano_newton(x0)
    f=F(x0)
    dx=np.linalg.solve(j,-f)
    xnew=x0+dx
    er=la.norm(xnew-x0,2)
    x0=np.copy(xnew)
    i+=1
    if i==maxiter:
        print('diverge')
        break

costo_total=costo(xnew)


print(f'Numero de iteraciones: {i}')
print('+++++++++++++++++++++++++++++++')
print('Solucion:')
print('+++++++++++++++++++++++++++++++')
print(f'La potencia óptima del mezclador es: {xnew[0]} KW')
print(f'La capacidad óptima del mezclador es: {xnew[1]} kg/lote')
print(f'El costo obtenido con las variables optimizadas es: {costo_total} $/año')
