# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 22:39:38 2022

@author: Roberto Gomez
"""
#Librerías
import numpy as np
from scipy.optimize import fsolve

#Defino el sistema de ecuaciones 
def F(x):
    T=378
    P=1000
    F=np.zeros(2,float)
    F[0]=(50-x[0])*(x[0]+x[1])*np.exp(16.1753-2948.78/(T-44.5633))-x[0]*(100-x[0]-x[1])*P
    F[1]=(50-x[1])*(x[0]+x[1])*np.exp(16.2655-3242.38/(T-47.1806))-x[1]*(100-x[0]-x[1])*P
    return F

#Coloco mis valores iniciales en el vector x0
x0=np.array([30,20])

solution=fsolve(F,x0)

print('La solución del sistema es:')
print('V1:',solution[0])
print('V2:',solution[1])

