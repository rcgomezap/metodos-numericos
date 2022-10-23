# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:51:11 2022

@author: Julian
"""

import sympy as sy
import matplotlib.pyplot as plt
import numpy as np

k=1
P=1
U=1
D=1
o=3.8e-14
x,y,t = sy.symbols('x y t')

def T(r):
    
    f1 = P*np.exp(-U*r)/(4*np.pi*D*r)
    f = o*f1/(4*np.pi*k*r)

    return f

def T1(x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    Tx=T(r)
    return Tx

def Tt(r,ti):
    f1 = P*sy.exp(-U*x)/(4*sy.pi*D*x)
    f2 = 2*(sy.erfc(x/sy.sqrt(4*k*t))*(o*f1/(4*sy.pi*k*x)))
    f3 = sy.lambdify([x,t],f2, "numpy")
    f=f3(r,ti)
    return f

rx=np.linspace(-30e-9,30e-9,100)
ry=rx
rz=rx
tiempo=np.linspace(1e-3,0.01,100)
dim=(100,100)
Z=np.zeros(dim)
y=np.zeros(100)
Ttiempo=np.zeros(100)

plt.figure(1)
X, Y = np.meshgrid(rx, ry)
Z=T1(X,Y,5e-9)
plt.pcolor(X, Y, Z)
plt.colorbar()
plt.xlabel('Eje x (m)')
plt.ylabel('Eje y (m)')
plt.show()

for i in range (0,100):
    Ttiempo[i]=Tt(1e-9,tiempo[i])
    
plt.figure(2)
plt.plot(tiempo,Ttiempo)
plt.xlabel('Tiempo (s)')
plt.ylabel('Cambio en la temperatura (°C)')