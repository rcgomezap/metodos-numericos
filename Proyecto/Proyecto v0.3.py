# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:51:11 2022

@author: Roberto - Julian
"""

import sympy as sy
import matplotlib.pyplot as plt
import numpy as np

N=10
Nt=50
Nnano=8e13
R=0.1
k=1
P=1
U=1
D=1
o=3.8e-14
x,y,t = sy.symbols('x y t')


def Tt(r,ti):
    f1 = P*sy.exp(-U*x)/(4*sy.pi*D*x)
    f2 = 2*(sy.erfc(x/sy.sqrt(4*k*t))*(o*f1/(4*sy.pi*k*x)))
    f3 = sy.lambdify([x,t],f2, "numpy")
    f=Nnano*f3(r,ti)
    return f

def T1(x,y,z,t):
    r=np.sqrt(x**2+y**2+z**2)
    Tx=Tt(r,t)
    return Tx

rx=np.linspace(-R,R,N)
ry=rx
rz=rx
tiempo=np.linspace(1e-5,1,Nt)
dim=(N,N,Nt)
Z=np.zeros(dim)

for i3 in range(0,Nt):
    for i2 in range (0,N):
        for i1 in range (0,N):
                    Z[i1,i2,i3]=T1(rx[i1],ry[i2],0.05,tiempo[i3])
                    print(f'{i1} {i2} {i3}')

        

plt.figure(1)
X, Y = np.meshgrid(rx, ry)
plt.pcolor(X, Y, Z[:,:,2])
plt.colorbar()
plt.xlabel('Eje x (m)')
plt.ylabel('Eje y (m)')
plt.show()

plt.plot(tiempo,Z[4,4,:])

for i in range (0,N):

    plt.figure(3)
    X, Y = np.meshgrid(rx, ry)
    plt.pcolormesh(X, Y, Z[:,:,i], cmap="plasma", shading=("gouraud"))
    plt.colorbar()
    plt.xlabel('Eje x (m)')
    plt.ylabel('Eje y (m)')
    plt.clim([0, 14])
    plt.show()