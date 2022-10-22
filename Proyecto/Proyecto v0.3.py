# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:51:11 2022

@author: Roberto - Julian
"""

from math import erfc
import matplotlib.pyplot as plt
import numpy as np

N=10
Nt=50
Nnano=8e6
R=4e-3
K=1
pnir=1
o=3.8e-14
ua=2.2
usp=1220
ueff=np.sqrt(3*ua*(ua+usp))
Doptic=1/(3*(usp+ua))



def Tt(r,ti):
    f1 = (pnir*np.exp(-ueff*r))/(4*np.pi*Doptic*r)
    f2 = 2*(erfc(r/np.sqrt(4*K*ti))*(o*f1/(4*np.pi*K*r)))
    return Nnano*f2
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
                    Z[i1,i2,i3]=T1(rx[i1],ry[i2],0,tiempo[i3])
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
    # plt.clim([0, 14])
    plt.show()