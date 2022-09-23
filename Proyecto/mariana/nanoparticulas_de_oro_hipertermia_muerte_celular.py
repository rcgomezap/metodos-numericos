# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 16:25:06 2021

@author: USUARIO
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


T=np.array([37.0,60.0,100.0])


tk = 40.5
kfs = 3.33e-3
kb = 7.77e-3

#EDO a resolver
rangos=np.array([3000,2000,500])
def fode(Z, t, iteracion):
    
    A, D = Z
    kf = kfs*np.exp(T[iteracion]/tk)*(1-A)    
    
    return[-(kf)*A+kb*(1-A-D), kf*(1-A-D)]

for i in range(len(T)):
    y0 = [1, 0.01]    
    t = np.linspace(0, rangos[i], 50)
    
    sol = odeint(fode, y0, t, args=tuple([i]))

    V=1-sol[:,0]-sol[:,1]
    
    plt.figure(i)
    plt.plot(t, sol[:, 0], 'purple', label='A')
    plt.plot(t, sol[:, 1], 'black', label='D')
    plt.plot(t,V,'orange',label='V')
    plt.legend(loc='best')
    plt.title(f'T={T[i]} °C',fontweight="bold")
    plt.xlabel('Tiempo(s)',fontweight="bold")
    plt.ylabel('Fracción celular',fontweight="bold")
    plt.tight_layout()
    plt.grid(True)
    
#EDO a para la muerte celular lenta

def fode2(t):
    D=t
    kss=0.316e-3
    Ds=0.208
    ks=kss*D*(1-D)*(D-Ds)**2
    return ks
 
t = np.linspace(0, 1, 500)

plt.figure(4)
plt.plot(t, fode2(t), 'g')
plt.legend(loc='best')
plt.title('Muerte celular lenta',fontweight="bold")
plt.xlabel('Fracción celulas muertas',fontweight="bold")
plt.ylabel('Ks',fontweight="bold")
plt.tight_layout()
plt.grid(True)
plt.show()



