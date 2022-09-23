# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 08:26:45 2021

@author: USUARIO
"""

import numpy as np
import matplotlib.pyplot as plt
#-------------------------------- CONSTANTES ----------------------------------

rho=1000 #densidad del agua [kg/m3]

Cp=3400 #capacidad calorífica [J/kg-K]

Q=0.008e3 #Energia del laser (W/m3)
 
OD=0.2

#------------------------------- PARÁMETROS -----------------------------------

L=0.1 #Grosor de la pared
n=10
T0=0
T1s=37+273.25
T2s=45+273.15
dx=L/n
alpha=3e-7
t_final=60+273.15
dt=0.1

x=np.linspace(dx/2,L-dx/2,n)

T=np.ones(n)*T0 #Crea un vector de uno para luego multiplicarlo por la T0
dTdt=np.empty(n) #Se define un vector vacio

t=np.arange(0,t_final,dt) #Se crea un vector con el tiempo

for j in range(1,len(t)):
    
    for i in range(1,n-1):
        dTdt[i] = alpha*((T[i+1]-2*T[i]+T[i-1])/dx**2) + Q/(rho*Cp)*(1-10**(-OD))
    dTdt[0]= alpha*((T[1]-2*T[0]+T1s)/dx**2) + Q/(rho*Cp)*(1-10**(-OD))
    dTdt[n-1]=alpha*((T2s-2*T[n-1]+T[n-2])/dx**2) + Q/(rho*Cp)*(1-10**(-OD))
    T= T+ dTdt*dt
    plt.figure(1)
    plt.plot(x,T)
    plt.axis([0,L,0,50])
    plt.xlabel('Distancia (m)')
    plt.ylabel('Temperatura(°C)')
    plt.show
    plt.pause(1e-20)