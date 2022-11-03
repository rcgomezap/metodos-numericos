# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 08:33:18 2022

@author: Roberto Gomez, Juan Fernandez
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def parabola(a0,a1,a2,x):
    y=a0+a1*x+a2*x**2
    return y

data=pd.read_csv('datos.csv', sep=',')
xdata=data['x'].to_numpy()
ydata=data['y'].to_numpy()

x=np.linspace(xdata[0],xdata[-1],100)

n=len(xdata)

sumax=sum(xdata)
sumax2=sum(xdata**2)
sumax3=sum(xdata**3)
sumax4=sum(xdata**4)
sumay=sum(ydata)
sumaxy=sum(xdata*ydata)
sumax2y=sum(xdata**2*ydata)

A=np.array([[n,sumax,sumax2],[sumax,sumax2,sumax3],[sumax2,sumax3,sumax4]])
B=np.array([[sumay],[sumaxy],[sumax2y]])

[[a0],[a1],[a2]]=np.linalg.solve(A,B)

#Validación 2: curvefit
popt,pcov=curve_fit(parabola,xdata,ydata)
a,b,c=popt




plt.figure(1)
plt.scatter(xdata,ydata,label='Datos experimentales', color='r')
plt.plot(x,parabola(a0,a1,a2,x),label='Ajuste con minimos cuadrados', color='g')
plt.legend()
plt.xlabel('Distancia (m)')
plt.ylabel('ALtura (m)')
plt.title('Ajuste con minimos cuadrados')
plt.grid()
plt.show()

plt.figure(2)
plt.plot(x,parabola(a0,a1,a2,x),label='Ajuste con minimos cuadrados', color='g')
plt.plot(x,parabola(a,b,c,x),label='Ajuste con Curvefit', color='r')