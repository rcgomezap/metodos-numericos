# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#pendulo simple - anim
#resuelto con RK4
#modelo de ecuacion diferencial de orden superior

#lib
import numpy as np
import matplotlib.pyplot as plt
from math import pi

#funcion a resolver
def teta(t,y):
    dydt=np.array([y[1],-g/(m*L)*np.sin(y[0])])
    return dydt

def rk4(t,y,h,f):
    t=t+h
    
    k1=h*f(t,y)
    k2=h*f(t*h/2,y+k1/2)
    k3=h*f(t*h/2,y+k2/2)
    k4=h*f(t+h,y+k3)
    
    y=y+(k1+2*k2+2*k3+k4)/6
    return t,y

#main

#parametros fisicos
L=0.6 #longitud (m)
g=9.8 #gravedad (m/s2)
m=1.0 #masa (kg)
a=0.5 #coef de arrastre


#parametros de tiempo
ti=0.0 #tiempo inicial (s)
tf=20.0 #tiempo final (s)
n=100 #numero de datos
h=(tf-ti)/n-1 #paso
#declarar vectores
tiempo=np.zeros(n)
y=np.zeros([n,2],float) #matriz solucion
solution=np.zeros([n,3],float)


#condiciones iniciales
teta_inicial=np.radians(45) #posicion angular inicial
tetap_inicial=np.radians(0) #velocidad angular inicial

y[0,0]=teta_inicial
y[0,1]=tetap_inicial


#ciclo
i=1
while i<n:
    tiempo[i],y[i,:]=rk4(tiempo[i-1],y[i-1,:], h, teta)
    i=i+1
    
solution[:,0]=tiempo
solution[:,1]=y[:,0]
solution[:,2]=y[:,1]

print(solution)

    