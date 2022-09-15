#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:00:15 2021

@author: macbookpro
"""

#programa del pendulo simple
#resuelto con rk4
# modelo con ecuacion diferencial de orden superior

import numpy as np
import matplotlib.pyplot as plt
from math import pi


#funcion  a resolver (vectorial)
def fode(t,y):
    
    #dydt=np.array([y[1],-gra/(L*m)*np.sin(y[0])])  #sin rozamiento
    dydt=np.array([y[1],-gra/(L*m)*np.sin(y[0])-(a/m)*y[1]])  #con rozamiento
    
    return dydt


def mRK4(t,y,h,f):
    
    t=t+h
    
    k1=h*f(t,y)
    k2=h*f(t+h/2,y+k1/2)
    k3=h*f(t+h/2,y+k2/2)
    k4=h*f(t+h,y+k3)
    y=y+(k1+2*k2+2*k3+k4)/6
    
    return t,y


########### main  #################3

#parametros fisicos
L=0.6 # longitud pendulo en m
gra=9.8 #m/s2
m=1. #masa en kg
a=0.5 #coeficiente de arrastre


#parametros del tiempo
ti=0.  #tiempo inicial en s
tf=20.  #tiempo final en s
n=100 #numero de datos
h=(tf-ti)/(n-1)  #paso del tiempo

#declaracion vectores
tiempo=np.zeros(n)
y=np.zeros([n,2],float)  #matriz


#condiciones iniciales
initial_theta=np.radians(30) #posicion angular inicial en rad

initial_thetap=np.radians(0)  # velocidad angular inicial en rad/s

y[0,0]=initial_theta

y[0,1]=initial_thetap

#calculo por RK4
for i in range(1,n):
    
    tiempo[i],y[i,0:2]=mRK4(tiempo[i-1],y[i-1,0:2],h,fode)



    
#grafico
plt.figure(1)

plt.subplot(2,1,1)
plt.plot(tiempo,y[:,0]*180/pi,'o-',color='r',label='posicion angular')
plt.legend()
plt.grid()

plt.subplot(2,1,2)
plt.plot(tiempo,y[:,1]*180/pi,'--',color='b',label='velocidad angular')
plt.legend()
plt.grid()


# #animacion
# dim=len(tiempo)
# plt.figure(2)
# for i in range(0,dim):
#     plt.xlim(-1,1)
#     plt.ylim(-1,1)
    
#     posx=L*np.cos(pi/2-y[i,0])
#     posy=-L*np.sin(pi/2-y[i,0])
    
#     datax,datay=[0,posx],[0,posy]
#     plt.plot(datax,datay,color='r')
    
#     titulo='theta:'+str(y[i,0]*180/pi)
#     plt.title(titulo)
    
#     plt.grid()
#     plt.pause(0.001)#detener el código unos segundos









    
    