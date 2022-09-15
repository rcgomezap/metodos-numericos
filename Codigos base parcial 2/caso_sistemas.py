 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 17:25:53 2022

@author: usuario
"""
from IPython import get_ipython
ipython = get_ipython()

import numpy as np
import matplotlib.pyplot as plt 


from IPython import get_ipython
ipython = get_ipython()

def fode(x,y):

    dydx=np.array([-0.5*y[0],4-0.3*y[1]-0.1*y[0]])
    return dydx


def meuler(xi,yi,h,f):
    x=xi+h
    y=yi+h*f(xi,yi)
    
    return x,y

####################
#MAIN

#Parametros iniciales
xi=0.0
xf=2.0
n=5 #numero de datos

h=(xf-xi)/(n-1) #paso

#condicion inicial
yi=np.array([4,6])


solution=np.zeros((n,3))


solution[0,0]=xi

solution[0,1]=yi[0]

solution[0,2]=yi[1]



for i in range(1,n):
    x,y=meuler(xi,yi,h,fode)
    solution[i,:]=[x,y[0],y[1]]
    xi=x
    yi=y
    


np.savetxt("resultado.txt",solution)
    
#reporte resultados
ipython.magic('%clear') #para borrar consola

print("Solución de sistema ODE por Euler")
print("+++++++++++++++++++++++++++++++++++")
print(solution)

np.savetxt('resultadoseuler.txt',solution)


#grafico
plt.figure(1)
plt.plot(solution[:,0],solution[:,1],'*-',color='r',label='y1')
plt.plot(solution[:,0],solution[:,2],'o-',color='b',label='y2')
plt.xlabel('x')
plt.ylabel('y')
plt.title('paso h:'+str(h))
plt.legend(loc='best')