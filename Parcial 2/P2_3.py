# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:19:57 2022

@authors:   Roberto Carlos Gomez Araque_ID: 000423542
            Juan Camilo Fernandez Angarita_ID: 000424705
"""
# importamos librerias
import numpy as np 

#matrices 
J=np.array([[98.,9.,2.,1.,0.5],[11.,118.,9.,4.,0.88],[27.,27.,85.,8.,2.],[1.,3.,17.,142.,25.],[2.,4.,7.,17.,118.]])
F=np.array([[0.1100],[0.2235],[0.2800],[0.3000],[0.1400]])

np.diag(J)
D=np.diag(np.diag(J))
L=np.tril(J,k=-1)
U=np.triu(J,k=1)

#desarrollo

iter=0
tolerancia=1e-10
error=1.0

x1=np.array([[0.],[0.],[0.],[0.],[0.]])

n=len(F) #con este vamos a calcular el tamano del vector
x2=np.zeros(n)

while error>=tolerancia:
    cj=np.matmul(-np.linalg.inv(L+D),(U))
    
    dj=np.matmul(np.linalg.inv(D+L),F)
    
    x2=np.matmul(cj,x1)+dj
    
    error=np.linalg.norm(x2-x1,2)  #este sera el error que tenemos con la euclidea
    iter+=1
    
    x1=x2
    

print(f'Con {iter} interaciones, la solución del sistema de ecuaciones es: ')
print(x2)

