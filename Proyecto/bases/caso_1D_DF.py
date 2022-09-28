#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: valencia
"""
from IPython import get_ipython
ipython = get_ipython()
ipython.magic('%clear')

import numpy as np
import matplotlib.pyplot as plt




# Define the dimensions of the domain
L = 1

# Define the number of divisions (11 for a cell size of dx = 0.1)
n = 11

# Define the cell size, set to 0.1.
dx = L/(n-1)  

# Physical parameters:
V = 1

# Boundary conditions Dirichlet
u0 = 0
u10 = 1

# Initialize the matrix for the solution
A=np.zeros([n-2,n-2],float)

B = np.zeros(n-2)

# Define the elements of diagonals in M

Dm=2./(dx**2)
Dup=-1./(dx**2) + V/(2.*dx)
Dlo=- 1./(dx**2)-V/(2.*dx)

# Assemble Matrix M Nodes 2 to 8
for i in range(1, n-3):
    A[i,i] = Dm
    A[i,i-1] =  Dlo
    A[i,i+1] = Dup

#Node 1
A[0,0] = Dm
A[0,1] = Dup

#Node 9
A[n-3, n-3] = Dm
A[n-3, n-4] = Dlo

# Assemble the vector B
B[0] = -Dlo*u0 #node 1
B[n-3] = -Dup*u10 #node 9

#solution
solution = np.linalg.solve(A, B) #numerical solution

# exact solution
x = np.linspace(0,L,n)
uanalitica = (np.exp(x) - 1)/(np.e -1)

# Plot solution
plt.figure
unum=np.hstack([u0, solution, u10])#assemble the vector u
plt.plot(x, unum, label = 'FD')
plt.plot(x, uanalitica, 'o', label = 'Analitica')
plt.title('1-D FINITE DIFFERENCE PROBLEM WITH DIRECT SOLUTION')
plt.xlabel('x')
plt.ylabel('u')
plt.legend(loc='upper left')
plt.show()