# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 20:19:04 2022

@author: rober
"""
#librerias
import numpy as np
import matplotlib.pyplot as plt
import sympy as sy

T=5
R=8.314
Ea=20
A=2

t = sy.Symbol('t')
resultado = float(sy.integrate(A*sy.exp(-Ea/R*T),(t,1,4)).evalf())
print(resultado)