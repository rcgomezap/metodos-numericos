# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:04:34 2022

@author: rober
"""

n=8
N=0

def pos(i,j,k):
    pos=i+n*j+n**2*k
    return pos

def pos1(N):
    k=N//n**2
    j=N%n**2//n
    i=N%n**2%n
    return i,j,k

for k in range(n):
    for j in range(n):
        for i in range(n):
            if pos(i,j,k) == pos(pos1(N)[0],pos1(N)[1],pos1(N)[2]):
                print('si')
                x=1
            else:
                print('no')
                y=0
            N+=1