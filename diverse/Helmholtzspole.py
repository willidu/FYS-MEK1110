# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 20:54:07 2022

@author: Administrator
"""


import numpy as np
import matplotlib.pyplot as plt

#Variables
N= 330
R= 0.7#m
I= 1.00 #A
mu = 1.26*10**-6

numsteps = 20
x=np.linspace(-20,20,numsteps)

def helmholtzspole(x, a):
    B = (N*mu*I)/(2*R) * ((1 + ((x-a/2)**2)/(R**2))**(-3/2) + (1 + ((x+a/2)**2)/(R**2))**(-3/2))
    return B

with open('b1.txt', 'w') as file:
    result = helmholtzspole(x, 2*R)
    for e in result:
        file.write(f'{e}\n')

with open('b2.txt', 'w') as file:
    result = helmholtzspole(x, R)
    for e in result:
        file.write(f'{e}\n')

with open('b3.txt', 'w') as file:
    result = helmholtzspole(x, R/2)
    for e in result:
        file.write(f'{e}\n')

B1 = helmholtzspole(x, 2*R)
B2 = helmholtzspole(x, R)
B3 = helmholtzspole(x, R/2)

#Plot
plt.figure('helmholtzspole')
plt.plot(x, B1, label='a=2R')
plt.plot(x, B2, label='a=R')
plt.plot(x, B3, label='a=R/2')
plt.legend()
plt.show()

'--------------'

def solenoide(z, R=0.05, mu=0):
    t1 = np.arctan(R/z)
    return t1
z = np.arage(0, 100)
solenoide(z)

from spolefil import solenoide