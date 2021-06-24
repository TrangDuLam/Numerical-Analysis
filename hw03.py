#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW03 Resistor Network
# ID : 106061121, Name : Yulan Chuang
# Date : Apr. 25, 2021

import numpy as np
from ee4070 import *  
from RV import *
import time 
import sys, os 

def circuit_system(n) :

    g = n/2000 #conductance

    k = n + 1  #node each side
    md = k**2 # shape of the system matrix
    gnd = int(n * k + n/2) #ground point denotation

    A = np.zeros((md, md), dtype = float)
    b = np.zeros((md,), dtype = float)

    for i in range(0, md) :  #for each node
        if (i+1) % k != 0 :  #horizontal direction
            Radd(A, i, i+1, g)
        if i < k * (k-1) :  #vertical direction
            Radd(A, i, i+k, g)
    
    VsrcSym(md, A, b, 0, 1)  #voltage source
    VsrcSym(md, A, b, gnd, 0)  #ground node

    return A, b


def main() :

    n = int(input("Number of resistors per side : "))

    print("N = ", n)

    g = n/2000

    dim = n + 1 #nodes amount each side
    md = dim**2 #total node amount

    ne = n # northeast node index
    w = int(dim * n / 2) #west node index
    e = w + dim - 1 #east node index

    A, b = circuit_system(n)
    
    t0 = time.process_time()  #record the execution time
    x = linSol(md, A, b)
    t = (time.process_time()-t0) #taking the average ti

    current = g*(2 - x[1] - x[n+1])

    print("Vw = ", x[w], "V")  #print the west node
    print("Vne = ",x[ne], "V") #print the northeast node
    print("Ve = ",x[e],"V") #print the east node
    print("Req = ", 1/current,"Ohms") #print the equivalent resistance
    print("CPU time =", t, "seconds") #print the execution time




main()