#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW05 Conjugate Gradient Method
# ID : 106061121, Name : Yulan Chuang
# Date : Apr. 25, 2021

import numpy as np    
import time   
from ee4070 import *
from RV import *
import sys, os   


def circuit_system(n) :

    g = n/2000 #conductance

    k = n + 1  #node each side
    md = k**2 # shape of the system matrix
    gnd = int(n * k + n/2) #ground point denotation

    A = np.zeros((md, md), dtype = float)
    b = np.zeros((md,), dtype = float)

    for i in range(0, md) : #for each node
        if (i+1) % k != 0 :  #horizontal direction
            Radd(A, i, i+1, g)
        if i < k * (k-1) : #vertical direction
            Radd(A, i, i+k, g)
    
    VsrcSym(md, A, b, 0, 1) #source voltage
    VsrcSym(md, A, b, gnd, 0) #ground node

    return A, b

def main() :
    n = int(input("Number of resistors per side : ")) #input n

    print("N = ", n)

    g = n / 2000

    dim = n + 1 #nodes amount each side
    md = dim**2 #total node amount

    ne = n # northeast node index
    w = int(dim * n / 2) #west node index
    e = w + dim - 1 #east node index

    ci , b = circuit_system(n)

    t0 = time.process_time()
    CG(md,ci, b)
    t = (time.process_time() - t0) #average operation time
    

    round = CG(md,ci, b)[1] #return the iteration time
    volt = CG(md, ci, b)[2] #return the solution 
    current = g*(2-volt[1]-volt[n+1]) #return the total current

    print("Node w : ", volt[w], "V")
    print("Node ne : ", volt[ne],"V")   
    print("Node e : ", volt[e],"V")
    print("Equivalent Resistance : ", 1/current,"ohm")
    print("Iteration : ", round, "times")
    print("CPU time : ", t)

main()