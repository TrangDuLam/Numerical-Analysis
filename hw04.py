#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW05 Linear Iterative Method
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
    
    VsrcSym(md, A, b, 0, 1) #source voltage
    VsrcSym(md, A, b, gnd, 0) #ground node

    return A, b


def main() :

    n = int(input("Number of resistors per side : "))

    print("N = ", n)


    dim = n + 1 #nodes amount each side
    md = dim**2 #total node amount

    ci, b = circuit_system(n)


    #Jacobi part
    t0 = time.process_time()  #record the execution time
    Jacobi(md, ci, b, enorm = norm1)
    t = (time.process_time()-t0) #taking the average time
    print("Jacobi 1 norm : ", t)

    t0 = time.process_time()  #record the execution time
    Jacobi(md, ci, b)
    t = (time.process_time()-t0) #taking the average time
    print("Jacobi 2 norm : ", t)

    t0 = time.process_time()  #record the execution time
    Jacobi(md, ci, b, enorm = normInf)
    t = (time.process_time()-t0) #taking the average time
    print("Jacobi Inf norm : ", t)


    #Gauss-Seidel part 
    t0 = time.process_time()  #record the execution time
    GS(md, ci, b, enorm = norm1)
    t = (time.process_time()-t0) #taking the average time
    print("GS 1 norm : ", t)

    t0 = time.process_time()  #record the execution time
    GS(md, ci, b)
    t = (time.process_time()-t0) #taking the average time
    print("GS 2 norm : ", t)

    t0 = time.process_time()  #record the execution time
    GS(md, ci, b, enorm = normInf)
    t = (time.process_time()-t0) #taking the average time
    print("GS Inf norm : ", t)

    #SGS part 
    t0 = time.process_time()  #record the execution time
    SGS(md, ci, b, enorm = norm1)
    t = (time.process_time()-t0) #taking the average time
    print("SGS 1 norm : ", t)

    t0 = time.process_time()  #record the execution time
    SGS(md, ci, b)
    t = (time.process_time()-t0) #taking the average time
    print("SGS 2 norm : ", t)

    t0 = time.process_time()  #record the execution time
    SGS(md, ci, b, enorm = normInf)
    t = (time.process_time()-t0) #taking the average time
    print("SGS Inf norm : ", t)




main()