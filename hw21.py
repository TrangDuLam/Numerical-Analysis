#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW21 Linear Solution Methods
# ID : 106061121, Name : Yulan Chuang
# Date : Jun. 19, 2021

import numpy as np
from scipy import linalg
from ee4070 import *
import sys, os, time


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
'''

    scipy no heading // i.e. not to type scipy.linalg

    test case 1 :  lu 
    code :

       //SciPy
       lu_fac = linalg.lu_factor(A)
       y = linalg.lu_solve(lu_fac, b) 
       //ee4070
       y = linSol(matrix dimension, A, b)


    test case 2 : direct solve 
        //SciPy
        y = np.linalg.solve(A, b)
        y = linalg.solve(A, b)
    
    test case 3 : inverse matrix
        //NumPy 
        A_inv = np.linalg.inv(A)
        y = np.matmul(A_inv, b)
        //SciPy 
        A_inv = linalg.inv(A)
        y = np.matmul(A_inv, b)

    test case 4 : Cholesky 
        //SciPy 
        cho_fac = linalg.cho_factor(A)
        y = linalg.cho_solve(cho_fac, b)

'''


def main() :


    side = np.array([10, 20, 40, 50, 60, 80, 100, 120], dtype = int)

    eff = [] #record the time complexity

    #direct solve

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        y = np.linalg.solve(A, b)
        t = time.process_time() - t0

        eff.append(t)

    print("Time for NumPy direct solve", eff)
    eff.clear()

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        y = linalg.solve(A, b)
        t = time.process_time() - t0

        eff.append(t)

    print("Time for SciPy direct solve", eff)
    eff.clear()

    #LU 

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        lu_fac = linalg.lu_factor(A)
        y = linalg.lu_solve(lu_fac, b) 
        t = time.process_time() - t0

        eff.append(t)

    print("Time for SciPy LU solve", eff)
    eff.clear()

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        y = linSol((i+1)**2, A, b)
        t = time.process_time() - t0

        eff.append(t)

    print("Time for ee4070 LU solve", eff)
    eff.clear()

    #inverse matrix

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        A_inv = np.linalg.inv(A)
        y = np.matmul(A_inv, b)
        t = time.process_time() - t0

        eff.append(t)

    print("Time for NumPy inverse solve", eff)
    eff.clear()

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        A_inv = linalg.inv(A)
        y = np.matmul(A_inv, b)
        t = time.process_time() - t0

        eff.append(t)

    print("Time for SciPy inverse solve", eff)
    eff.clear()

    #cholesky

    for i in side :

        A, b = circuit_system(i)
        
        t0 = time.process_time()
        cho_fac = linalg.cho_factor(A)
        y = linalg.cho_solve(cho_fac, b)
        t = time.process_time() - t0

        eff.append(t)

    print("Time for Cholesky solve", eff)
    eff.clear()


main()
