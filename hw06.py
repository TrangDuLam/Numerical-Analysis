# EE4070 Numerical Analysis
# HW05 Matrix Condition Number
# ID : 106061121, Name : Yulan Chuang
# Date : Apr. 13, 2021


import numpy as np    
import time   
from ee4070 import *
from RV import *
import sys, os   

#error 1 : np.abs(v_tmp - v)
#error 2 : norm2(n, (q_tmp-q))
#error 3 : norm2(n, r)
#error 4 : norm2(n, r)/(q @ A)

def EVpwr(n, A, q0, iterMax = 10**6, tol = 10**(-9)) :
    
    q = q0
    v = 0
    iter = 0
    
    for iter in range(iterMax) :
        
        q_tmp = A @ q
        q_tmp = q_tmp/norm2(n, q_tmp)
        v_tmp = q_tmp @ A @ q_tmp

        r = A @ q_tmp - v_tmp * q_tmp

        if ( norm2(n, r) < tol) :
            return (True, iter+1, v)
        else :
            v = v_tmp
            q = q_tmp
            pass

    return (False, iter+1, v)
    

def EVipwr(n, A, q0, iterMax = 10**6, tol = 10**(-9)) :

    q = q0
    m = 0
    iter = 0

    for iter in range(iterMax) :

        z = linSol(n, A, q) #combination of fwdSub and bwdSub
        q_tmp = z/norm2(n, z) #normalization
        m_tmp = q_tmp @ A @ q_tmp #eigenvalue renew

        r = A @ q_tmp - m_tmp * q_tmp

        if ( norm2(n, r) < tol) :
            return (True, iter+1, m)
        else :
            m = m_tmp
            q = q_tmp
            pass

    
    
    return (False, iter+1, m)

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
    
    VsrcSym(md, A, b, 0, 0.003)  #voltage source
    VsrcSym(md, A, b, gnd, 0) #ground node

    A[0, 0] = 0.003  #modification
    A[gnd, gnd] = 0.003  #modification

    return A, b


def main() :

    n = int(input("Number of resistors per side : ")) #input n

    dim = n + 1 #nodes amount each side
    md = dim**2 #total node amount

    print("Number of variables: ", md)

    ci = circuit_system(n)[0]
    q0 = np.ones((md, ), dtype = float)

    t0 = time.process_time()
    lambda1, iter1 = EVpwr(md, ci, q0)[2], EVpwr(md, ci, q0)[1]
    t = time.process_time() - t0
    print("Lambda-1: ", lambda1)
    print("iter-1: ", iter1)
    print("CPU time: ", t, " s")

    t0 = time.process_time()
    lambda2, iter2 = EVipwr(md, ci, q0)[2], EVipwr(md, ci, q0)[1]
    t = time.process_time()-t0
    print("Lambda-n: ", lambda2)
    print("iter-n: ", iter2)
    print("CPU time: ", t, " s")

    print("Matrix condition number: ", lambda1/lambda2)

main()