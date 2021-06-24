#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW05 Matrix Condition Number
# ID : 106061121, Name : Yulan Chuang
# Date : Apr. 13, 2021


import numpy as np    
import time   
from ee4070 import *
import sys, os   


def error(A) : #error calculation

    n = np.shape(A)[0]
    A = np.abs(A) #taking the absolute value

    max = A[1][0] #max by default

    for i in range(1, n) : #max checking
        if max < A[i][i-1] : 
            max = A[i][i-1]

    
    return max

def EVqr(n, A, maxIter = 10**6, tol = 10**(-9)) :

    A = A.astype(float) #declare datatype
    T = np.copy(A) #copy A to T

    iter = 0

    for iter in range(maxIter) :
        Q, R = QR(n, T) #do qr decomposition
        T_tmp = R @ Q #R times Q
        err = error(T_tmp) #take error
        ev = np.diagonal(T_tmp) #take the eigenvalues

        if(err < tol) : #find the convergence
            return  (True, iter, ev) #values return
            
        else :
            T = T_tmp # renew T

    

    return (False, iter, ev) # value return

def EVqrShifted(n, A, mu = 0.5, maxIter = 10**6, tol = 10**(-9)):

    A = A.astype(float) #declare datatype
    T = np.copy(A) #copy A to T

    iter = 0

    for iter in range(maxIter) :
        Q, R = QR(n, T - mu*np.identity(n)) #do qr decomposition
        T_tmp = R @ Q + mu*np.identity(n) #R times Q
        err = error(T_tmp) #take error
        ev = np.diagonal(T_tmp) #take the eigenvalues

        if(err < tol) : #find the convergence
            return  (True, iter, ev) #values return
            
        else :
            T = T_tmp # renew T

    
    return (False, iter, ev) # value return

def readMat(n) :

    mat = np.empty((n,n), dtype=float)
    for i in range(n) :
        t = input()
        t1 = t.split()
        for j in range(n) :
            mat[i, j] = float(t1[j])
    
    return mat

def main():
    
    
    string = input()
    n = int(string)
    A = readMat(n)
    print(n)
    t0 = time.process_time()
    print(EVqr(n, A))
    t = time.process_time() - t0
    print("CPU time : ", t, "s")

    t0 = time.process_time()
    print(EVqrShifted(n, A))
    t = time.process_time() - t0
    print("CPU time : ", t, "s")

    


main()