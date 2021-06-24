

import csv
import numpy as np 
import time
import sys 

def norm(k) :
    dp = np.dot(k, k)
    n = np.sqrt(dp)

    return n

def GS(A) : #01   

    m, n = A.shape
    G = np.zeros((m, n))
    G[:, 0] = A[:, 0]/(norm(A[:, 0]))
    for k in range(1,n) :
        for i in range(k) : 
            G[:, k] = G[:, k] + (np.dot(A[:, k], G[:, i])*G[:, i])/(np.dot(G[:, i], G[:, i]))

        G[:, k] = (A[:, k] - G[:, k])/(norm(A[:, k] - G[:, k]))

    return G

def MGS1(A) :
    m, n = A.shape
    G = np.zeros((m, n))
    G[:, 0] = A[:, 0]/(norm(A[:, 0]))
    
    for k in range(1,n) :
        G[:, k] = A[:, k].copy()
        for i in range(k) :
            G[:, k] = G[:, k] - (np.dot(G[:, k], G[:, i])*G[:, i])/(np.dot(G[:, i], G[:, i]))
            G[:, k] = G[:, k]/(norm(G[:, k]))
        
    return G

def MGS2(A) :
    m, n = A.shape
    G = np.zeros((m, n))
    G[:, 0] = A[:, 0]/(norm(A[:, 0]))

    for k in range(1,n) :
        G[:, k] = A[:, k].copy()
        for i in range(k):
            G[:, k] = G[:, k] - (np.dot(G[:, k], G[:, i])/np.dot(G[:, i], G[:, i]))*G[:, i]
            G[:, k] = G[:, k]/(norm(G[:, k]))
            
    return G

def MGS3(A) :
    m, n = A.shape
    G = np.zeros((m, n))

    for i in range(n):
        G[:, i] = A[:, i].copy()
        for j in range(n-1) :
            alpha = np.dot(G[:, j], G[:, j])
            for k in range(j+1, n) :
                G[:, k] = G[:, k] - (np.dot(G[:, k], G[:, j]))/(alpha)*G[:, j]
    
    for col in range(n):
        G[:, col] = (G[:, col])/(norm(G[:, col]))

    return G

def sigma(D) :
    m, n =np.shape(D)
    s = 0

    for i in range(n) :
        for j in range(n) :
            if j != i :
                s = s + D[i, j]**2
    s = np.sqrt(s)

    return s

def prep(adr) : 

    data_raw = list()
    cwb_filename = adr

    with open(cwb_filename, newline='') as csvfile:
        mycsv = csv.DictReader(csvfile)
        header = mycsv.fieldnames
        for row in mycsv :
            data_raw.append(row[header[0]])

    for i in range(len(data_raw)) :
        data_raw[i] = data_raw[i].split()

    data = np.array(data_raw).astype(np.int)

    return data

def main() :


    address = r'E:\VS_Code_Stations\Nunerical_Analysis\m3.dat'
    #address = str(input())
    mtx = prep(address)

    t0 = time.process_time_ns()
    g0 = MGS1(mtx)
    t1 = time.process_time_ns() - t0
    print("Time for GS = ", t1)
    delta0 = np.transpose(g0) @ g0
    print(sigma(delta0))

    t0 = time.process_time_ns()
    g1 = MGS1(mtx)
    t1 = time.process_time_ns() - t0
    print("Time for MGS1 = ", t1)
    delta1 = np.transpose(g1) @ g1
    print(sigma(delta1))

    t0 = time.process_time()
    g2 = MGS2(mtx)
    t1 = time.process_time() - t0
    print("Time for MGS2 = ", t1)
    delta2 = np.transpose(g2) @ g2
    print(sigma(delta2))
    
    '''
    t0 = time.process_time()
    g3 = MGS3(mtx)
    t1 = time.process_time() - t0
    print("Time for MGS3 = ", t1)
    delta3 = np.transpose(g3) @ g3
    print(sigma(delta3))
    '''
main()