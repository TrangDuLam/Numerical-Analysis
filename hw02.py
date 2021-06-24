import numpy as np
import csv
import time
import sys, os
from ee4070 import *

def main():

    cwb_filename = str(input()) #input the filenam

    data_raw = list()

    with open(cwb_filename, newline='') as csvfile:
        mycsv = csv.DictReader(csvfile)  #read in the data file
        header = mycsv.fieldnames
        for row in mycsv :
            data_raw.append(row[header[0]]) #take the rows into another list

    for i in range(len(data_raw)) :
        data_raw[i] = data_raw[i].split() #split the element 

    #process of seperating matrix A and b
    #Since dim(A) = dim(b)
    w = len(data_raw)  
    n = int(w/2) # take the dimension of 2 matrices

    A = data_raw[:n]  #taking the matrix A
    A = np.array(A).astype(np.float) #storing the element as float type
    b = data_raw[n:] #taking the matrix b
    b = np.array(b).astype(np.float) #storing the element as float type

    t0 = time.process_time()  #record the execution time
    for i in range(100) : #do the solution process
       y = fwdSubs(n, A, b)
       bwdSubs(n, A, y)

    t = (time.process_time()-t0)/100 #taking the average time
    print(t) #print the time


main()
    

    