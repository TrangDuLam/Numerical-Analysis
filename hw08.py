#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW08 Lagrange Interpolation
# ID : 106061121, Name : Yulan Chuang
# Date : May 4, 2021

import numpy as np 
from matplotlib import pyplot as plt  


def Lagrange(x, xdata, ydata) :  #(x-axis for fn generating, interpolation x, interpolation y)
                                 # len. xdata = len. ydata

    fi = np.zeros(len(x))  #initialize interpolation function

    for j in range(len(x)) : #for every point in x-axis
    
        for i in range(len(xdata)) : #for every interpolation data point x
            L = 1  #setting L as 1

            for k in range(len(xdata)) :  
                if i != k :   #forevery x_i not equal to x_k
                    L *= (x[j] - xdata[k])/(xdata[i] - xdata[k])  #Lagrange value renew

            fi[j] += ydata[i]*L   #renewing


    return fi  #return the interpolation result

def MaxErr(y, lagi) :  #checking the maximum error in range(550, 700)

    y_sl = np.abs(y[75:226])  #550 - 475 = 75
    lagi_sl = np.abs(lagi[75:226]) #700 - 475 + 1 =226

    mxerr = np.max(np.abs(y_sl - lagi_sl))

    return mxerr

def main() :
    
    filename = "f301.dat"  #import the original waveform
    x, y = np.genfromtxt(filename, skip_header=1, unpack=True, dtype=float)  #read data from the original waveform

    L = list() #list for interpolation data
    xi = list() #x point of interpolation data
    yi = list() #y point of interpolation data

    while True : #input the interpolation data
        try :
            t = input()
            t1 = t.split()
            L.append(t1)
        except EOFError :  #checking the end of file
            break

    for i in range(1, len(L)) :  #data seperation
        xi.append(float(L[i][0]))
        yi.append(float(L[i][1]))

    xi = np.array(xi).astype(np.float) #declaring data type 
    yi = np.array(yi).astype(np.float) #declaring data type
    

    lagi = Lagrange(x, xi, yi) #do lagrange

    print("Max error of whole range = ", np.max(np.abs(y - lagi)))
    print("Max error in range 550 to 700 = ", MaxErr(y, lagi))

    
    plt.plot(x, y, label = 'Original waveform')
    plt.scatter(xi, yi, c = 'g')
    plt.plot(x, lagi, label = 'Interpolation function')
    plt.legend()
    plt.savefig("interp.png")

main()
