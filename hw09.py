#!/usr/bin/env python3
# EE4070 Numerical Analysis
# HW09 Spline Interpolation
# ID : 106061121, Name : Yulan Chuang
# Date : May 11, 2021


import numpy as np  
from matplotlib import pyplot as plt 
import bisect


def splineM(N, X, Y) :
    
    h = np.diff(X)


    mu = [h[i] / (h[i] + h[i + 1]) for i in range(N - 2)] + [0] #using the attribute of python list to generate mu
    mu = np.array(mu, dtype = float)  #reform as numpy array
    twos = [2] * N #diagonal are 2s
    twos = np.array(twos, dtype = float)
    lam = [0] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(N - 2)] #using the attribute of python list to generate lamda
    lam = np.array(lam, dtype = float)  #reform as numpy array

    d = [0] + [6 * ((Y[i + 1] - Y[i]) / h[i] - (Y[i] - Y[i - 1]) / h[i - 1]) / (h[i] + h[i-1]) for i in range(1, N - 1)] + [0] ##using the attribute of python list to generate d
    
    lam_p = np.append(lam, float(0)) #append the boundary condition of lamda_n
    d_p = np.zeros((N, ), dtype= float) #declare the temporary vector d_p
    M = np.zeros((N, ), dtype= float) #declare the moment vector



    #Thomas algorithm process
    lam_p[0] = lam[0] / twos[0]
    d_p[0] = d[0] / twos[0]

    for i in range(1, N) :
        lam_p[i] = lam_p[i] / (twos[i] - lam_p[i - 1] * mu[i - 1])
        d_p[i] = (d[i] - d_p[i - 1] * mu[i - 1]) / (twos[i] - lam_p[i - 1] * mu[i - 1])

    M[-1] = d_p[-1]
    for i in range(N - 2, -1, -1):
        M[i] = d_p[i] - lam_p[i] * M[i + 1]

    return M

def spline(N, X, Y, M, x) :

    h = np.diff(X) #take the h

    coefficients = [[(M[i+1]-M[i])*h[i]*h[i]/6, M[i]*h[i]*h[i]/2, (Y[i+1] - Y[i] - (M[i+1]+2*M[i])*h[i]*h[i]/6), Y[i]] for i in range(N-1)] #using the attribute of python list to generate interpolation coefficients

    def polt_spline(val) : #define a function specified for the plotting process
        idx = min(bisect.bisect(X, val)-1, N-2) #search where the x value should be classified upon built-in binary search function  
        z = (val - X[idx]) / h[idx] #assign a new variable for convenience
        C = coefficients[idx] #return the interpolation corresponding to the given index  
        return (((C[0] * z) + C[1]) * z + C[2]) * z + C[3] #return the interpolation value

    return polt_spline #return the plot function as final output

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

    N = len(xi)

    mmt = splineM(N, xi, yi) #attain the moment vector
    result = spline(N, xi, yi, mmt, x) #generating the interpolation function

    yspline = [result(v) for v in x] #computing the interpolation value upon the function attained by spline function
    yspline = np.array(yspline) #declare the data type
    error = np.abs(y - yspline) #calculate the error 

    print("max error = ", max(error)) #return the max error 

    plt.plot(x, y, label = 'Original waveform') #plot the original waveform
    plt.scatter(xi, yi, c = 'r') #plot the supporting points
    plt.plot(x, yspline, label = 'Interpolation function') #plot the interpolation function 
    plt.plot(x, error, '--', label = 'error') #plot the error on the graph
    plt.legend() 
    plt.savefig("interp.png")

main()