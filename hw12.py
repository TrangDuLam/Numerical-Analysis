
# EE4070 Numerical Analysis
# HW12 RLC circuit
# ID : 106061121, Name : Yulan Chuang
# Date : Jun. 7, 2021

import numpy as np 
from matplotlib import pyplot as plt
from ee4070 import *

def fwdEuler(h) : #h : time step
    V = 1
    R = 1
    Rs = 0.1
    C = 1
    Cs = 0.1
    L = 1
    v1 = [0]
    v2 = [0]
    v3 = [0]
    iL = [0]

    t = 0 

    while t <= 10 : 

        v1_next = v1[-1] + (h/(Rs*Cs)) * (V - v1[-1]) - (h/Cs) * iL[-1]
        v1.append(v1_next)
    
        iL_next = iL[-1] - (h/L) * (v3[-1] - v2[-1])
        iL.append(iL_next)

        v2_next = v1[-1] - iL[-1] * R
        v2.append(v2_next)

        v3_next = v3[-1] + (h/C)*iL[-2]
        v3.append(v3_next)

        t += h 

    return v1, v2, v3, iL

def bwdEuler(h) :
    V = 1
    R = 1
    Rs = 0.1
    C = 1
    Cs = 0.1
    L = 1
    v1 = [0]
    v2 = [0]
    v3 = [0]
    iL = [0]

    t = 0 

    while t <= 10 : 

        cir = np.zeros((4, 4), dtype = float)
        var = np.zeros((4,), dtype = float)

        cir[0, 0] = 1/h + 1/(Cs*Rs)
        cir[0, 3] = 1/Cs
        var[0] = V/(Rs*Cs) + v1[-1]/h

        cir[1, 1] = -1/L
        cir[1, 2] = 1/L
        cir[1, 3] = 1/h
        var[1] = iL[-1]/h

        cir[2, 0] = 1
        cir[2, 1] = -1
        cir[2, 3] = -R

        cir[3, 2] = 1/h
        cir[3, 3] = -1/C
        var[3] = v3[-1]/h

        y = linSol(4, cir, var)

        v1.append(y[0])
        v2.append(y[1])
        v3.append(y[2])
        iL.append(y[3])
    
        t += h

    return v1, v2, v3, iL

def trapezoidal(h) :
    
    V = 1
    R = 1
    Rs = 0.1
    C = 1
    Cs = 0.1
    L = 1
    v1 = [0]
    v2 = [0]
    v3 = [0]
    iL = [0]

    t = 0 

    while t <= 10 : 

        cir = np.zeros((4, 4), dtype = float)
        var = np.zeros((4,), dtype = float)

        #variables assignment
        cir[0, 0] = 1 + h/(2*Cs*Rs)
        cir[0, 3] = h/(2*Cs)
        var[0] = (1 - (h/(2*Rs*Cs))) * v1[-1] + (h*V/(Rs*Cs)) - (h/(2*Cs))*iL[-1] 

        cir[1, 1] = -h/(2*L)
        cir[1, 2] = h/(2*L)
        cir[1, 3] = 1
        var[1] = iL[-1] + h/(2*L)*v2[-1] - h/(2*L)*v3[-1]

        cir[2, 0] = 1
        cir[2, 1] = -1
        cir[2, 3] = -R

        cir[3, 2] = 1
        cir[3, 3] = -h/(2*C)
        var[3] = v3[-1] + h/(2*C)*iL[-1]

        #solving
        y = linSol(4, cir, var)

        v1.append(y[0])
        v2.append(y[1])
        v3.append(y[2])
        iL.append(y[3])
    
        t += h 

    return v1, v2, v3, iL

def main() : 

    h = 0.01
    v1, v2, v3, iL = bwdEuler(h)

    print("v1 max : ", max(v1), " min : ", min(v1))
    print("v2 max : ", max(v2), " min : ", min(v2))
    print("v3 max : ", max(v3), " min : ", min(v3))
    print("iL max : ", max(iL), " min : ", min(iL))

    time_interval = int((10/h)) + 2

    time_scope = np.linspace(0, 10, time_interval)

    volt = plt.figure(1)
    plt.plot(time_scope, v1, label = "v1")
    plt.plot(time_scope, v2, label = "v2")
    plt.plot(time_scope, v3, label = "v3")
    plt.xlabel("time (sec)")
    plt.ylabel("voltage (V)")
    plt.legend()
    plt.savefig("volt_be.png")

    current = plt.figure(2)
    plt.plot(time_scope, iL, label = "iL")
    plt.xlabel("time (sec)")
    plt.ylabel("current (A)")
    plt.savefig("iL_be.png")

    


main()
