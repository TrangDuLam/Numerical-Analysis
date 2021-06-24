
# EE4070 Numerical Analysis
# Project Gear's Method
# ID : 106061121, Name : Yulan Chuang
# Date : Jun. 16, 2021

import numpy as np 
from matplotlib import pyplot as plt 
from ee4070 import *
import sys, os

def bwd_initial(h, order) :

    '''
    to solve the initial conditions for each variable and store them in the list
    For i-th order Gear, we need the initial conditions of t = 0 to t = (i-1)*h
    '''

    V = 1.0
    R = 1.0
    Rs = 0.1
    C = 1.0
    Cs = 0.1
    L = 1.0


    # initial list
    v1 = [0]
    v2 = [0]
    v3 = [0]
    iL = [0]
    
    t = 0

    while t <= h*(order - 2) : # time stepping

        cir = np.zeros((4, 4), dtype = float)
        var = np.zeros((4,), dtype = float)
        '''
        Using the linear equations to solve the value of the next step t + h
        All the ODE could be described into the difference equations

        '''

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

        #solving and value renewing

        v1.append(y[0])
        v2.append(y[1])
        v3.append(y[2])
        iL.append(y[3])
    
        t += h

    return v1, v2, v3, iL 

def Gear2(h) : 

    V = 1.0
    R = 1.0
    Rs = 0.1
    C = 1.0
    Cs = 0.1
    L = 1.0

    t = h #the starting time

    v1, v2, v3, iL = bwd_initial(h, 2) #initial condition 

    '''
    x(t+h) = a_1 * x(t) + a_2 * x(t-h) + a_3 * h * f(t+h)

    as t = h is the first round iteration

    '''

    a1 = float(4/3)
    a2 = float(-1/3)
    a3 = float(2*h/3)

    while t <= 20 :

        cir = np.zeros((4, 4), dtype = float)
        var = np.zeros((4,), dtype = float)

        cir[0, 0] = 1 + a3/(Cs*Rs)
        cir[0, 3] = a3/Cs
        var[0] = a1 * v1[-1] + a2 * v1[-2] + (a3 * V) / (Rs *Cs)

        cir[1, 1] = -a3/L
        cir[1, 2] = a3/L
        cir[1, 3] = 1
        var[1] = a1*iL[-1] + a2*iL[-2]

        cir[2, 0] = -1
        cir[2, 1] = 1
        cir[2, 3] = R

        cir[3, 2] = 1
        cir[3, 3] = -a3/C
        var[3] = a1*v3[-1] + a2*v3[-2]

        y = linSol(4, cir, var)

        v1.append(y[0])
        v2.append(y[1])
        v3.append(y[2])
        iL.append(y[3])

        t += h

    return v1, v2, v3, iL

def Gear3(h) :

    V = 1.0
    R = 1.0
    Rs = 0.1
    C = 1.0
    Cs = 0.1
    L = 1.0

    t = 2*h #staring time

    '''
    x(t+h) = a_1 * x(t) + a_2 * x(t-h) + a_3 * x(t-2h) + a_4*h*f(t+h)

    as t = 2h is the first round iteration

    '''

    v1, v2, v3, iL = bwd_initial(h, 3) #initial condition

    a1 = float(18/11)
    a2 = float(-9/11)
    a3 = float(2/11)
    a4 = float(6*h/11)

    while t <= 20 :

        cir = np.zeros((4, 4), dtype = float)
        var = np.zeros((4,), dtype = float)

        cir[0, 0] = 1 + a4/(Rs*Cs)
        cir[0, 3] = a4 / Cs
        var[0] = a1*v1[-1] + a2*v1[-2] + a3*v1[-3] + (a4*V)/(Rs*Cs)

        cir[1, 1] = -a4/L
        cir[1, 2] = a4/L
        cir[1, 3] = 1
        var[1] = a1*iL[-1] + a2*iL[-2] + a3*iL[-3]

        cir[2, 0] = -1
        cir[2, 1] = 1
        cir[2, 3] = R

        cir[3, 2] = 1
        cir[3, 3] = -a4/C
        var[3] = a1*v3[-1] + a2*v3[-2] + a3*v3[-3]
        
        y = linSol(4, cir, var)

        v1.append(y[0])
        v2.append(y[1])
        v3.append(y[2])
        iL.append(y[3])

        t += h

    return v1, v2, v3, iL

def Gear4(h) :

    V = 1.0
    R = 1.0
    Rs = 0.1
    C = 1.0
    Cs = 0.1
    L = 1.0

    t = 3*h

    '''
    x(t+h) = a_1 * x(t) + a_2 * x(t-h) + a_3 * x(t-2h) + a_4 * x(t-3h) + a_5 * h * f(t+h)

    as t = 3h is the first round iteration

    '''

    v1, v2, v3, iL = bwd_initial(h, 4)

    a1 = float(48/25)
    a2 = float(-36/25)
    a3 = float(16/25)
    a4 = float(-3/25)
    a5 = float(12*h/25)

    while t <= 20 :

        cir = np.zeros((4, 4), dtype = float)
        var = np.zeros((4,), dtype = float)

        cir[0, 0] = 1 + a5/(Rs*Cs)
        cir[0, 3] = a5 / Cs
        var[0] = a1*v1[-1] + a2*v1[-2] + a3*v1[-3] + a4*v1[-4] + (a5*V)/(Rs*Cs)

        cir[1, 1] = -a5/L
        cir[1, 2] = a5/L
        cir[1, 3] = 1
        var[1] = a1*iL[-1] + a2*iL[-2] + a3*iL[-3] + a4*iL[-4]

        cir[2, 0] = -1
        cir[2, 1] = 1
        cir[2, 3] = R

        cir[3, 2] = 1
        cir[3, 3] = -a5/C
        var[3] = a1*v3[-1] + a2*v3[-2] + a3*v3[-3] + a4*v3[-4]
        
        y = linSol(4, cir, var)

        v1.append(y[0])
        v2.append(y[1])
        v3.append(y[2])
        iL.append(y[3])

        t += h

    return v1, v2, v3, iL


def main() : 

    h = 0.1
    v1, v2, v3, iL = Gear4(h)

    print("v1 max : ", max(v1), " min : ", min(v1))
    print("v2 max : ", max(v2), " min : ", min(v2))
    print("v3 max : ", max(v3), " min : ", min(v3))
    print("iL max : ", max(iL), " min : ", min(iL))

    time_scope = np.linspace(0, 20, len(v1))

    volt = plt.figure(1)
    plt.plot(time_scope, v1, label = "v1")
    plt.plot(time_scope, v2, label = "v2")
    plt.plot(time_scope, v3, label = "v3")
    plt.xlabel("time (sec)")
    plt.ylabel("voltage (V)")
    plt.legend()
    plt.savefig("vlpj_order_h.png")

    current = plt.figure(2)
    plt.plot(time_scope, iL, label = "iL")
    plt.xlabel("time (sec)")
    plt.ylabel("current (A)")
    plt.savefig("ilpj_order_h.png")

    


main()