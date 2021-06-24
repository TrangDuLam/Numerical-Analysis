
# Q3, 106061121, 莊裕嵐
# Solving ODE

import numpy as np
from math import *
from ee4070 import *
from matplotlib import pyplot as plt

def fwdEuler(h) :

    y1 = [1]
    y2 = [0]
    
    t = 0

    while t <= 3.5 : 

        y1_new = y1[-1] + h*y2[-1]
        y2_new = y2[-1] + h*(-0.1 * y2[-1] - y1[-1])

        #update altogether
        y1.append(y1_new)
        y2.append(y2_new)

        t += h
    
    return y1, y2

def bwdEuler(h) :

    y1 = [1]
    y2 = [0]

    t = 0

    while t <= 35 :

        A = np.zeros((2, 2), dtype = float)
        b = np.zeros((2,), dtype = float)

        A[0, 0] = 1
        A[0, 1] = -h
        A[1, 0] = h
        A[1, 1] = 1 + 0.1*h

        b[0] = y1[-1]
        b[1] = y2[-1]

        ynew = linSol(2, A, b)

        y1.append(ynew[0])
        y2.append(ynew[1])

        t += h

    return y1, y2

def trapezoidal(h) :

    y1 = [1]
    y2 = [0]

    t = 0

    while t <= 35 :

        A = np.zeros((2, 2), dtype = float)
        b = np.zeros((2,), dtype = float)

        A[0, 0] = 1
        A[0, 1] = -h/2
        A[1, 0] = h/2
        A[1, 1] = 1 + 0.1*h / 2

        b[0] = y1[-1] + h*y2[-1]/2
        b[1] = y2[-1] + h*(-0.1*y2[-1] - y1[-1])/2

        ynew = linSol(2, A, b)

        y1.append(ynew[0])
        y2.append(ynew[1])

        t += h

    return y1, y2




def main() :

    h = 0.1

    y1, y2 = fwdEuler(h)

    print("Forward Euler Method :")
    print("y1(0) = {:g}, y1(1) = {:g}, y1(2) = {:g}, y1(3) = {:g}".format(y1[0], y1[10], y1[20], y1[30]))
    print("y2(0) = {:g}, y2(1) = {:g}, y2(2) = {:g}, y2(3) = {:g}".format(y2[0], y2[10], y2[20], y2[30]))

    y1.clear()
    y2.clear()

    y1, y2 = bwdEuler(h)
    print("Backward Euler Method :")
    print("y1(0) = {:g}, y1(1) = {:g}, y1(2) = {:g}, y1(3) = {:g}".format(y1[0], y1[10], y1[20], y1[30]))
    print("y2(0) = {:g}, y2(1) = {:g}, y2(2) = {:g}, y2(3) = {:g}".format(y2[0], y2[10], y2[20], y2[30]))

    y1.clear()
    y2.clear()

    y1, y2 = trapezoidal(h)
    print("Trapezoidal Euler Method :")
    print("y1(0) = {:g}, y1(1) = {:g}, y1(2) = {:g}, y1(3) = {:g}".format(y1[0], y1[10], y1[20], y1[30]))
    print("y2(0) = {:g}, y2(1) = {:g}, y2(2) = {:g}, y2(3) = {:g}".format(y2[0], y2[10], y2[20], y2[30]))





main()

