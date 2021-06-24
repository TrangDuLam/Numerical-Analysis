
# Q2, 106061121, 莊裕嵐
# Solving nonlinear systems

import numpy as np
from math import *
from ee4070 import *

def F01(x) :

    F_arr = np.empty((2,), dtype = float)
    F_arr[0] = x[0]**2 + 3*x[0] - x[1] - 3.81
    F_arr[1] = -x[0] + 3*x[1]**2 + x[1] -1.07 
    return F_arr

def J01(x) :

    J = np.empty((2, 2), dtype = float)
    J[0, 0] = 2*x[0] + 3
    J[0, 1] = -1
    J[1, 0] = -1
    J[1, 1] = 6*x[1] + 1

    return J

def F02(x) :

    F_arr = np.empty((2,), dtype = float)
    F_arr[0] = x[0]**3 + x[0] - x[1] + 1.375
    F_arr[1] = -x[0] + x[1]**3 + 2 * x[1] - 11.5

    return F_arr

def J02(x) :

    J = np.empty((2,2), dtype = float)
    J[0, 0] = 3 * x[0]**2 + 1
    J[0, 1] = -1
    J[1, 0] = -1
    J[1, 1] = 3 * x[1]**2 + 2

    return J

def F03(x) :

    F_arr = np.empty((2,), dtype = float)
    F_arr[0] = np.sqrt(x[0] + 1) - x[1] - 0.60384
    F_arr[1] = np.sqrt(x[1] + 2) - x[0] - 0.943168

    return F_arr

def J03(x) :

    J = np.empty((2,2), dtype = float)
    J[0, 0] = 1/(2 * np.sqrt(x[0] + 1))
    J[0, 1] = -1
    J[1, 0] = -1
    J[1, 1] = 1/(2 * np.sqrt(x[1] + 2))

    return J

def F04(x) :

    F_arr = np.empty((2,), dtype = float)
    F_arr[0] = 1/(x[0] + 1) - x[1] + 5/3
    F_arr[1] = -x[0] + 2/(x[1] + 1) + 4/3

    return F_arr

def J04(x) :

    J = np.empty((2,2), dtype = float)
    J[0, 0] = 1 / ((x[0] + 1)**2)
    J[0, 1] = -1
    J[1, 0] = -1
    J[1, 1] = -2 / ((x[1] + 1)**2)

    return J

def F05(x) :

    F_arr = np.empty((2,), dtype = float)
    F_arr[0] = np.sin(x[0]) - x[1] + 0.208793
    F_arr[1] = -x[0] + np.cos(x[1]) + 0.646404

    return F_arr

def J05(x) :

    J = np.empty((2,2), dtype = float)
    J[0, 0] = np.sin(x[0])
    J[0, 1] = -1
    J[1, 0] = -1
    J[1, 1] = -np.cos(x[1])

    return J


def main() :

    x0 = np.ones((2,), dtype = float)
    sol01 = Newton_N_dim(x0, 2, F01, J01)
    print("x = {:g}, y = {:g}".format(sol01[0], sol01[1]))

    x0 = np.ones((2,), dtype = float)
    sol02 = Newton_N_dim(x0, 2, F02, J02)
    print("x = {:g}, y = {:g}".format(sol02[0], sol02[1]))

    x0 = np.ones((2,), dtype = float)
    sol03 = Newton_N_dim(x0, 2, F03, J03)
    print("x = {:g}, y = {:g}".format(sol03[0], sol03[1]))

    x0 = np.ones((2,), dtype = float)
    sol04 = Newton_N_dim(x0, 2, F04, J04)
    print("x = {:g}, y = {:g}".format(sol04[0], sol04[1]))

    x0 = np.ones((2,), dtype = float)
    sol05 = Newton_N_dim(x0, 2, F05, J05)
    print("x = {:g}, y = {:g}".format(sol05[0], sol05[1]))


main()

