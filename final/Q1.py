
# Q1, 106061121, 莊裕嵐
# Solving nonlinear equations

import numpy as np
from math import *
from ee4070 import *

def fun01(x) : return np.exp(x) + np.exp(np.exp(x)) + np.exp(np.exp(np.exp(x))) - 20 
def fun02(x) : return np.log(x) + np.log(np.log(x)) + np.log(np.log(np.log(x))) - 3 
def fun03(x) : return np.log(x) + np.log(x) * np.log(x) + np.log(x) * np.log(x) * np.log(x) - 5
def fun04(x) : 
    return (1/(np.log(x))) + (1/(np.log(x) * np.log(x))) + (1/(np.log(x)*np.log(x)*np.log(x))) - 0.5
def fun05(x) : return np.sin(x) + np.sin(np.sin(x)) + np.sin(np.sin(np.sin(x))) - 0.9
def fun06(x) : return np.cos(x) + np.cos(np.cos(x)) + np.cos(np.cos(np.cos(x))) - 2




def main() :

    x1 = 0
    root1 = Newton_root(x1, fun01)
    print("answer = {:g}".format(root1))

    x2 = 5 #1st = 0.1, 2nd = 1
    root2 = Newton_root(x2, fun02)
    print("answer = {:g}".format(root2))

    x3 = 1 
    root3 = Newton_root(x3, fun03)
    print("answer = {:g}".format(root3))

    x4 = 3 #1st = 0.1, 2nd = 1
    root4 = Newton_root(x4, fun04)
    print("answer = {:g}".format(root4))

    x5 = 0
    root5 = Newton_root(x5, fun05)
    print("answer = {:g}".format(root5))

    x6 = np.pi/2 #1st = 0, 2nd = pi
    root6 = Newton_root(x6, fun06)
    print("answer = {:g}".format(root6))

main()
    