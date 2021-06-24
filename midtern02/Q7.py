#!/usr/bin/env python3
# Q7, 106061121, 莊裕嵐
# solving nonlinear equation
# Please solve the following equation:
#
# log(x) * log(x) = 1.8
#

import numpy as np 
from ee4070 import *

def function(x) : return np.log(x) * np.log(x) -1.8

def main() :

    x0 = 0.5  #initial guessing 
    '''
    first try = 1 //failed => derivative = 0
    second try = 0 //failed => invalid value 
    third try = 0.5 //success
    forth try = 0.01 //success

    '''

    zero = Newton_root(x0, function, epsilon = 10 ** -9)

    print(zero)

main()
