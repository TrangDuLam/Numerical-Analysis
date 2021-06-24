#!/usr/bin/env python3
# Q1, 106061121, 莊裕嵐
# polynomial roots
# Given a polynomial below, please find all its roots
#
# P(x) = x**3 -0.8 x**2 -0.35 x +0.15
#

import numpy as np 
from ee4070 import * 

def main() : 

    p = np.array([0.15, -0.35, -0.8, 1], dtype = complex)
    # permutation from low to high

    R = roots_sol(p, epsilon = 10 ** -9) #return sorted list

    print(np.round(R, 9))

    '''
    using np.round function to avoid 0.99999999~ issue
    The default tolerance of root solution function in ee4070 is 10 ** -12 
    '''


main()
