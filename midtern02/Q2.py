#!/usr/bin/env python3
# Q2, 106061121, 莊裕嵐
# polynomial roots
# Given a polynomial below, please find all its roots
#
# P(x) = x**4 -1.2 x**3 -1.95 x**2 +0.55 x +0.3
#

import numpy as np 
from ee4070 import * 

def main() : 

    p = np.array([0.3, 0.55, -1.95, -1.2, 1], dtype = complex)

    R = roots_sol(p, epsilon = 10 ** -9)

    print(np.round(R, 9))


main()

