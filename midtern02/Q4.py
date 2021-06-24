#!/usr/bin/env python3
# Q4, 106061121, 莊裕嵐
# polynomial roots
# Given a polynomial below, please find all its roots
#
# P(x) = x**6 -3.6 x**5 +2.21 x**4 +2.254 x**3 -0.348 x**2 -0.376 x -0.048
#
import numpy as np 
from ee4070 import * 

def main() : 

    p = np.array([-0.048, -0.376, -0.348, 2.254, 2.21, -3.6, 1.], dtype = complex)

    R = roots_sol(p, epsilon = 10 ** -9)

    print(np.round(R, 9))


main()
