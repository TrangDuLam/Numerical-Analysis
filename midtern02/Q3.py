#!/usr/bin/env python3
# Q3, 106061121, 莊裕嵐
# polynomial roots
# Given a polynomial below, please find all its roots
#
# P(x) = x**5 -5.3 x**4 +8.74 x**3 -6.132 x**2 +1.908 x -0.216
#
import numpy as np 
from ee4070 import * 

def main() : 

    p = np.array([-0.216, 1.908, -6.132, 8.74, -5.3, 1], dtype = complex)

    R = np.array(roots_sol(p, epsilon = 10 ** -12))

    print(np.round(R, 9))


main()