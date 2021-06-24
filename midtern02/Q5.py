#!/usr/bin/env python3
# Q5, 106061121, 莊裕嵐
# solving nonlinear equation
# Please solve the following equation:
#
# exp(x) + exp(2*x) = 3
#
import numpy as np 
from ee4070 import * 

def function(x) : return np.exp(x) + np.exp(2*x) - 3

def main() :

    x0 = 1

    zero = Newton_root(x0, function, epsilon = 10 **-9) #default tolerance = 10 ** -12 in ee4070.py

    print(zero)

main()