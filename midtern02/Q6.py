#!/usr/bin/env python3
# Q6, 106061121, 莊裕嵐
# solving nonlinear equation
# Please solve the following equation:
#
# exp(x) + log(x) = 0.9
#
import numpy as np 
from ee4070 import *

def function(x) : return np.exp(x) + np.log(x) - 0.9

def main() :

    x0 = 1

    zero = Newton_root(x0, function, epsilon = 10 ** -9 )

    print(zero)

main()
