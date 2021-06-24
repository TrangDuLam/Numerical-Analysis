#!/usr/bin/env python3
# Q8, 106061121, 莊裕嵐
# solving nonlinear equation
# Please solve the following equation:
#
# exp(x) * log(x) = -2.7
#

import numpy as np
from ee4070 import *

def function(x) : return np.exp(x) * np.log(x) + 2.7

def main() :

    x0 = 1

    zero = Newton_root(x0, function, epsilon = 10 ** -9)

    print(zero)

main()
