#!/usr/bin/env python3
# Q10, 106061121, 莊裕嵐
# solving nonlinear equation
# Please solve the following equation:
#
# log(1.0 / x) = -1.9
#

import numpy as np
from ee4070 import *

def function(x) : return np.log(1.0 / x) + 1.9

def main() :

    x0 = 1
    '''
    first try = 1 //success
    second try = 0.8 //success
    '''

    zero = Newton_root(x0, function, epsilon = 10 ** -9)

    print(zero)

main()