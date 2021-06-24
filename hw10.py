
# EE4070 Numerical Analysis
# HW10 Numerical Integration
# ID : 106061121, Name : Yulan Chuang
# Date : May 18, 2021


import numpy as np 

def function(x) :
    return np.exp(x) + np.exp(-x)

def first_order_nt(Y, h) :

    w = np.array([1/2, 1/2])
    
    i = 0
    integral = 0
    while(i < len(Y) - 1) :
        integral += np.sum(Y[i:i+2] * w)
        i += 1

    integral *= h


    return integral

def second_order_nt(Y, h) :
    
    w = np.array([1/3, 4/3, 1/3])
    
    i = 0
    integral = 0
    while(i < len(Y) - 2) :
        integral += np.sum(Y[i:i+3] * w)
        i += 2

    integral *= h


    return integral

def third_order_nt(Y, h) :
    
    w = np.array([3/8, 9/8, 9/8, 3/8])
     
    i = 0
    integral = 0
    while(i < len(Y) - 3) :
        integral += np.sum(Y[i:i+4] * w)
        i += 3

    integral *= h


    return integral

def fourth_order_nt(Y, h) :
    
    w = np.array([14/45, 64/45, 24/45, 64/45, 14/45])
     
    i = 0
    integral = 0
    while(i < len(Y) - 4) :
        integral += np.sum(Y[i:i+5] * w)
        i += 4

    integral *= h


    return integral

def sixth_order_nt(Y, h) :

    w = np.array([41/140, 216/140, 27/140, 272/140, 27/140, 216/140, 41/140])
     
    i = 0
    integral = 0
    while(i < len(Y) - 6) :
        integral += np.sum(Y[i:i+7] * w)
        i += 6

    integral *= h


    return integral

def main() :

    N = int(input("N = ")) #given by assignment

    exact_integral = np.exp(2) - np.exp(-2) #exact value of the integral

    x_i = np.linspace(0, 2, N + 1)  #add 1  (start, end, portion + 1) interval seperation

    y_i = function(x_i) #return the function value
    h = 2 / N #find the length of interval


    integral_1 = first_order_nt(y_i, h)
    print("Numerical of the first order : ",integral_1)
    print("Approximate error : ",np.abs(exact_integral - integral_1))

    integral_2 = second_order_nt(y_i, h)
    print("Numerical of the second order : ", integral_2)
    print("Approximate error : ", np.abs(exact_integral - integral_2))

    integral_3 = third_order_nt(y_i, h)
    print("Numerical of the third order : ", integral_3)
    print("Approximate error : " ,np.abs(exact_integral - integral_3))

    integral_4 = fourth_order_nt(y_i, h)
    print("Numerical of the fourth order : ", integral_4)
    print("Approximate error : " ,np.abs(exact_integral - integral_4))

    integral_6 = sixth_order_nt(y_i, h)
    print("Numerical of the sixth order : ", integral_6)
    print("Approximate error : " , np.abs(exact_integral - integral_6))


main()
