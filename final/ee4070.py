import numpy as np
import sys 
import os

import bisect #spline interpolation need

def norm1(n, x) :
    
    return np.sum(np.abs(x))

def norm2(n, x) :

    return np.sqrt(x @ x)

def normInf(n, x) : 

    return np.abs(x).max()

def luFact(n, A) :  #(dimension, target matrix)

  U = A.copy()  #copy the matrix A as U
  L = np.array(np.identity(n)) #let the matrix L be the identity matrix
  #accroding to the format of L


  for i in range(n):  #for every row in the matrix
    for k in range(i+1, n) : #for each element in the given row
      c = - (float(U[k, i])/U[i, i])  #take the ration of the diagonal element
      U[k, i:] = U[k, i:] + c*U[i, i:]  #do the addition in the row element in U
      L[k:, i] = L[k:, i] - c*L[k:, k]  #do the subtraction in the row element in L

  return L, U

def fwdSubs(n, A, b) : #forward substitution
  
  b = np.array(b)  #declare the data type as array

  y = b.copy()  #copy the matrix b
  L = luFact(n, A)[0] #call the matrix L

  for i in range(n) : #for every element in the matrix b (i-th row in matrix L)
     for k in range(0, i) :  #for every element i-th row 
       y[i] = y[i] - L[i, k]*y[k] #subtraction 

  
  return y

def bwdSubs(n, A, y) :

  y = np.array(y) #declare the data type as array

  x = y.copy() #copy the matrix y
  U = luFact(n, A)[1] #call the matrix U

  for i in range(n-1, -1, -1) : #for every element in the matrix b (i-th row in matrix L)
                                #-1 : the backward operation 
    for k in range(n-1, i, -1) : #for every element i-th row
      x[i] -= U[i, k]*x[k]  #subtraction

    x[i] = float(x[i])/U[i, i]  #division to the diagonal element

  return x

def linSol(n, A, b) :

    y = fwdSubs(n, A, b)
    x = bwdSubs(n, A, y)

    return x

def Gram_Schimdt(A) :
    m, n = A.shape
    G = np.copy(A)
    G[:, 0] = G[:, 0]/(norm2(n, G[:, 0]))

    for k in range(1,n) :
        for i in range(k-1):
            G[:, k] = G[:, k] - (np.dot(G[:, k], G[:, i])*G[:, i])/(np.dot(G[:, i], G[:, i]))
        
        G[:, k] = G[:, k]/(norm2(n, G[:, k]))


    return G

def Jacobi(n, A, b, maxIter = 10**6, tol = 10**(-7), enorm = norm2) :

    x = np.zeros(n)  #initial value of x
    x_tmp = x.copy() #copy
    answer = list([False]) #set the return value 

    D = np.diag(A) #extract the diagonal elements
    R = A - np.diagflat(D) #subtract the input matrix A by the diagonal elements

    for itr in range(maxIter) :

        for i in range(n) :
            x_tmp[i] = (b[i] - R[i, :] @ x)/D[i] #computing process

        error = enorm(n, x_tmp - x) #computing the error

        if(error < tol) : #convergence reaches
          answer[0] = True 
          answer.append(x_tmp)
          break
        
        else : #renew the matrix x 
            x = np.copy(x_tmp)
            pass


    return answer

def GS(n, A, b, maxIter = 10**6, tol = 10**(-7), enorm = norm2):     # Gauss-Seidel 

    x = np.zeros(n) #initial value of x
    x_tmp = np.copy(x) #copy
    answer = list([False]) #set the return value as False initially

    D = np.diag(A) #extract the diagonal elements
    R = A - np.diagflat(D) #subtract the input matrix A by the diagonal elements

    for iter in range(maxIter) :

        for i in range(n) :  #computing process
            sigma = 0 

            for j in range(i):
                sigma += R[i, j]*x_tmp[j]
            for j in range(i+1, n):
                sigma += R[i, j]*x[j]
        
            x_tmp[i] = (b[i] - sigma)/D[i]

        error = enorm(n, x_tmp - x) #computing the error
        
        if(error < tol) : #convergence reaches
            answer[0] = True
            answer.append(x_tmp)
            break
        
        else :
            x = np.copy(x_tmp) #renew the matrix x
            pass

    return answer

def SGS(n, A, b, maxIter = 10**6, tol = 10**(-7), enorm = norm2): 

    x = np.zeros((n,)) #initial value of x
    x_tmp = np.copy(x) #copy
    answer = list([False]) #set the return value as False initially


    for iter in range(maxIter) :
        
        for i in range(n) : #computing process : forward SGS
            m = x * A[i, :]
            sigma = np.sum(m)
            sigma -= m[i]
            x_tmp[i] = (b[i] - sigma)/A[i][i]

        for i in reversed(range(n)) :
            m = x * A[i, :]
            sigma = np.sum(m)
            sigma -= m[i]
            x_tmp[i] = (b[i] - sigma)/A[i][i]

        error = enorm(n, x_tmp - x)  #computing the error
        
        if(error < tol) : #convergence reaches
            answer[0] = True
            answer.append(x_tmp)
            break
        
        else :
            x = np.copy(x_tmp) #renew the matrix x
            pass

    return answer

def CG(n, A, b, maxIter = 10**6, tol = 10**(-7)):  #Conjugate Gradient Decend

    ans = list([False]) #assume the convergence as 'false'

    #initial condition
    x = np.zeros((n,)) #initial condition
    r = b - A @ x #calculate r0
    p = np.copy(r)  #equality
    nold = np.dot(r, r) #calculate the initial norm

    iter = 0 #denote the iteration time

    for iter in range(maxIter) : #
        Ap = A @ p #computing A times p ar first
        alpha = nold/(np.dot(p, Ap)) #computing alpha
        x = x + alpha * p #renew x
        r = r - alpha * Ap #renew r
        nnew = np.dot(r, r) #computing new norm
        err = np.sqrt(nnew/n) #computing the error

        if err < tol :
            ans[0] = True #if the condition reaches
            ans.append(iter+1) #append the total iteration time
            ans.append(x) #append the final solution
            #the foramt assigned by the instruction
            break #convergence reached

        p = r + (nnew/nold)*p #renew p
        nold = nnew #renew norm


    return ans #return final answer

def QR(n, A) :

    A = A.astype(float) #declare datatype
    Q = np.copy(A) #copy A to Q
    R = np.zeros((n, n)) #declare a zero matrix

    r00 = norm2(n, Q[:, 0]) #norm of column vector 0
    R[0, 0] = r00 #value assigned
    Q[:, 0] = Q[:, 0]/r00 #normalization

    for j in range(1, n) : # row loop

        for i in range(0, j) : 
            rij = Q[:, i] @ Q[:, j] # inner product
            R[i, j] = rij # value assigned
            Q[:, j] = Q[:, j] - rij*Q[:, i] # value computing
        

        rjj = norm2(n, Q[:, j]) #take norm of the column vector j
        R[j, j] = rjj #value assigned
        Q[:, j] = Q[:, j]/rjj # normalization

    return Q, R


def Lagrange(x, xdata, ydata) :  #(x-axis for fn generating, interpolation x, interpolation y)
                                 # len. xdata = len. ydata

    fi = np.zeros(len(x))  #initialize interpolation function

    for j in range(len(x)) :
    
        for i in range(len(xdata)) :
            L = 1

            for k in range(len(xdata)) :
                if i != k :
                    L *= (x[j] - xdata[k])/(xdata[i] - xdata[k])

            fi[j] += ydata[i]*L


    return fi


def splineM(N, X, Y) :
    
    h = np.diff(X)


    mu = [h[i] / (h[i] + h[i + 1]) for i in range(N - 2)] + [0] #using the attribute of python list to generate mu
    mu = np.array(mu, dtype = float)  #reform as numpy array
    twos = [2] * N #diagonal are 2s
    twos = np.array(twos, dtype = float)
    lam = [0] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(N - 2)] #using the attribute of python list to generate lamda
    lam = np.array(lam, dtype = float)  #reform as numpy array

    d = [0] + [6 * ((Y[i + 1] - Y[i]) / h[i] - (Y[i] - Y[i - 1]) / h[i - 1]) / (h[i] + h[i-1]) for i in range(1, N - 1)] + [0] ##using the attribute of python list to generate d
    
    lam_p = np.append(lam, float(0)) #append the boundary condition of lamda_n
    d_p = np.zeros((N, ), dtype= float) #declare the temporary vector d_p
    M = np.zeros((N, ), dtype= float) #declare the moment vector



    #Thomas algorithm process
    lam_p[0] = lam[0] / twos[0]
    d_p[0] = d[0] / twos[0]

    for i in range(1, N) :
        lam_p[i] = lam_p[i] / (twos[i] - lam_p[i - 1] * mu[i - 1])
        d_p[i] = (d[i] - d_p[i - 1] * mu[i - 1]) / (twos[i] - lam_p[i - 1] * mu[i - 1])

    M[-1] = d_p[-1]
    for i in range(N - 2, -1, -1):
        M[i] = d_p[i] - lam_p[i] * M[i + 1]

    return M

def spline(N, X, Y, M, x) :

    h = np.diff(X) #take the h

    coefficients = [[(M[i+1]-M[i])*h[i]*h[i]/6, M[i]*h[i]*h[i]/2, (Y[i+1] - Y[i] - (M[i+1]+2*M[i])*h[i]*h[i]/6), Y[i]] for i in range(N-1)] #using the attribute of python list to generate interpolation coefficients

    def polt_spline(val) : #define a function specified for the plotting process
        idx = min(bisect.bisect(X, val)-1, N-2) #search where the x value should be classified upon built-in binary search function  
        z = (val - X[idx]) / h[idx] #assign a new variable for convenience
        C = coefficients[idx] #return the interpolation corresponding to the given index  
        return (((C[0] * z) + C[1]) * z + C[2]) * z + C[3] #return the interpolation value

    return polt_spline #return the plot function as final output


# ****   Integral tools ****

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


# *** polynomial roots ****

def derivatives(function, point) :  

    h = 10 ** -12
    return (function(point + h) - function(point - h))/(2* h)

def chord(a, b, function, IterMax = 10 ** 6, epsilon = 10 ** -12) : 

    k = 0
    err = 1 + epsilon
    x = b
    g = (function(b) - function(a))/(b - a)

    while k < IterMax and err > epsilon : 
        x -= function(x)/g
        err = np.abs(function(x))

        k += 1

    return x

def regular_falsi(a, b, function, IterMax = 10 ** 6, epislon = 10 ** -12) : 

    k = 0 
    err = 1 + epislon 

    while k < IterMax and err > epislon :
        x = a - function(a) * ((b -  a)/ (function(b) - function(a)))

        if function(x) * function(a) <= 0 : 
            b = x
        else : 
            a = x 

        err = np.abs(function(x))
        k += 1

    return x

def secant(initial, x_0, function, IterMax = 10 ** 6, epislon = 10 ** -12) :  #x_-1, x_0

    k = 0
    err = 1 + epislon

    x_m1 = initial
    x = x_0 

    while k < IterMax and err > epislon :
        x -= function(x)*((x - x_m1)/(function(x) - function(x_m1)))

        k += 1
        err = np.abs(function(x))

    return x

def Newton_root(x, function, IterMax = 10 ** 6, epsilon = 10** -12) :

    k = 0 
    err = 1 + epsilon

    while k < IterMax and err > epsilon : 
        x -= (function(x))/(derivatives(function, x))

        k += 1
        err = np.abs(function(x))

    return x

def roots_sol(A, IterMax = 10 ** 6, epsilon = 10 ** -12) : # A : coefficient array of polynomial

    #coefficient permutation is a_0, a_1 ~ a_n (nth order polynomial)

    '''
    order of function = n 
    len(A) = n - 1

    b : -1 ~ n-1 => len(b) = n + 2 for indexing  
    c : -1 ~ n-2 => len(c) = n + 1 for indexing
    '''
    n = len(A) - 1 #order
    roots = list()
    b = np.empty((n + 2,), dtype = complex)
    c = np.empty((n + 1,), dtype = complex)
    x0 = np.max(np.abs(A)) + 1  #initial guess upon the Cauchy's rule

    x = x0

    while n >= 1 :

        err = 1 + epsilon 
        k = 0 

        while err >= epsilon and k < IterMax :
            b[n] = A[n]
            c[n - 1] = b[n]

            for j in range(n - 1, -1, -1) :   #(-1 ~ n-2) => (0 ~ n - 1)
                b[j] = A[j] + x * b[j + 1]
            for j in range(n - 2, -1, -1) :   #(-1 ~ n-3) => (0 ~ n - 2)
                c[j] = b[j + 1] + x * c[j + 1]

            f = b[0]
            f_p = c[0]
            x -= f/f_p
            err = np.abs(f)
            k += 1
        
        roots.append(x)

        for i in range(n) : 
            A[i] = b[i+1]

        #x = roots[-1]
        

        n -= 1

    
    return roots


def Newton_N_dim(x0, N, F_arr, dF_mat, IterMax = 10 ** 6, epsilon = 10 ** -12) :

    k = 0
    err = 1 + epsilon
    ''' 
        vector of F and the derivative dF should be pre-defined in the main function
        

    '''

    x = x0
    F = F_arr(x)
    df = dF_mat(x)

    while k < IterMax and err > epsilon :

        x_diff = linSol(N, df, -F)
        x += x_diff
        F = F_arr(x)
        df = dF_mat(x)
        err = norm2(2, x_diff)
        k +=1

    return x

def Radd(A, rp, rn, g):
	''' Add resistor with conductance g to circuit matrix A
        rp: positive node; rn: negative node '''
	A[rp, rp] += g
	A[rn, rn] += g
	A[rp, rn] -= g
	A[rn, rp] -= g

def Vsrc(n, A, b, vp, V):
	''' Add grounded votage source to circuit matrix A and RHS b
        vp: positive node for the voltage source
        V: voltage value
		This creates a asymmetric matrix'''
	for i in range(n):
		A[vp, i] = 0
	A[vp, vp] = 1
	b[vp] = V

def VsrcSym(n, A, b, vp, V):
	''' Add grounded votage source to circuit matrix A and RHS b
        vp: positive node for the voltage source
        V: voltage value
		This creates a symmetric matrix'''
	Vsrc(n, A, b, vp, V)
	for i in range(n):
		if A[i, vp] != 0 and i != vp:
			b[i] -= A[i, vp] * V
			A[i, vp] = 0
