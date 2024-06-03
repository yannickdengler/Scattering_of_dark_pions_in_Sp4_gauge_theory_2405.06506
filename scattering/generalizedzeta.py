import numpy as np
import scipy.integrate as scpint

def y_elem_Z3(order):
    """
    Returns y in integers in 3 dimensions
    """
    y = []
    for i in [x for x in range(-order, order+1)]:
        for j in [x for x in range(-order, order+1) if abs(x)+abs(i) < order]:
            for k in [x for x in range(-order, order+1) if abs(x)+abs(i)+abs(j) < order]:
                y.append((i,j,k))
    return y

def y_squared(order):
    res = []
    for y in y_elem_Z3(order):
        res.append(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])
    return res
    
def term1(q2, order=10):
    """
    Calculation of the first term
    """
    res = 0
    for y_2 in y_squared(order):
        res += np.exp(-(y_2-q2))/(y_2-q2)
    return res

def factorial(x):
    """
    Factorial of x
    """
    res = 1
    for i in range(1, x+1):
        res = res*i
    return res

def term2(q2, order=23):
    """
    Calculation of the second term
    """
    res = 0
    for k in range(order):
        res += (q2**k)/(factorial(k)*(k-0.5))
    return np.pi**(3/2.)*res

def exp_3(u, order):
    res = 0
    for y_2 in y_squared(order):
        res += np.exp(-np.pi**2*y_2/u)
    return res-1

def e_u_q2(u, q2, order):
    return np.exp(u*q2)*(np.pi/u)**(3/2.)*exp_3(u, order)

def term3(q2, order=7):
    """
    Calculation of the third term
    """
    return scpint.quad(func=e_u_q2, a=0, b=1, args=(q2, order,))[0]

def Zeta(q2, orders=(10,23,7)):
    """
    The Zeta function calculated from q2, orders is a scale for the precision. To which orders are the terms calculated.
    """
    return (term1(q2, orders[0]) + term2(q2, orders[1]) + term3(q2, orders[2]))/(np.sqrt(4*np.pi))