#This script compute the Floquet Magnus expansion at order n
from stt_qtt import*
from oprts import*
import qutip as qt 
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from fd_interaction import*

#Hamiltonian for each half period are defined as H1 and H2
X = -0.5j*Hp
Y = -0.5j*Hm

#first term of the expansion
P1 = tau*(X + Y)
    
#P is the only term that has to be calculated explicitly in this program
def P(i, X, Y, tau):
    P_aux = 0.0
    for k in range(0, i + 1):
        P_aux += ((X**(i-k))*(Y**k))/(math.factorial(i-k)*math.factorial(k))
    return (tau**i)*P_aux

def Q(i, j):
    if j <= 2 or j >= i:
        return 0.0
    else:
        Qsum = 0.0
        for k in range(1, i - j + 2):
            Qsum += Q(k, 1)*Q(i - j,k - 1)
    return Qsum

def R(i):
    if i == 1:
        return P1 
    else:
        Qtot = 0.0
        for k in range(2, i + 1):
            Qtot += - Q(i,k)/math.factorial(k)      
        return P(i, X, Y, tau) - Qtot 

    