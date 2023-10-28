import numpy as np
import math
from scipy.special import factorial

def Q(beta, alpha, R, Q_values):
    if beta <= 2 or beta >= alpha:
        return 0.0
    else:
        result = 0.0
        for gamma in range(1, alpha - beta + 2):
            result += R[gamma - 1] * Q(beta - 1, alpha - gamma, R, Q_values)
        return result

def R(alpha, R_values, Q):
    if alpha == 1:
        return R_values[0]
    else:
        result = R[alpha - 1]
        beta_values = np.arange(2, alpha + 1)
        # Compute factorials for each value in beta_values
        factorials = factorial(beta_values)
        result -= np.sum(1 / factorials * Q(beta_values, alpha, R_values, Q_values))
        return result

# Define the value of alpha
alpha = 5

# Initialize the R and Q arrays
R_values = np.zeros(alpha)
Q_values = np.zeros((alpha, alpha))

# Set P_1 = R_1
P_1 = 1.0  # You can change this to your desired value
R_values[0] = P_1

# Calculate R_alpha
result_R_alpha = R(alpha, R_values, Q_values)

print(f'R_{alpha} = {result_R_alpha}')
