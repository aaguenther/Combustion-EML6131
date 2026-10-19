import math as m
import numpy as np
from matplotlib import pyplot as plt

T1 = 1000 # K
T2 = 2000 # K
k = 1.380649e-23 # J/K
N_A = 6.02214e+23 # mol^-1
# H properties
M_H = 1.008 # kg/kmol
# O2 properties
M_O2 = 31.999 # kg/kmol

m_H = (M_H/1000) / N_A # kg
m_O2 = (M_O2/1000) / N_A # kg

# Reduced mass

m_r = (m_H * m_O2) / (m_H + m_O2) # kg

# Beta values

beta_T1 = m.sqrt(m_r / (2 * k * T1))
beta_T2 = m.sqrt(m_r / (2 * k * T2))

# Create array for values of c_r that exceed E*

c_r = np.linspace(11969.3, 100000, 1000)

# Distribution function

def dist_function(beta, c):
    return (4 / m.pi**0.5) * beta**3 * c**2 * m.exp(-beta**2 * c**2)

f_cr_T1 = []
f_cr_T2 = []

for i in c_r:
    f_cr_T1.append(dist_function(beta_T1, i))
    f_cr_T2.append(dist_function(beta_T2, i))

# Probability of pairs occupying energy states greater than E* at each temp given by

F_T1 = np.trapezoid(f_cr_T1, c_r)
F_T2 = np.trapezoid(f_cr_T2, c_r)
print(F_T1)
print(F_T2)
