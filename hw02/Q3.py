import math as m
import numpy as np
from matplotlib import pyplot as plt

# Part a

# GIVEN

T1 = 300 # K
T2 = 2000 # K
k = 1.380649e-23 # J/K
N_A = 6.02214e+23 # mol^-1
# He properties
M_He = 4.003 # kg/kmol
# N2 properties
M_N2 = 28.013 # kg/kmol

# Molecular masses

m_He = (M_He/1000) / N_A # kg
m_N2 = (M_N2/1000) / N_A # kg

beta_N2_T1 = (m_N2 / (2 * k * T1)) ** 0.5
beta_N2_T2 = (m_N2 / (2 * k * T2)) ** 0.5
beta_He_T1 = (m_He / (2 * k * T1)) ** 0.5
beta_He_T2 = (m_He / (2 * k * T2)) ** 0.5

def dist_function(beta, c_prime):
    return (4 / m.pi**0.5) * beta**3 * c_prime**2 * m.exp(-beta**2 * c_prime**2)

# Create an array for c_prime

c_prime = np.linspace(0, 10000, 10000)

# Create an empty array for the distribution function

f_c_prime_N2_T1 = []
f_c_prime_N2_T2 = []
f_c_prime_He_T1 = []
f_c_prime_He_T2 = []

for i in c_prime:
    f_c_prime_N2_T1.append(dist_function(beta_N2_T1, i))
    f_c_prime_N2_T2.append(dist_function(beta_N2_T2, i))
    f_c_prime_He_T1.append(dist_function(beta_He_T1, i))
    f_c_prime_He_T2.append(dist_function(beta_He_T2, i))

# Most probable thermal speeds:

c_m_prime_N2_T1 = 1/beta_N2_T1
print(c_m_prime_N2_T1)
c_m_prime_N2_T2 = 1/beta_N2_T2
print(c_m_prime_N2_T2)
c_m_prime_He_T1 = 1/beta_He_T1
print(c_m_prime_He_T1)
c_m_prime_He_T2 = 1/beta_He_T2
print(c_m_prime_He_T2)

# Create plots

# Create lines for c'
plt.plot(c_prime, f_c_prime_N2_T1, label="N2 @ 300 K", color = "blue")
plt.plot(c_prime, f_c_prime_N2_T2, label="N2 @ 2000 K", color = "blue", linestyle = "--")
plt.plot(c_prime, f_c_prime_He_T1, label="He @ 300 K", color = "red")
plt.plot(c_prime, f_c_prime_He_T2, label="He @ 2000 K", color = "red", linestyle = "--")

# Add points for cm'
plt.scatter(c_m_prime_N2_T1, dist_function(beta_N2_T1, c_m_prime_N2_T1), color = "blue")
plt.scatter(c_m_prime_N2_T2, dist_function(beta_N2_T2, c_m_prime_N2_T2), color = "blue")
plt.scatter(c_m_prime_He_T1, dist_function(beta_He_T1, c_m_prime_He_T1), color = "red")
plt.scatter(c_m_prime_He_T2, dist_function(beta_He_T2, c_m_prime_He_T2), color = "red")

plt.xlabel("c'")
plt.ylabel("f'_c")
plt.title("Distribution Function vs Thermal Speed")
plt.legend()
plt.show()

# Part c


# Numerical solutions for mean and RMS thermal speed

mean_thermal_speed_N2_T1 = np.trapezoid(c_prime * f_c_prime_N2_T1, c_prime)
print(f"N2 at 300 K - Numerical mean speed: {mean_thermal_speed_N2_T1:.2f} m/s")

rms_thermal_speed_N2_T1 = m.sqrt(np.trapezoid(c_prime**2 * f_c_prime_N2_T1, c_prime))
print(f"N2 at 300 K - Numerical RMS speed: {rms_thermal_speed_N2_T1:.2f} m/s")

mean_thermal_speed_N2_T2 = np.trapezoid(c_prime * f_c_prime_N2_T2, c_prime)
print(f"N2 at 2000 K - Numerical mean speed: {mean_thermal_speed_N2_T2:.2f} m/s")

rms_thermal_speed_N2_T2 = m.sqrt(np.trapezoid(c_prime**2 * f_c_prime_N2_T2, c_prime))
print(f"N2 at 2000 K - Numerical RMS speed: {rms_thermal_speed_N2_T2:.2f} m/s")

mean_thermal_speed_He_T1 = np.trapezoid(c_prime * f_c_prime_He_T1, c_prime)
print(f"He at 300 K - Numerical mean speed: {mean_thermal_speed_He_T1:.2f} m/s")

rms_thermal_speed_He_T1 = m.sqrt(np.trapezoid(c_prime**2 * f_c_prime_He_T1, c_prime))
print(f"He at 300 K - Numerical RMS speed: {rms_thermal_speed_He_T1:.2f} m/s")

mean_thermal_speed_He_T2 = np.trapezoid(c_prime * f_c_prime_He_T2, c_prime)
print(f"He at 2000 K - Numerical mean speed: {mean_thermal_speed_He_T2:.2f} m/s")

rms_thermal_speed_He_T2 = m.sqrt(np.trapezoid(c_prime**2 * f_c_prime_He_T2, c_prime))
print(f"He at 2000 K - Numerical RMS speed: {rms_thermal_speed_He_T2:.2f} m/s")


# Closed-form solutions for mean and RMS thermal speed

mean_thermal_speed_closed_N2_T1 = 2 / (m.pi**0.5 * beta_N2_T1)
print(f"N2 at 300 K - Analytical mean speed: {mean_thermal_speed_closed_N2_T1:.2f} m/s")

rms_thermal_speed_closed_N2_T1 = (1.5**0.5) / beta_N2_T1
print(f"N2 at 300 K - Analytical RMS speed: {rms_thermal_speed_closed_N2_T1:.2f} m/s")

mean_thermal_speed_closed_N2_T2 = 2 / (m.pi**0.5 * beta_N2_T2)
print(f"N2 at 2000 K - Analytical mean speed: {mean_thermal_speed_closed_N2_T2:.2f} m/s")

rms_thermal_speed_closed_N2_T2 = (1.5**0.5) / beta_N2_T2
print(f"N2 at 2000 K - Analytical RMS speed: {rms_thermal_speed_closed_N2_T2:.2f} m/s")

mean_thermal_speed_closed_He_T1 = 2 / (m.pi**0.5 * beta_He_T1)
print(f"He at 300 K - Analytical mean speed: {mean_thermal_speed_closed_He_T1:.2f} m/s")

rms_thermal_speed_closed_He_T1 = (1.5**0.5) / beta_He_T1
print(f"He at 300 K - Analytical RMS speed: {rms_thermal_speed_closed_He_T1:.2f} m/s")

mean_thermal_speed_closed_He_T2 = 2 / (m.pi**0.5 * beta_He_T2)
print(f"He at 2000 K - Analytical mean speed: {mean_thermal_speed_closed_He_T2:.2f} m/s")

rms_thermal_speed_closed_He_T2 = (1.5**0.5) / beta_He_T2
print(f"He at 2000 K - Analytical RMS speed: {rms_thermal_speed_closed_He_T2:.2f} m/s")



