import math as m

# GIVEN

P = 3.2e-5 # (kPa)
T = 195 # (K)
L_char = 0.5 # (m)
d = 4.17e-10 # (m)

# Constants

R_u = 8314.46 # (J/kmol * K)
M_N2 = 28.013 # (kg/kmol)
k_b = 1.3806488e-23 # (J/K)

mean_thermal_speed = m.sqrt((8 * R_u * T) / (m.pi * M_N2))
mean_relative_velocity = m.sqrt(2) * mean_thermal_speed
collision_cross_section = m.pi * d**2
number_density = (P * 1000) / (k_b * T)

mean_free_path = mean_thermal_speed / (mean_relative_velocity * collision_cross_section * number_density)
Kn = mean_free_path / L_char

if Kn < 0.1:
    print("Continuum satisfied")
    print("Kn =", Kn)
else:
    print("Continuum not satisfied")
    print("Kn =", Kn)
