import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from src.thermo_lookup import get_property

g_H2 = get_property('H2', 1000, 'delta_gf')
g_C2H4 = get_property('C2H4', 1000, 'delta_gf')
g_CH4 = get_property('CH4', 1000, 'delta_gf')

delta_G = 2 * g_H2 +  g_C2H4 - 2 * g_CH4

Kp = np.exp((-delta_G)/(.008314*1000))

print("Equillibrium Coefficient = ", Kp)