import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from src.thermo_lookup import get_property
from scipy.optimize import brentq

# Stoichiometric case...

# Energy balance...

H_reactants = -250.31

def H_products(T_AD):
    return (
          8 * get_property('CO2', T_AD, 'delta_hf')
        + 9 * get_property('H2O', T_AD, 'delta_hf')
        + 47 * get_property('N2', T_AD, 'delta_h_sensible')
        + 8 * get_property('CO2', T_AD, 'delta_h_sensible')
        + 9 * get_property('H2O', T_AD, 'delta_h_sensible')
          )

def f(T_AD):
    return H_products(T_AD) - H_reactants

# Use Brent's method to solve for T_AD...

root = brentq(f, 2000, 3000)

print("Adiabatic Flame Temperature @ phi = 1:", root,"K")


# Lean case...

# Energy balance...

H_reactants = -250.31

def H_products(T_AD):
    return (
          8 * get_property('CO2', T_AD, 'delta_hf')
        + 9 * get_property('H2O', T_AD, 'delta_hf')
        + 47 * get_property('N2', T_AD, 'delta_h_sensible')
        + 8 * get_property('CO2', T_AD, 'delta_h_sensible')
        + 9 * get_property('H2O', T_AD, 'delta_h_sensible')
        + 12.5 * get_property('O2', T_AD, 'delta_h_sensible')
    )

def f(T_AD):
    return H_products(T_AD) - H_reactants

# Use Brent's method to solve for T_AD...

root = brentq(f, 2000, 3000)

print("Adiabatic Flame Temperature @ phi = 0.5:", root,"K")

