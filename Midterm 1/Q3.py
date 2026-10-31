from scipy.optimize import brentq

# Energy balance...

LHV = 1460.8e3 # J/mol

def H_products(T_AD):
    return (
        + 16.732 * (20.79 * (T_AD - 298))
        + 2.35 * (37.13 * (T_AD - 298))
        + 3.36 * (33.59 * (T_AD - 298))
        + 0.893 * (29.38 * (T_AD - 298))
          )

def f(T_AD):
    return H_products(T_AD) - LHV

# Use Brent's method to solve for T_AD...

root = brentq(f, 1000, 4000)

print("Adiabatic Flame Temperature @ phi = 0.8:", root,"K")