import math as m
from scipy.optimize import brentq

def f(alpha):
    return 0.330 - (alpha / (1 - alpha)) * m.sqrt((0.5 * alpha) / (0.5 * alpha + 1))

root = brentq(f, 0+1e-6, 1-1e-6)

print("Dissociation Fraction (alpha)", root)