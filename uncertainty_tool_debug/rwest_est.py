# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import os
import sys
import pickle
import copy
import numpy as np
import rmgpy.chemkin
import rmgpy
import rmgpy.tools.uncertainty
import rmgpy.kinetics.uncertainties

import random

import rmgpy.kinetics
import matplotlib.pyplot as plt
# %matplotlib inline

import scipy
import scipy.misc
from scipy.misc import derivative

sys.path.append('/home/moon/autoscience/reaction_calculator/database')
import database_fun
# -







# +
## SEVY LOOK HERE
x = 0.5
def f(x, a, b):
    return  x / (x + np.exp(-a) + np.exp(-b))


def f(x, a, b):
    return  x / np.minimum(a, b)


def f(x, a, b):
    return  x / (0.5 * np.float_power(a, 2.0) + 0.5 * np.float_power(b, 2.0))


a, da = 1.3, 0.7
b, db = 1.2, 1.0
#######

dist_a = scipy.stats.norm(a, da) # assume the da is a standard deviation
dist_b = scipy.stats.norm(b, db)

def estimate(N):
    samples_a = dist_a.rvs(N)
    samples_b = dist_b.rvs(N)
    y = f(x, samples_a, samples_b)
    print(f"mean: {y.mean()}, std: {y.std()}")
    return y.mean(), y.std()

Ns = np.logspace(1, 6, 15).astype(int)
means = []
stds = []
for N in Ns:
    mean, std = estimate(N)
    means.append(mean)
    stds.append(std)
plt.semilogx(Ns, means)
plt.ylabel('mean')
plt.show()
plt.plot(Ns, stds)
plt.ylabel('std')
plt.semilogx()

# +
df_da = derivative(lambda p: f(x, p, b), a)
df_db = derivative(lambda p: f(x, a, p), b)
y = f(x, a, b)
dy = np.sqrt((df_da * da)**2 + (df_db * db)**2)
dy # analytical

print(f"analytical: {y}, std: {dy}")
print(f"numerical: {means[-1]}, std: {stds[-1]}")

# -


