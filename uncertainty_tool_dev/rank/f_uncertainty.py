# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline

f_values = np.array([0.265, 0.286, 0.181, 0.205, 0.600, 0.370, 0.403, 0.427,
                     0.555, 0.541, 0.502, 0.262, 0.568, 0.400, 0.400, 0.244,
                     0.476, 0.700, 0.500, 0.300, 0.300, 0.500, 0.700])

plt.hist(f_values, 32)

multipliers = np.float_power(10.0, f_values)

lnks_uniform = np.log(multipliers) / np.sqrt(3.0)

lnks = np.log(multipliers) / 1.96

plt.hist(lnks, alpha=0.5)
plt.hist(lnks_uniform, alpha=0.5)

# +
f = 1.0

print(f'f={f}')

multiplier = np.float_power(10.0, f)
print(f'kmax/k0={multiplier}')

lnk_uniform = np.log(multiplier) / np.sqrt(3.0)
print(f'uniform delta ln k = {lnk_uniform}')

normal_95 = np.log(multiplier) / 1.96
print(f'normal, 95% cutoff ln k = {normal_95}')

normal_3sigma = f * np.log(10.0) / 3
print(f'normal, 3sigma cutoff ln k = {normal_3sigma}')
# -

np.log10(1.1)

np.log10(1.4)

np.log10(1.5)

1.0 * np.log(10) / 3.0

np.log10(1.9)

# +
barrier_uncertainty = 3.1 / 2.0 * 4184 # J/mol


# barrier_uncertainty = 13020 # J/mol


# barrier_uncertainty = 12 * 4184

T_ref = 1000.0

R = 8.314


# k = A e ^ - Ea / RT

np.log10(np.exp( barrier_uncertainty / (R * T_ref)))

# -

# # $$k_0 = A e^{-\frac{E_a}{RT}}$$
#
# Here we use $2\sigma$ for the maximum barrier energy
# # $$k_\text{max} = A e^{-\frac{E_a - \Delta E_a}{RT}}$$
#
# # $$k_\text{max} = e^{\frac{\Delta E_a}{RT}} A e^{-\frac{E_a}{RT}}$$
#
# So the multiplicative factor is $e^{\frac{\Delta E_a}{RT}}$

# +
print(np.exp(barrier_uncertainty / (R * T_ref)))  # multiplier

print(np.log10(np.exp(barrier_uncertainty / (R * T_ref)))) # f

# +
# BACs



# -

# # $$k_0 = A e^{-\frac{E_a}{RT}}$$
#
# Here we use $2\sigma$ for the maximum barrier energy
# # $$k_\text{max} = A e^{-\frac{E_a - 2\sigma_{ E_a}}{RT}}$$
#
# # $$\frac{k_\text{max}}{k_0} = e^{\frac{2\sigma_{ E_a}}{RT}}$$
#
# So the multiplicative factor is $e^{\frac{2 \sigma_{E_a}}{RT}}$





# +
# Baulch

logks = np.array([0.2, 0.2, 0.3, 0.3, 0.15, 0.3, 0.35, 0.25, 0.5, 0.4, 0.25, 0.3, 0.5, 0.25, 0.1, 0.4, 0.3, 0.2, 0.5, 0.2, 0.4, 0.5, 0.25, 0.5, 0.2, 0.4, 0.225, 0.25, 0.65, 0.3, 0.2, 0.2, 0.6, 0.3, 0.3, 0.3, 0.3, 0.5, 0.3, 0.2, 0.6, 0.4, 0.75, 0.55, 0.6, 0.8333, 0.3667, 0.5, 0.15, 0.3167, 0.2, 0.4, 0.3, 0.2, 0.3, 0.3, 0.6, 0.65, 0.2167, 0.3, 0.3, 0.3, 0.2, 0.2, 0.45, 0.45, 0.3, 0.4, 0.55, 0.35, 0.2, 0.4, 0.25, 0.4, 0.3, 0.4, 0.5, 0.15, 0.4, 0.2, 0.3, 0.5, 0.2, 0.2, 0.3, 0.25, 0.4, 0.3, 0.2, 0.35, 0.275, 0.15, 0.3, 1.0, 0.385, 0.3, 0.55, 0.5, 0.2, 0.25, 0.385, 0.55, 0.5, 0.625, 0.55, 0.35, 0.2667, 0.3, 0.4, 0.35, 0.25, 0.3667, 0.475, 0.075, 0.25, 0.325, 0.15, 0.3, 0.25, 0.2, 0.4, 0.15, 0.125, 0.225, 0.3667, 0.3, 0.75, 0.5, 0.2, 0.3, 0.225, 0.2, 1.0, 0.3, 0.5, 0.5, 0.5, 0.5, 0.7, 0.7, 0.52, 0.55, 0.5, 0.5, 0.5, 0.4, 0.3, 0.2, 0.25, 0.2, 0.25, 0.2, 0.2, 0.3, 0.225, 0.15, 0.3, 0.3, 0.2, 0.3, 0.6, 0.5, 0.3, 0.25, 0.3, 0.2, 0.5, 0.3, 0.2, 0.3, 0.4, 0.4, 0.55, 0.3, 0.3, 0.35, 0.25, 0.3, 0.2, 0.2, 0.3, 0.5, 0.25, 0.2, 0.3, 0.3, 0.4, 0.25, 0.2, 0.3, 0.2, 0.3, 0.3, 0.3, 0.5, 0.3, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.7, 0.7, 0.5, 0.55, 0.4, 0.6, 0.4, 0.5, 0.3, 0.35, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.3, 0.325, 0.3667, 0.3, 0.3, 0.2, 0.3, 0.4, 0.3, 0.35, 0.5, 0.4, 0.4, 0.3, 0.3, 0.35, 0.5, 0.5, 0.55, 0.7, 0.7, 0.5, 0.5, 0.325, 0.35])
# -

logks

# +
multipliers = np.exp( np.sqrt(3.0) * logks)

# multipliers = np.exp(3 * logks)
f = np.log10(multipliers)
# -

np.median(f)

plt.hist(f)

u95s = np.array([
    0.0119502868068834,
    0.0358508604206501,
    0.12189292543021,
    0.152963671128107,
    0.1,
    0.3,
    0.3,
    0.2,
    0.669216061185468,
    0.669216061185468,
    0.669216061185468,
    2.19885277246654,
    3.11185468451243,
    1.13766730040153,
    1.23087954110899,
    2.075,
    1.19502868068834,
    0.956022944550669,
    17.5,
    18.125,
    10,
    8.5,
    6.25,
    21,
    7.425,
    15.5,
    12.25,
    3.975,
    18.325
])
ranks = np.array([
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    4,
    5,
    5,
    5,
    5,
    5,
    5,
    6,
    6,
    6,
    6,
    6,
    7,
    7,
    7,
    7,
    7,
    8
])

# # Thermo justification plot

# +
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

cs_thermo = [colors[i - 1] for i in ranks]
ax = plt.gca()
ax.spines[ 'top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.scatter(ranks, u95s, color=cs_thermo)
plt.yscale('log')
plt.xlabel('RMG Rank', fontsize=16)
plt.ylabel('2$\sigma$ (kcal/mol)', fontsize=16)


# plt.axhline(y=0.1, color=colors[0], label='Rank=1', zorder=0, linewidth=0.7)
# plt.axhline(y=0.2, color=colors[1], label='Rank=2', zorder=0, linewidth=0.7)
# plt.axhline(y=0.5, color=colors[2], label='Rank=3', zorder=0, linewidth=0.7)
# plt.axhline(y=1.5, color=colors[3], label='Rank=4', zorder=0, linewidth=0.7)
# plt.axhline(y=3.1, color=colors[4], label='Rank=5', zorder=0, linewidth=0.7)
# plt.axhline(y=8.0, color=colors[5], label='Rank=6', zorder=0, linewidth=0.7)
# plt.axhline(y=12.0, color=colors[6], label='Rank=7', zorder=0, linewidth=0.7)
# plt.axhline(y=18.0, color=colors[7], label='Rank=8', zorder=0, linewidth=0.7)
xlims = plt.xlim()
xticks = plt.xticks()

plt.hlines(xmin=0.5, xmax=1.5, y=0.1, color=colors[0], label='Rank=1', zorder=0, linewidth=0.7)
plt.hlines(xmin=1.5, xmax=2.5, y=0.2, color=colors[1], label='Rank=2', zorder=0, linewidth=0.7)
plt.hlines(xmin=2.5, xmax=3.5, y=0.5, color=colors[2], label='Rank=3', zorder=0, linewidth=0.7)
plt.hlines(xmin=3.5, xmax=4.5, y=1.5, color=colors[3], label='Rank=4', zorder=0, linewidth=0.7)
plt.hlines(xmin=4.5, xmax=5.5, y=3.1, color=colors[4], label='Rank=5', zorder=0, linewidth=0.7)
plt.hlines(xmin=5.5, xmax=6.5, y=8.0, color=colors[5], label='Rank=6', zorder=0, linewidth=0.7)
plt.hlines(xmin=6.5, xmax=7.5, y=12.0, color=colors[6], label='Rank=7', zorder=0, linewidth=0.7)
plt.hlines(xmin=7.5, xmax=8.5, y=18.0, color=colors[7], label='Rank=8', zorder=0, linewidth=0.7)

step_color = 'gray'

plt.vlines(x=1.5, ymin=0.1, ymax=0.2, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=2.5, ymin=0.2, ymax=0.5, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=3.5, ymin=0.5, ymax=1.5, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=4.5, ymin=1.5, ymax=3.1, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=5.5, ymin=3.1, ymax=8.0, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=6.5, ymin=8.0, ymax=12.0, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=7.5, ymin=12.0, ymax=18.0, color=step_color, zorder=0, linewidth=0.7)


plt.xlim(xlims)
plt.xticks(xticks[0], xticks[1])


plt.legend()
# -

len(f_ranks)

10 ** 0.05

# # Kinetics justification plot

# +
fs = np.array([0.0791812460476248, 0.27, 0.301029995663981, 0.0791812460476248, 0.113943352306837, 0.176091259055681, 0.477121254719662, 0.278753600952829, 0.431363764158987, 0.278753600952829, 0.431363764158987, 0.204119982655925, 0.278753600952829, 0.437115254338034, 0.218557627169017, 0.698970004336019, 0.437115254338034, 1.46609174689932, 0.698970004336019, 0.301029995663981, 1, 1, 1])
f_ranks = np.array([1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9])

T_ref = 1000.0
f_from_u95s = np.log10(np.exp(u95s * 4184 / (8.314 * T_ref)))

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

cs = [colors[i - 1] for i in f_ranks]
ax = plt.gca()
ax.spines[ 'top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.scatter(f_ranks, fs, color=cs)

plt.scatter(ranks, f_from_u95s, marker='*', color=cs_thermo)

plt.yscale('log')
plt.xlabel('RMG Rank', fontsize=16)
plt.ylabel('$f$', fontsize=16)

xlims = plt.xlim()
xticks = plt.xticks()

plt.hlines(xmin=0.5, xmax=1.5, y=0.05, color=colors[0], label='Rank=1', zorder=0, linewidth=0.7)
plt.hlines(xmin=1.5, xmax=2.5, y=0.1, color=colors[1], label='Rank=2', zorder=0, linewidth=0.7)
plt.hlines(xmin=2.5, xmax=3.5, y=0.2, color=colors[2], label='Rank=3', zorder=0, linewidth=0.7)
plt.hlines(xmin=3.5, xmax=4.5, y=0.4, color=colors[3], label='Rank=4', zorder=0, linewidth=0.7)
plt.hlines(xmin=4.5, xmax=5.5, y=0.5, color=colors[4], label='Rank=5', zorder=0, linewidth=0.7)
plt.hlines(xmin=5.5, xmax=6.5, y=0.7, color=colors[5], label='Rank=6', zorder=0, linewidth=0.7)
plt.hlines(xmin=6.5, xmax=7.5, y=1.0, color=colors[6], label='Rank=7', zorder=0, linewidth=0.7)
plt.hlines(xmin=7.5, xmax=8.5, y=1.5, color=colors[7], label='Rank=8', zorder=0, linewidth=0.7)
plt.hlines(xmin=8.5, xmax=9.5, y=2.0, color=colors[8], label='Rank=9', zorder=0, linewidth=0.7)
step_color = 'gray'

plt.vlines(x=1.5, ymin=0.05, ymax=0.1, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=2.5, ymin=0.1, ymax=0.2, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=3.5, ymin=0.2, ymax=0.4, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=4.5, ymin=0.4, ymax=0.5, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=5.5, ymin=0.5, ymax=0.7, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=6.5, ymin=0.7, ymax=1.0, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=7.5, ymin=1.0, ymax=1.5, color=step_color, zorder=0, linewidth=0.7)
plt.vlines(x=8.5, ymin=1.5, ymax=2.0, color=step_color, zorder=0, linewidth=0.7)


plt.xlim(xlims)
plt.xticks(xticks[0], xticks[1])
plt.legend()

# -


