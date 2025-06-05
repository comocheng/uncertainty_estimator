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

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
from scipy.optimize import fsolve





# +
# Read in the mechanism
model_file = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/cantera/chem_annotated.yaml'
# model_file = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241112/mincat/cantera/chem_annotated.yaml'
model_file = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241112/og_lib/cantera/chem_annotated.yaml'



gas = ct.Solution(model_file, "gas")
surf = ct.Interface(model_file, "surface1", [gas])
# -

gas.species()

# +
# gas.species()
# -



# +
# Define reactor conditions
CO_0 = 1.0  # kmol / m^3
O2_0 = 1.0  # kmol / m^3
Ar_0 = 1.0  # kmol / m^3
REACTOR_VOLUME = 1.0  # m^3
REACTOR_TEMPERATURE = 1500  # K
REACTOR_PRESSURE = 100000.0  # 1 bar = 100000 Pa
MAX_SIMULATION_TIME = 1.0
CONCENTRATIONS = {
    'CH4(2)': 0.1,
    'O2(3)': 0.1,
}
CONCENTRATIONS['Ar'] = 1.0 - CONCENTRATIONS['CH4(2)'] - CONCENTRATIONS['O2(3)']


x_CH4 = CONCENTRATIONS['CH4(2)']
x_O2 = CONCENTRATIONS['O2(3)']
x_Ar = CONCENTRATIONS['Ar']


# initialize T and P
gas.TPX = REACTOR_TEMPERATURE, REACTOR_PRESSURE, CONCENTRATIONS
surf.TP = REACTOR_TEMPERATURE, REACTOR_PRESSURE

volume_flow = 1.0

# Catalyst settings (really need to double check these)
catalyst_weight = 4.24e-3
cat_site_per_wt = 5*61.67*1e-6*1e3 # [mol/kg] 1e-6mol/micromole, 1000g/kg
site_density = (
    surf.site_density * 1000
)  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
cat_area = (catalyst_weight * cat_site_per_wt) / site_density  # [m^3]
surf.coverages = "X(1):1.0"

# +
gas_reactor = ct.IdealGasReactor(gas)
gas_reactor.volume = REACTOR_VOLUME
surface_reactor = ct.ReactorSurface(surf, gas_reactor, A=cat_area)

# set up mass flow controllers
inlet = ct.Reservoir(gas)
exhaust = ct.Reservoir(gas)
# -

# N = 10000
# delta_t = 1e-2
t_end = 1e-6
t = np.linspace(0, t_end, 1001)
delta_t = t[1] - t[0]
# t = np.arange(0, N * delta_t, delta_t)
X_cov = np.zeros(len(t))
CO_cov = np.zeros(len(t))
O_cov = np.zeros(len(t))
CO2_cov = np.zeros(len(t))
gas_history = []
surf_history = []

for i in range(0, len(t)):
    X_cov[i] = surf.coverages[surf.species_names.index('X(1)')]
    CO_cov[i] = surf.coverages[surf.species_names.index('COX(23)')]
    O_cov[i] = surf.coverages[surf.species_names.index('OX(25)')]
    CO2_cov[i] = surf.coverages[surf.species_names.index('CO2X(22)')]
    gas_history.append(gas.X)
    surf_history.append(surf.coverages)
    surf.advance_coverages(delta_t)

gas.species_names.index('CH4(2)')



# +
CH4_conc = [x[gas.species_names.index('CH4(2)')] for x in gas_history]
O2_conc = [x[gas.species_names.index('O2(3)')] for x in gas_history]

OX_cov = [x[surf.species_names.index('OX(25)')] for x in surf_history]
# -

plt.plot(t, CH4_conc, label='CH4')
plt.plot(t, O2_conc, label='O2')
plt.legend()

plt.plot(t, OX_cov, label='OX')
         

# i = 8
for i in range(0, 10):
    plt.plot(t, [x[i] for x in surf_history], label=surf.species_names[i])
plt.legend(loc='lower right')
plt.yscale('log')
plt.ylabel('coverage')
plt.show()


