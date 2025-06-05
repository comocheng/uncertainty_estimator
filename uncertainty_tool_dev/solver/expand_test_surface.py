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
# simplest example available

# +
import rmgpy.chemkin
import os

import cantera as ct
import numpy as np
import logging

import rmgpy.solver  # import SimpleReactor, TerminationTime, SurfaceReactor
import rmgpy.quantity  #import Quantity
import rmgpy.rmg.listener # import SimulationProfileWriter, SimulationProfilePlotter
import rmgpy.rmg.settings  #import ModelSettings, SimulatorSettings

import matplotlib
import matplotlib.pyplot as plt
# %matplotlib inline

# logging.basicConfig(level=logging.DEBUG)
import rmgpy.constants as constants
from rmgpy.kinetics import SurfaceArrhenius, StickingCoefficient
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
# from rmgpy.solver.surface import SurfaceReactor
from rmgpy.species import Species
from rmgpy.thermo import ThermoData, NASA, NASAPolynomial


# +
sevy_mech = '/home/moon/nitridation/fe110_20241206/'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_chemkin_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')

species_list_gas, reaction_list_gas = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict)
species_list_surf, reaction_list_surf = rmgpy.chemkin.load_chemkin_file(s_chemkin_surface, s_dict)

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_chemkin_surface)

output_directory = 'surface_csvs'

os.makedirs(output_directory, exist_ok=True)
# -

len(reaction_list_surf)

len(reaction_list_gas)

reaction_list_surf[1].kinetics

surf.reactions()[1].rate.pre_exponential_factor

surf.reactions()[1].rate.temperature_exponent

gas.reactions()

species_list[9]

# +
core_species = species_list
edge_species = []
core_reactions = reaction_list
edge_reactions = []

T = 800.0
# T = 1000.0
P_initial = 1.0e5

# NH3 mech
initial_mole_fractions = {
    species_list[4]: 0.5,  # O2
    species_list[3]: 0.5,  # NH3
}
initial_surface_coverages = {
    species_list[9]: 1.0
}


rxn_system = rmgpy.solver.SurfaceReactor(
    T,
    P_initial,
    n_sims=1,
    initial_gas_mole_fractions=initial_mole_fractions,
    initial_surface_coverages=initial_surface_coverages,
    surface_volume_ratio=(1.0, "m^-1"),
    surface_site_density=(2.72e-9, "mol/cm^2"),
    termination=[],
)
# in chemkin, the sites are mostly occupied in about 1e-8 seconds.

# +
# rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions)


rxn_system.initialize_model(core_species, core_reactions, edge_species, edge_reactions,
                            surface_species=[x for x in species_list if x.contains_surface_site()],
                            surface_reactions=[x for x in reaction_list if x.is_surface_reaction()]
                           )
tlist = np.logspace(-15, -6, 501, dtype=float)

# -

# Integrate to get the solution at each time point
t = []
y = []
ys_rmg = []
reaction_rates = []
species_rates = []
ys_rmg.append(rxn_system.core_species_concentrations)
t.append(rxn_system.t)
# You must make a copy of y because it is overwritten by DASSL at
# each call to advance()
y.append(rxn_system.y.copy())
reaction_rates.append(rxn_system.core_reaction_rates.copy())
species_rates.append(rxn_system.core_species_rates.copy())
print("time: ", t)
print("moles:", y)
print("reaction rates:", reaction_rates)
print("species rates:", species_rates)
for t1 in tlist:
    rxn_system.advance(t1)
    t.append(rxn_system.t)
    # You must make a copy of y because it is overwritten by DASSL at
    # each call to advance()
    y.append(rxn_system.y.copy())
    reaction_rates.append(rxn_system.core_reaction_rates.copy())
    species_rates.append(rxn_system.core_species_rates.copy())
    ys_rmg.append(rxn_system.core_species_concentrations)
# # Convert the solution vectors to np arrays
# t = np.array(t, float)
# y = np.array(y, float)
# reaction_rates = np.array(reaction_rates, float)
# species_rates = np.array(species_rates, float)
# V = constants.R * rxn_system.T.value_si * np.sum(y) / rxn_system.P_initial.value_si



# +

# plot results
for i in range(len(y[0])):
    if core_species[i].contains_surface_site():
        continue
    y_hs = [x[i] for x in y]
    plt.plot(t, y_hs, label=str(core_species[i]))
plt.xlabel('time (s)')
plt.ylabel('concentration')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()

for i in range(len(y[0])):
    if not core_species[i].contains_surface_site():
        continue
#     y_hs = [x[i] for x in y]
    y_hs = [x[i] for x in ys_rmg]
    plt.plot(t, y_hs, label=str(core_species[i]))
plt.xlabel('time (s)')
plt.ylabel('concentration')
# plt.yscale('log')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()

# -







# +
mech_yaml = '/home/moon/nitridation/fe110_20241206/chem_annotated.yaml'


surf = ct.Interface(mech_yaml, 'surface1')
gas = surf.adjacent['gas']

T = 800.0
# T = 1000.0
P = 1e5

# nitridation
initial_mole_fractions = 'NH3(2):0.5, O2(3):0.5'
initial_surface_coverages = 'X(1):1.0'

# initial_mole_fractions = '[CH3]: 1.0'  # test data
# initial_surface_coverages = '[Pt]:1.0'

gas.TPX = T, P, initial_mole_fractions
surf.TP = T, P
surf.coverages = initial_surface_coverages

catalyst_weight = 4.24e-3
cat_site_per_wt = 5*61.67*1e-6*1e3 # [mol/kg] 1e-6mol/micromole, 1000g/kg
site_density = (
    surf.site_density * 1000
)  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
cat_area = (catalyst_weight * cat_site_per_wt) / site_density  # [m^3]

gas_reactor = ct.IdealGasReactor(gas, energy='off')
# gas_reactor = ct.IdealGasReactor(gas)
gas_reactor.volume = 1.0

# gas_reactor = ct.IdealGasConstPressureReactor(gas, energy='off')
# opting not to set catalyst area here... we'll see how that goes


surface_reactor = ct.ReactorSurface(surf, gas_reactor)

net = ct.ReactorNet([gas_reactor])



# #########


# upstream = ct.Reservoir(gas, name='upstream')

# # # create a reservoir for the reactor to exhaust into. The composition of
# # # this reservoir is irrelevant.
# downstream = ct.Reservoir(gas, name='downstream')

# V = 30 # velocity m/s (setting very high values allows for effectively fixed gas-phase composition)

# mass_flow_rate = V * gas.density * 1.0 * 1.0
# # # The mass flow rate into the reactor will be fixed by using a
# # # MassFlowController object.
# m = ct.MassFlowController(upstream, gas_reactor, mdot=mass_flow_rate)

# v = ct.PressureController(gas_reactor, downstream, master=m, K=1e-5)
# ##############

history = []
gas_history = []

ts_ct = [net.time]
ys_ct = [gas.X]
ys_surf = [surf.coverages]
y_conc = [surf.concentrations]
while net.time < 1e-6:
    net.step()
    ts_ct.append(net.time)
    ys_ct.append(gas.X)
    ys_surf.append(surf.coverages)
    y_conc.append(surf.concentrations)
# plot results
for i in range(len(ys_ct[0])):
    y_hs = [x[i] for x in ys_ct]
    plt.plot(ts_ct, y_hs, label=gas.species_names[i])
plt.xlabel('time (s)')
plt.ylabel('concentration')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()

for i in range(len(ys_surf[0])):
    y_hs = [x[i] for x in ys_surf]
    plt.plot(ts_ct, y_hs, label=surf.species_names[i])
plt.xlabel('time (s)')
plt.ylabel('concentration')
# plt.yscale('log')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()


# +
for i in range(len(y[0])):
    if not core_species[i].contains_surface_site():
        continue
    y_hs = [x[i] for x in ys_rmg]
    plt.plot(t, y_hs, label=str(core_species[i]))
plt.xlabel('time (s)')
plt.ylabel('concentration')
# plt.yscale('log')
plt.legend(bbox_to_anchor=(1.1, 1.05))


for i in range(len(ys_surf[0])):
    y_hs = [x[i] *1e3 for x in y_conc]
    plt.plot(ts_ct, y_hs, label=surf.species_names[i], linestyle='dashed')
plt.xlabel('time (s)')
plt.ylabel('concentration')

plt.legend(bbox_to_anchor=(1.1, 1.05))
# plt.show()

# plt.yscale('log')
# plt.ylim([1e-18, 1e-3])
# plt.show()

# +
# ys_rmg
# -

y_conc

tlist


