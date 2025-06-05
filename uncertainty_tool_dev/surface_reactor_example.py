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
# trying to figure out how to use the SurfaceReactor

# +
import rmgpy.chemkin
import os

import logging

import rmgpy.solver  # import SimpleReactor, TerminationTime, SurfaceReactor
import rmgpy.quantity  #import Quantity
import rmgpy.rmg.listener # import SimulationProfileWriter, SimulationProfilePlotter
import rmgpy.rmg.settings  #import ModelSettings, SimulatorSettings
# T = Quantity(T)
# P = Quantity(P)

logging.basicConfig(level=logging.DEBUG)

# from rmgpy.rmg.listener import SimulationProfileWriter, SimulationProfilePlotter
# from rmgpy.rmg.settings import ModelSettings, SimulatorSettings

# +
# sevy_mech = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin'
# s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
# s_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
# s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
# species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)
# -

sevy_mech = '/home/moon/rmg/RMG-Py/test/rmgpy/test_data/chemkin/chemkin_py/surface/'
s_chemkin = os.path.join(sevy_mech, 'chem-gas.inp')
s_surface = os.path.join(sevy_mech, 'chem-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)

# +
# test data
initial_mole_fractions = {
    species_list[0]: 0.8,  # Ar
    species_list[5]: 0.1,  # O2
    species_list[4]: 0.1,  # CH4
}
initial_surface_coverages = {
    species_list[6]: 1.0
}

# my mech
# initial_mole_fractions = {
#     species_list[0]: 0.8,  # Ar
#     species_list[4]: 0.1,  # O2
#     species_list[3]: 0.1,  # CH4
# }
# initial_surface_coverages = {
#     species_list[22]: 1.0
# }

T = (1000, 'K')
P = (1, 'bar')
termination_time = (1.0, 'ms')
termination = [rmgpy.solver.TerminationTime(rmgpy.quantity.Quantity(termination_time))]

surface_volume_ratio = (1600, 'm^-1')  # Emily says she made this up. I respect that
surface_site_density = (2.7200E-09, 'mol/cm^2')


reaction_system = rmgpy.solver.SurfaceReactor(
    T,
    P,
    initial_mole_fractions,
    initial_surface_coverages,
    surface_volume_ratio,
    surface_site_density,
    n_sims=1,
    termination=termination,
    sensitive_species=[species_list[3]]
)


# -

reaction_system.initialize_model(
    species_list,
    reaction_list,
    [],
    [],
    surface_species=[x for x in species_list if x.contains_surface_site()],
    surface_reactions=[x for x in reaction_list if x.is_surface_reaction()],
)




reaction_system.t

for i in range(10):
    print(f'Step {i},\tt={reaction_system.t}')
    result = reaction_system.step(1e-3)
    print(f'\t{result}')

reaction_system.t

reaction_system.advance(1e-12)

reaction_system.t

reaction_system.advance(2e-12)

reaction_system.t

reaction_system.advance(1e-6)

# +
output_directory = 'test_surf'
os.makedirs(output_directory, exist_ok=True)
sens_worksheet = []
reaction_system_index = 0
# reaction_system.attach(SimulationProfileWriter(
#     output_directory, reaction_system_index, species_list))
# reaction_system.attach(SimulationProfilePlotter(
#     output_directory, reaction_system_index, species_list))

simulator_settings = rmgpy.rmg.settings.SimulatorSettings()  # defaults

model_settings = rmgpy.rmg.settings.ModelSettings()  # defaults
model_settings.tol_move_to_core = 0.1
model_settings.tol_interrupt_simulation = 1.0
model_settings.tol_keep_in_edge = 0.0

set_results = reaction_system.simulate(
    core_species=species_list,
    core_reactions=reaction_list,
    edge_species=[],
    edge_reactions=[],
    surface_species=[x for x in species_list if x.contains_surface_site()],
    surface_reactions=[x for x in reaction_list if x.is_surface_reaction()],
    model_settings=model_settings,
    simulator_settings=simulator_settings,
    sensitivity=False,
#     sens_worksheet=sens_worksheet,
)
# -







reaction_system.step(0.0012)





reaction_system.t



terminated, False, invalid_objects, surface_species, surface_reactions, self.t, conversion


5112 / 36



[species_list[3]]

species_list[22]

species_list[3]


