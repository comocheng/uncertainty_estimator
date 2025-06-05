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
import rmgpy.chemkin
import os

import numpy as np

import logging

import rmgpy.util
import rmgpy.solver  # import SimpleReactor, TerminationTime, SurfaceReactor
import rmgpy.quantity  #import Quantity
import rmgpy.rmg.listener # import SimulationProfileWriter, SimulationProfilePlotter
import rmgpy.rmg.settings  #import ModelSettings, SimulatorSettings

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

# +
surface = True

if surface:
    sevy_mech = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin'
    s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
    s_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
    s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)
    output_directory = 'surface_csvs'
else:
    sevy_mech = '/home/moon/rmg/RMG-Py/examples/rmg/superminimal/chemkin'
    s_chemkin = os.path.join(sevy_mech, 'chem_annotated.inp')
    s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict)
    reaction_list = [reaction_list[0], reaction_list[1], reaction_list[95], reaction_list[96]]
    output_directory = 'simple_csvs'

os.makedirs(output_directory, exist_ok=True)
# -



# +
T = (1000, 'K')
P = (1, 'bar')
termination_time = (0.5, 'ms')
termination = [rmgpy.solver.TerminationTime(rmgpy.quantity.Quantity(termination_time))]
sensitivity_threshold=1e-3

if surface:
    initial_gas_mole_fractions = {
        species_list[0]: 0.131246 / 0.207992,  # Ar
        species_list[4]: 0.03488 / 0.207992,  # O2
        species_list[3]: 1.0 - (0.131246 / 0.207992) - 0.03488 / 0.207992,  # CH4
    }
    initial_surface_coverages = {
        species_list[22]: 1.0
    }

    surface_volume_ratio = (1600, 'm^-1')  # Emily says she made this up. I respect that
    surface_site_density = (0.7200E-09, 'mol/cm^2')
    sensitive_species = [species_list[3]]
    
else:
    initial_mole_fractions={
        species_list[4]: 0.67,  # H2
        species_list[5]: 0.33,  # O2
    }
    sensitive_species = [species_list[4]]
# -



if surface:
    reaction_system = rmgpy.solver.SurfaceReactor(
        T,
        P,
        initial_gas_mole_fractions,
        initial_surface_coverages,
        surface_volume_ratio,
        surface_site_density,
        n_sims=1,  # default of simple reactor
        termination=termination,
        sensitive_species=sensitive_species,
        sensitivity_threshold=sensitivity_threshold,
    )
else:
    reaction_system = rmgpy.solver.SimpleReactor(
        T=T,
        P=P,
        initial_mole_fractions=initial_mole_fractions,
        termination=termination,
        sensitive_species=sensitive_species,
        sensitivity_threshold=sensitivity_threshold
    )

# +
# Create the csv worksheets for logging sensitivity
rmgpy.util.make_output_subdirectory(output_directory, 'solver')
sens_worksheet = []
reaction_system_index = 0
for spec in reaction_system.sensitive_species:
    csvfile_path = os.path.join(output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1, spec.index))
    sens_worksheet.append(csvfile_path)

reaction_system.attach(rmgpy.rmg.listener.SimulationProfileWriter(
    output_directory, reaction_system_index, species_list))
reaction_system.attach(rmgpy.rmg.listener.SimulationProfilePlotter(
    output_directory, reaction_system_index, species_list))

simulator_settings = rmgpy.rmg.settings.SimulatorSettings()  # defaults

model_settings = rmgpy.rmg.settings.ModelSettings()  # defaults
model_settings.tol_move_to_core = 0.1
model_settings.tol_interrupt_simulation = 1.0
model_settings.tol_keep_in_edge = 0.0
# -

reaction_system.simulate(
    core_species=species_list,
    core_reactions=reaction_list,
    edge_species=[],
    edge_reactions=[],
    surface_species=[x for x in species_list if x.contains_surface_site()],
    surface_reactions=[x for x in reaction_list if x.is_surface_reaction()],
    model_settings=model_settings,
    simulator_settings=simulator_settings,
    sensitivity=True,
    sens_worksheet=sens_worksheet
)

np.zeros(3) - np.ones(4)

dir(reaction_list[47].kinetics)


