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
# script to explain the uncertainty tool

# +
import os
import rmgpy.chemkin
import rmgpy.tools.uncertainty

import importlib
importlib.reload(rmgpy.tools.uncertainty)
# -





# pick an annotated chemkin file to analyze
chemkin_file = '/home/moon/uncertainty_estimator/uncertainty_tool_dev/ethane_limit_families/chemkin/chem_annotated.inp'
dict_file = '/home/moon/uncertainty_estimator/uncertainty_tool_dev/ethane_limit_families/chemkin/species_dictionary.txt'


# make an uncertainty object
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations')


# +
# load the reactions and species (directly from the file)
uncertainty.load_model(chemkin_file, dict_file)

# or indirectly by passing the species/reaction lists through memory
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(
    species_list=species_list,
    reaction_list=reaction_list,
    output_directory='uncertainty_calculations'
)


# +
# Load the database - this takes a while because the averaging up apparently needs Julia
# make sure all of the libraries/families match what was used to generate the mechanism
thermo_libraries = ['primaryThermoLibrary', 'BurkeH2O2']
reaction_libraries = ['BurkeH2O2inN2']
kinetics_families = [
    'Disproportionation',
    'H_Abstraction',
    'intra_H_migration',
    'R_Recombination',
    'Intra_Disproportionation',
]

uncertainty.load_database(
    thermo_libraries=thermo_libraries,
    kinetics_families=kinetics_families,
    reaction_libraries=reaction_libraries,
)


# +
# Extract sources from model fills out reaction sources dict and the species thermo sources dict
uncertainty.extract_sources_from_model()

# can be accessed like
print(uncertainty.species_sources_dict[uncertainty.species_list[0]])
print()
print(uncertainty.reaction_sources_dict[uncertainty.reaction_list[0]])
# uncertainty.compile_all_sources appears to be completely useless, might want to remove it from RMG

# -

# ## Kinetics Source Dictionary

for i in range(len(uncertainty.reaction_list)):
#     display
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    print(i, source.keys())
    if i > 10:
        break

# See what's in Rate Rules ones
for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Rate Rules' in source:
        
        print(i, source['Rate Rules'])
        

# See what's in the training ones
for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Training' in source:
        
        print(i, source['Training'])
        









# compile sources helps summarize where all the kinetics/thermo came from
uncertainty.compile_all_sources()
uncertainty.all_kinetic_sources['Rate Rules']['H_Abstraction']

# +
# Next, assign the parameter uncertainties
# Now you have to choose whether to include correlations or not
# We'll start without correlations

uncertainty.assign_parameter_uncertainties()

print(uncertainty.thermo_input_uncertainties[0])
print()

print(uncertainty.kinetic_input_uncertainties[0])
# -













uncertainty.database.thermo.libraries

uncertainty

uncertainty.all_kinetic_sources

uncertainty.assign_parameter_uncertainties()





uncertainty.kinetic_input_uncertainties

uncertainty.all_kinetic_sources

uncertainty.reaction_sources_dict[uncertainty.reaction_list[0]]

uncertainty.species_sources_dict


