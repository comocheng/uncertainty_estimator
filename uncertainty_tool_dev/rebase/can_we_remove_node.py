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
import re
import numpy as np
import rmgpy.tools.uncertainty

import matplotlib.pyplot as plt

# %matplotlib inline
# -

# pick an annotated chemkin file to analyze
chemkin_file = '/home/moon/uncertainty_estimator/uncertainty_tool_dev/ethane_limit_families/chemkin/chem_annotated.inp'
dict_file = '/home/moon/uncertainty_estimator/uncertainty_tool_dev/ethane_limit_families/chemkin/species_dictionary.txt'


# make an uncertainty object
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='parse_comment_rebase')
uncertainty.load_model(chemkin_file, dict_file)

# +
# Load the database
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

# -

uncertainty.extract_sources_from_model()
uncertainty.compile_all_sources()

uncertainty.assign_parameter_uncertainties(correlated=False)

for i in range(len(uncertainty.species_list)):
    print(uncertainty.species_list[i], uncertainty.thermo_input_uncertainties[i])

for i in range(len(uncertainty.reaction_list)):
    family = getattr(uncertainty.reaction_list[i], 'family', 'PDEP')
    print(i, uncertainty.reaction_list[i], family, uncertainty.kinetic_input_uncertainties[i])

uncertainty.reaction_sources_dict[uncertainty.reaction_list[29]]



print(uncertainty.reaction_list[19].kinetics)

uncertainty.reaction_sources_dict[uncertainty.reaction_list[22]]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[22]]['Rate Rules'][1]['template']

source_dict = uncertainty.reaction_sources_dict[uncertainty.reaction_list[22]]['Rate Rules'][1]

exact = source_dict['exact']
rule_weights = [ruleTuple[-1] for ruleTuple in source_dict['rules']]
training_weights = [trainingTuple[-1] for trainingTuple in source_dict['training']]

node = uncertainty.reaction_sources_dict[uncertainty.reaction_list[22]]['Rate Rules'][1]['template'][0]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[22]]['Rate Rules'][1]['rules'][0][0].data

uncertainty.reaction_sources_dict[uncertainty.reaction_list[50]]['Rate Rules'][1]['rules'][0][0].data.comment

np.log(np.float_power(10.0, 10.0))

np.log10(np.exp(23.23))

rootnode = source_dict['template'][0].parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent

rootnode.data_count

source_dict['rules'][0]

node.label

std_dev_regex = r'Total Standard Deviation in ln\(k\): (\d+\.\d+)?'
n_train_regex = r'BM rule fitted to (\d+)? training reactions'

mytext = uncertainty.database.kinetics.families['Disproportionation'].rules.entries[node.label][0].data.comment

mytext

std_dev_match = re.search(std_dev_regex, mytext)
std_dev = float(std_dev_match[1])

n_train_match = re.search(n_train_regex, mytext)
n_train = int(n_train_match[1])

std_dev

n_train




