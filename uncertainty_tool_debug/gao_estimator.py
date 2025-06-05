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

import os
import sys
import pickle
import rmgpy.chemkin
import rmgpy.tools.uncertainty

# +
# Load the model

# Must use annotated chemkin file
chemkin_file = 'ethane_model/chem_annotated.inp'
dict_file = 'ethane_model/species_dictionary.txt'

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file)


# +
# Run Gao estimation of input parameters (takes a long time to load database)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations')
uncertainty.load_model(chemkin_file, dict_file)

uncertainty.load_database(
    thermo_libraries=['DFT_QCI_thermo', 'primaryThermoLibrary'],
    kinetics_families='default',
    reaction_libraries=[],
)
uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()
# -

with open('rxn_sources.pickle', 'wb') as f:
    pickle.dump(uncertainty.reaction_sources_dict, f)
with open('rxn_uncertainties.pickle', 'wb') as f:
    pickle.dump(uncertainty.kinetic_input_uncertainties, f)

for i in range(0, 50):
    if 'Rate Rules' in uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]].keys():
        print(i, uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]])
        print()



# +
# Go through and summarize kinetic sources?

for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Library' in source:
        pass
    if 'PDep' in source:
        pass
    if 'Training' in source:
        pass # this is usually (always?) an exact match
#         print(source)
    if 'Rate Rules' in source:
#         if source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['rules']:
#             print(source)
        if source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['template']:
            print(source)
#         print(source)
# -

source['Rate Rules']



uncertainty.kinetic_input_uncertainties

uncertainty.reaction_list[-2].kinetics

uncertainty.reaction_sources_dict[uncertainty.reaction_list[-2]]

uncertainty.kinetic_input_uncertainties[-2]

uncertainty.reaction_list[-2]



i = -2
node_name = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']
uncertainty.database.kinetics.families[uncertainty.reaction_list[i].family].rules.entries[node_name][0].long_desc
#                     std_dev_matches = re.search(r'Standard Deviation in ln\(k\): ([0-9]*.[0-9]*)', long_desc)
#                     std_dev = -1.0
#                     if std_dev_matches is not None:
#                         std_dev = float(std_dev_matches[1])


