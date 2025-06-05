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
import rmgpy.data.rmg
import rmgpy
import rmgpy.chemkin

rmgpy.settings["test_data.directory"]



a = None

k = a if a is not None else rmgpy.settings["test_data.directory"]

k



database = rmgpy.data.rmg.RMGDatabase()
database.load(
    path=os.path.join(rmgpy.settings["test_data.directory"], "testing_database"),
    thermo_libraries=["primaryThermoLibrary"],
    reaction_libraries=["GRI-Mech3.0"],
    kinetics_families=[
        "R_Recombination",
        "Disproportionation",
        "R_Addition_MultipleBond",
        "H_Abstraction",
        "intra_H_migration",
    ],
    testing=True,
    depository=False,
    solvation=False,
    surface=False,
)

# +
# Load an example chemkin file or two
chemkin1 = '/home/moon/autoscience/fuels/butane_20240501/chem_annotated.inp'
sp_dict1 = os.path.join(os.path.dirname(chemkin1), 'species_dictionary.txt')
species_list1, reaction_list1 = rmgpy.chemkin.load_chemkin_file(chemkin1, sp_dict1)


chemkin2 = '/home/moon/nitridation/fe110_shifted_O2_nobetavdw/chem_annotated-gas.inp'
sp_dict2 = os.path.join(os.path.dirname(chemkin2), 'species_dictionary.txt')
species_list2, reaction_list2 = rmgpy.chemkin.load_chemkin_file(chemkin2, sp_dict2)
# -

count = 0
for i, r in enumerate(reaction_list1 + reaction_list2):
    if not hasattr(r, 'family'):
        continue
    if r.family not in database.kinetics.families:
        continue
    
    src = database.kinetics.families[r.family].extract_source_from_comments(r)
    print(i, src)
    if count >= 10:
        break
    count += 1
    

database.kinetics.families

database.kinetics.families['Disproportionation']




