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
# I should use the mech diff tool I built earlier...

# +
import rmgpy.chemkin
import numpy as np
import os

import matplotlib.pyplot as plt
# %matplotlib inline

import rmgpy.tools.diffmodels

import rmgpy.rmg.model  # import ReactionModel
import rmgpy.rmg.output   #import save_diff_html
# -



# +
emily_mech = './cpox_rh_emily'
e_chemkin = os.path.join(emily_mech, 'chem_annotated-gas.inp')
e_surface = os.path.join(emily_mech, 'chem_annotated-surface.inp')
e_dict = os.path.join(emily_mech, 'species_dictionary.txt')
e_sp, e_rxn = rmgpy.chemkin.load_chemkin_file(e_chemkin, e_dict, surface_path=e_surface, read_comments=False)

sevy_mech = './cpox_rh_20241028'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
s_sp, s_rxn = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)
# -





common_species, unique_species1, unique_species2, common_reactions, unique_reactions1, unique_reactions2 = rmgpy.tools.diffmodels.execute(
    os.path.join(sevy_mech, 'chem_annotated.inp'),
    os.path.join(sevy_mech, 'species_dictionary.txt'),
    '',
    os.path.join(sevy_mech, 'chem_annotated.inp'),
    os.path.join(sevy_mech, 'species_dictionary.txt'),
    '',
    read_comments=False
)

len(unique_species1)

# +
model1 = rmgpy.rmg.model.ReactionModel()
model2 = rmgpy.rmg.model.ReactionModel()

emily_mech = './cpox_rh_emily'
e_chemkin = os.path.join(emily_mech, 'chem_annotated-gas.inp')
e_surface = os.path.join(emily_mech, 'chem_annotated-surface.inp')
e_dict = os.path.join(emily_mech, 'species_dictionary.txt')
e_sp, e_rxn = rmgpy.chemkin.load_chemkin_file(e_chemkin, e_dict, surface_path=e_surface, read_comments=False)

sevy_mech = './cpox_rh_20241028'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
s_sp, s_rxn = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)


model1.species, model1.reactions = rmgpy.chemkin.load_chemkin_file(
    e_chemkin, e_dict, read_comments=False,
    surface_path=e_surface)
model2.species, model2.reactions = rmgpy.chemkin.load_chemkin_file(
    s_chemkin, s_dict, read_comments=True,
    surface_path=s_surface)

common_reactions, unique_reactions1, unique_reactions2 = rmgpy.tools.diffmodels.compare_model_reactions(model1, model2)
common_species, unique_species1, unique_species2 = rmgpy.tools.diffmodels.compare_model_species(model1, model2)

outputDir = '.'
output_path = outputDir + 'diff.html'
rmgpy.rmg.output.save_diff_html(output_path, common_species, unique_species1, unique_species2, common_reactions, unique_reactions1,
                   unique_reactions2)

# -

len(model1.species)

len(model2.species)

unique_reactions2

e_surface


