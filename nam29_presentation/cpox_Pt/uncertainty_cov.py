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
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty

import matplotlib.pyplot as plt
# %matplotlib inline
# -

mech_dir = '/home/moon/uncertainty_estimator/nam29_presentation/cpox_Pt'
chemkin = os.path.join(mech_dir, 'chem_annotated-gas.inp')
surface = os.path.join(mech_dir, 'chem_annotated-surface.inp')
species_dict = os.path.join(mech_dir, 'species_dictionary.txt')

uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory=os.path.join(mech_dir, 'rmg_uncertainty'))

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin, species_dict, surface_path=surface)
uncertainty.species_list = species_list
uncertainty.reaction_list = reaction_list
# uncertainty.load_model(chemkin, species_dict, surface_path=surface)

# --------------- CAUTION!!! Databases here must match the ones used to generate the mechanism
thermo_libs = [
    'surfaceThermoPt111',
    'primaryThermoLibrary',
    'thermo_DFT_CCSDTF12_BAC',
    'DFT_QCI_thermo'
]
kinetic_libs = [
    'Surface/CPOX_Pt/Deutschmann2006_adjusted',
    'BurkeH2O2inArHe'
]

kinetics_families = [
    'default',
    'Surface_Adsorption_Single',
    'Surface_Adsorption_vdW',
    'Surface_Adsorption_Dissociative',
    'Surface_Dissociation',
    'Surface_Abstraction',
    'Surface_Dissociation_Double_vdW',
    'Surface_Dissociation_vdW',
    'Surface_Abstraction_vdW',
    'Surface_Dissociation_Beta',
    'Surface_Adsorption_Bidentate',
    'Surface_Bidentate_Dissociation',
    'Surface_Monodentate_to_Bidentate',
    'Surface_Dissociation_to_Bidentate', 
    'Surface_vdW_to_Bidentate',
    'Surface_Adsorption_Dissociative_Double',
    'Surface_Abstraction_Beta',
    # 'Surface_Abstraction_Beta_double_vdW',
    'Surface_Dissociation_Double',
    'Surface_Dissociation_Beta_vdW',
    'Surface_Abstraction_Beta_vdW',
    'Surface_Abstraction_Single_vdW',
]

uncertainty.load_database(
    thermo_libraries=thermo_libs,
    kinetics_families=kinetics_families,
    reaction_libraries=kinetic_libs,
)

# Get the different kinetic and thermo sources
uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties(correlated=True)

uncertainty.get_thermo_covariance_matrix()

uncertainty.get_kinetic_covariance_matrix()



plt.matshow(uncertainty.get_overall_covariance_matrix())
plt.colorbar()
plt.clim([0, 0.1])

U = uncertainty.get_overall_covariance_matrix()

len(uncertainty.species_list)

len(uncertainty.reaction_list)

for i in range(len(uncertainty.species_list)):
    for j in range(len(uncertainty.species_list), 250):
        if U[i,j] > 0:
            print(i, j)

for i in range(len(uncertainty.reaction_list)):
    print(type(uncertainty.reaction_list[i].kinetics))

uncertainty.reaction_sources_dict[uncertainty.reaction_list[152]]['Rate Rules'][1]['training']

for i in range(uncertainty.reaction_sources_dict[uncertainty.reaction_list[152]]['Rate Rules'][1]['training'][0][1].data)

uncertainty.reaction_sources_dict[uncertainty.reaction_list[152]]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[152]]['Rate Rules'][1]['training'][0][0].data

uncertainty.reaction_sources_dict[uncertainty.reaction_list[152]]['Rate Rules'][1]['training'][0][1].data


