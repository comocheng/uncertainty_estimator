# script to export mechanism uncertainty data to npy files

# Instructions:
# 1. change the database settings to match the input file
# 2. specify the mechanism path
# 3. make sure RMG-database is on the same commit as the one used to run RMG
# 4. make sure RMG-Py is on the correlated_uncertainty branch


import os
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty


# # Load the RMG Model

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
uncertainty.assign_parameter_uncertainties()

np.save(os.path.join(mech_dir, 'gao_reaction_uncertainty.npy'), uncertainty.kinetic_input_uncertainties)
np.save(os.path.join(mech_dir, 'gao_species_uncertainty.npy'), uncertainty.thermo_input_uncertainties)
