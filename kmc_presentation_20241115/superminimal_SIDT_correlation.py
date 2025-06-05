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
# script to get the correlations noted between SIDT nodes based on training reactions

# +
import os
import sys
import pickle
import copy
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty
import rmgpy.kinetics.uncertainties

import rmgpy.tools.canteramodel
import random

import rmgpy.kinetics
import matplotlib.pyplot as plt
# %matplotlib inline

# +
# Load the model

# Must use annotated chemkin file
chemkin_file = '/home/moon/uncertainty_estimator/kmc_presentation_20241115/ethane/chemkin/chem_annotated.inp'
dict_file = '/home/moon/uncertainty_estimator/kmc_presentation_20241115/ethane/chemkin/species_dictionary.txt'

# Run Gao estimation of input parameters (takes a long time to load database)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations')
uncertainty.load_model(chemkin_file, dict_file)

# +
# # Superminimal
# thermo_libs = ['primaryThermoLibrary']
# kinetics_families = [
#     'H_Abstraction',
#     'Disproportionation',
#     'R_Recombination',
#     'Birad_recombination',
#     'Birad_R_Recombination'
# ]

# Ethane + reaction libraries
thermo_libs = ['primaryThermoLibrary', 'BurkeH2O2']
reaction_libraries = ['BurkeH2O2inN2']
kinetics_families = 'default'

uncertainty.load_database(
    thermo_libraries=thermo_libs,
    kinetics_families=kinetics_families,
    reaction_libraries=reaction_libraries,
)


# -

uncertainty.species_list[5].thermo

uncertainty.reaction_list[0]

global g_param_engine
g_param_engine = rmgpy.tools.uncertainty.ThermoParameterUncertainty(
    dG_library=np.sqrt(1.5),
    dG_QM=np.sqrt(3.0),
    dG_GAV=np.sqrt(1.5),
    dG_group=np.sqrt(0.10)
)

uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties(g_param_engine=g_param_engine, correlated=True)

uncorrelated_Gs = uncertainty.thermo_input_uncertainties

variances = np.float_power(uncorrelated_Gs, 2.0)

uncertainty.thermo_input_uncertainties[2]

uncertainty.species_list[2]

uncertainty.species_sources_dict[uncertainty.species_list[2]]

np.diag(cov_G)

# +
# Run sensitivity

ethane = rmgpy.species.Species().from_smiles('CC')
C2H4 = rmgpy.species.Species().from_smiles('C=C')
mapping = rmgpy.tools.canteramodel.get_rmg_species_from_user_species([ethane, C2H4], uncertainty.species_list)

initial_mole_fractions = {mapping[ethane]: 1.0}
T = (1300, 'K')
P = (1, 'bar')
termination_time = (0.5, 'ms')
sensitive_species = [mapping[ethane], mapping[C2H4]]

uncertainty.sensitivity_analysis(
    initial_mole_fractions,
    sensitive_species,
    T,
    P,
    termination_time,
    sensitivity_threshold=1e-6
)

# -

output_results = uncertainty.local_analysis(
    sensitive_species,
    reaction_system_index=0,
    correlated=False,
)

output = {}
for sens_species in [sensitive_species[1]]:
    csvfile_path = os.path.join(uncertainty.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(1,
                                                                     sens_species.index))
    time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)
    
    # Assign uncertainties
    thermo_data_list = []
    reaction_data_list = []
    for data in data_list:
        if data.species:
            for species in uncertainty.species_list:
                if species.to_chemkin() == data.species:
                    index = uncertainty.species_list.index(species)
                    break
            else:
                raise Exception('Chemkin name {} of species in the CSV file does not match anything in the '
                                'species list.'.format(data.species))

            data.uncertainty = uncertainty.thermo_input_uncertainties[index]
            thermo_data_list.append(data)

        if data.reaction:
            rxn_index = int(data.index) - 1
            data.uncertainty = uncertainty.kinetic_input_uncertainties[rxn_index]
            reaction_data_list.append(data)
    

sens = [x.data[-1] for x in reaction_data_list]

sens



reaction_data_list[0].data[-1]

time.data[-1]

reaction_data_list[0].uncertainty

uncertainty.reaction_list[8]

uncertainty.kinetic_input_uncertainties

uncertainty.species_sources_dict[uncertainty.species_list[10]]['GAV']

g_param_engine.



# +
# try out the covariance of thermo data

# g_param_engine = rmgpy.tools.uncertainty.ThermoParameterUncertainty()

cov_G = np.zeros((len(uncertainty.species_list), len(uncertainty.species_list)))

# for i, species in enumerate(self.species_list):
#     for j, other_species in enumerate(self.species_list):
#         if i == j:
#             cov_G[i, j] = self.thermo_input_uncertainties[i] ** 2
#         else:
#             cov_G[i, j] = 0.0
for i in range(len(uncertainty.species_list)):
    for j in range(len(uncertainty.species_list)):
#         cov_G[i, j] = 
        
        for label_i, dpG_i in uncertainty.thermo_input_uncertainties[i].items():
            for label_j, dpG_j in uncertainty.thermo_input_uncertainties[j].items():
                if label_i != label_j:
                    continue  # not correlated unless the labels match
                if 'estimation' in label_i.lower():
                    # assume estimation errors are uncorrelated unless it's an identity - same species
                    if i == j:
                        cov_G[i, j] += dpG_i * dpG_j
                elif 'group' in label_i.lower():
                    degeneracy_i = dpG_i / g_param_engine.dG_group
                    degeneracy_j = dpG_j / g_param_engine.dG_group
                    cov_G[i, j] += degeneracy_i * degeneracy_j * (g_param_engine.dG_group ** 2.0)
                elif 'library' in label_i.lower():
                    degeneracy_i = dpG_i / g_param_engine.dG_library
                    degeneracy_j = dpG_j / g_param_engine.dG_library
                    cov_G[i, j] += degeneracy_i * degeneracy_j * (g_param_engine.dG_library ** 2.0)
                else:
                    raise NotImplementedError()

# -



cov_G

np.sqrt(0.1)

plt.matshow(cov_G)
plt.colorbar()

np.diag(cov_G)

np.float_power(uncorrelated_Gs, 2.0)

uncertainty.species_list[6]

uncertainty.thermo_input_uncertainties

len(uncertainty.reaction_list)

len(data_list)



data_list[2].label

time.data

len(time.data)





output[sensitive_species[0]][1]



# +
# plot the thermo covariance matrix?

for i in range(len(uncertainty.species_list)):
    print(i, uncertainty.thermo_input_uncertainties[i])
# -

uncertainty.species_sources_dict[uncertainty.species_list[16]]

uncertainty.species_list[16]

len(uncertainty.species_list)

for i in range(len(uncertainty.reaction_list)):
    print(i, uncertainty.kinetic_input_uncertainties[i])

uncertainty.reaction_sources_dict[uncertainty.reaction_list[59]]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[18]]


