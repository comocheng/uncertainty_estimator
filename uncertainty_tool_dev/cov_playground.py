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
# a notebook where i can work on developing the code that will export the uncertainty covariance maxtrices

import os
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty
import matplotlib.pyplot as plt

# %matplotlib inline

import importlib
importlib.reload(rmgpy.tools.uncertainty)
importlib.reload(rmgpy.data.thermo)


# -

# # Simple surface

# pick an annotated chemkin file to analyze
chemkin_file = '/home/moon/uncertainty_estimator/cpox_pt/toy_surface2/chemkin/chem_annotated-gas.inp'
chemkin_file_surface = '/home/moon/uncertainty_estimator/cpox_pt/toy_surface2/chemkin/chem_annotated-surface.inp'
dict_file = '/home/moon/uncertainty_estimator/cpox_pt/toy_surface2/chemkin/species_dictionary.txt'


# +
# make an uncertainty object
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_surf1')
# load the reactions and species (directly from the file)
uncertainty.load_model(chemkin_file, dict_file)

# or indirectly by passing the species/reaction lists through memory
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file, surface_path=chemkin_file_surface)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(
    species_list=species_list,
    reaction_list=reaction_list,
    output_directory='uncertainty_calculations'
)


# +
# Load the database - this takes a while because the averaging up apparently needs Julia
# make sure all of the libraries/families match what was used to generate the mechanism
thermo_libraries = ['surfaceThermoPt111_reduced', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo']
reaction_libraries = ['Surface/CPOX_Pt/Deutschmann2006_adjusted','BurkeH2O2inArHe']
kinetics_families = [
    'default',
    'Surface_Abstraction',
    'Surface_Abstraction_Single_vdW',
    'Surface_Abstraction_vdW',
    #'Surface_Addition_Single_vdW',
    #'Surface_Adsorption_Abstraction_vdW',
    'Surface_Adsorption_Bidentate',
    'Surface_Adsorption_Dissociative',
    'Surface_Adsorption_Dissociative_Double',
    'Surface_Adsorption_Double',
    'Surface_Adsorption_Single',
    'Surface_Adsorption_vdW',
    'Surface_Bidentate_Dissociation',
    'Surface_Dissociation',
    'Surface_Dissociation_Beta',
    'Surface_Dissociation_Double_vdW',
    'Surface_Dissociation_vdW',
    'Surface_Dissociation_to_Bidentate',
    #'Surface_DoubleBond_to_Bidentate',
    #'Surface_Dual_Adsorption_vdW',
    'Surface_EleyRideal_Addition_Multiple_Bond',
    'Surface_Migration',
    'Surface_vdW_to_Bidentate'
]

uncertainty.load_database(
    thermo_libraries=thermo_libraries,
    kinetics_families=kinetics_families,
    reaction_libraries=reaction_libraries,
)

# -







uncertainty.extract_sources_from_model()
uncertainty.compile_all_sources()


# +
uncertainty.assign_parameter_uncertainties(correlated=False)

for i in range(len(uncertainty.species_list)):
    print(uncertainty.species_list[i], uncertainty.thermo_input_uncertainties[i])



# +
plt.matshow()


# -

uncertainty.assign_parameter_uncertainties(correlated=True)
uncertainty.get_thermo_covariance_matrix()
uncertainty.get_kinetic_covariance_matrix()

plt.matshow(uncertainty.thermo_covariance_matrix)
plt.colorbar()
# plt.clim([0, 0.1])

# +
# show uncertainty without surface additions:

# make an uncertainty object
uncertainty0 = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_surf0')
# load the reactions and species (directly from the file)
uncertainty0.load_model(chemkin_file, dict_file)

# or indirectly by passing the species/reaction lists through memory
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file, surface_path=chemkin_file_surface)
uncertainty0 = rmgpy.tools.uncertainty.Uncertainty(
    species_list=species_list,
    reaction_list=reaction_list,
    output_directory='uncertainty_calculations'
)


# +
# Load the database - this takes a while because the averaging up apparently needs Julia
# make sure all of the libraries/families match what was used to generate the mechanism
thermo_libraries = ['surfaceThermoPt111_reduced', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo']
reaction_libraries = ['Surface/CPOX_Pt/Deutschmann2006_adjusted','BurkeH2O2inArHe']
kinetics_families = [
    'default',
    'Surface_Abstraction',
    'Surface_Abstraction_Single_vdW',
    'Surface_Abstraction_vdW',
    #'Surface_Addition_Single_vdW',
    #'Surface_Adsorption_Abstraction_vdW',
    'Surface_Adsorption_Bidentate',
    'Surface_Adsorption_Dissociative',
    'Surface_Adsorption_Dissociative_Double',
    'Surface_Adsorption_Double',
    'Surface_Adsorption_Single',
    'Surface_Adsorption_vdW',
    'Surface_Bidentate_Dissociation',
    'Surface_Dissociation',
    'Surface_Dissociation_Beta',
    'Surface_Dissociation_Double_vdW',
    'Surface_Dissociation_vdW',
    'Surface_Dissociation_to_Bidentate',
    #'Surface_DoubleBond_to_Bidentate',
    #'Surface_Dual_Adsorption_vdW',
    'Surface_EleyRideal_Addition_Multiple_Bond',
    'Surface_Migration',
    'Surface_vdW_to_Bidentate'
]

uncertainty0.load_database(
    thermo_libraries=thermo_libraries,
    kinetics_families=kinetics_families,
    reaction_libraries=reaction_libraries,
)

# -

uncertainty0.extract_sources_from_model()
uncertainty0.compile_all_sources()

uncertainty0.thermo_input_uncertainties[-1]

# +

g_un = rmgpy.tools.uncertainty.ThermoParameterUncertainty(
    dG_surf_lib=1.5,
    dG_ADS=0,
    dG_ADS_group=0
)
# uncertainty0.assign_parameter_uncertainties(correlated=True)
uncertainty0.assign_parameter_uncertainties(correlated=True, g_param_engine=g_un)
uncertainty0.get_thermo_covariance_matrix()
uncertainty0.get_kinetic_covariance_matrix()
# -

# %matplotlib inline
plt.matshow(uncertainty0.thermo_covariance_matrix)
plt.colorbar()
plt.clim([0, 0.1])



uncertainty.species_list[20]

for i in range(100):
    if uncertainty.reaction_list[i].is_surface_reaction():
        print(i)
        break







# +
families = [r.family for r in uncertainty.reaction_list]

surface_rule_counts = []
gas_rule_counts = []

# for family in uncertainty.database.kinetics.families:
for family in families:
    try:
        
        if family.startswith('Surface'):
            surface_rule_counts.append(len(uncertainty.database.kinetics.families[family].rules.entries))
        else:
            gas_rule_counts.append(len(uncertainty.database.kinetics.families[family].rules.entries))
    except KeyError:
        pass
# -

np.average(np.array(surface_rule_counts))

np.average(np.array(gas_rule_counts))

surface_rule_counts





gas_rule_counts



len(uncertainty.database.kinetics.families['CO_Disproportionation'].rules.entries)

len(uncertainty.database.kinetics.families['Surface_Adsorption_vdW'].rules.entries)

len(uncertainty.database.kinetics.families['H_Abstraction'].rules.entries)





for i in range(0, 100):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Rate Rules' in source:
        print(i, uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]])
        print()
                                            



print(r'\n')

uncertainty.species_list[24].thermo.comment

# +
uncertainty.species_list[24].thermo.comment

comment = uncertainty.species_list[24].thermo.comment
comment = comment.replace('\n', ' ')
comment
# -

assert 'Adsorption correction: + Thermo group additivity estimation:' in comment

uncertainty.database.thermo.extract_source_from_comments(uncertainty.species_list[24])

uncertainty.species_list[24].thermo.comment





uncertainty.thermo_input_uncertainties



# # Gas-phase

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

# compile sources helps summarize where all the kinetics/thermo came from
uncertainty.compile_all_sources()
uncertainty.all_kinetic_sources['Rate Rules']['H_Abstraction']

# +
# Next, assign the parameter uncertainties
# Now you have to choose whether to include correlations or not
# We'll start without correlations

uncertainty.assign_parameter_uncertainties(correlated=True)

print(uncertainty.thermo_input_uncertainties[0])
print()

print(uncertainty.kinetic_input_uncertainties[0])
# -

# # check thermo covariance matrix

uncertainty.assign_parameter_uncertainties(correlated=True)
plt.matshow(uncertainty.get_thermo_covariance_matrix())
plt.colorbar()




uncertainty.assign_parameter_uncertainties(correlated=True)
plt.matshow(uncertainty.get_thermo_covariance_matrix())
plt.colorbar()
plt.clim([0, 0.1])

# +
uncertainty.assign_parameter_uncertainties(correlated=True)
A = uncertainty.get_thermo_covariance_matrix()

uncertainty.assign_parameter_uncertainties(correlated=False)
B = uncertainty.get_thermo_covariance_matrix()

assert np.all(np.isclose(np.diag(np.diagonal(A)), B))
# -

# # check kinetics covariance matrix

for family in uncertainty.all_kinetic_sources['Rate Rules'].keys():
    print(family)

uncertainty.assign_parameter_uncertainties(correlated=True)
plt.matshow(uncertainty.get_kinetic_covariance_matrix())
plt.colorbar()


# +
# %matplotlib inline

plt.matshow(uncertainty.kinetic_covariance_matrix)
plt.colorbar()
plt.clim([0, 0.25])
# -

# uncertainty.assign_parameter_uncertainties(correlated=True)
# A = uncertainty.get_kinetic_covariance_matrix()
plt.matshow(A)
plt.colorbar()
plt.clim([0, 0.25])

# +
uncertainty.assign_parameter_uncertainties(correlated=True)
A = uncertainty.get_kinetic_covariance_matrix()

uncertainty.assign_parameter_uncertainties(correlated=False)
B = uncertainty.get_kinetic_covariance_matrix()

assert np.all(np.isclose(np.diag(np.diagonal(A)), B))

# +
# tell me about all the off diagonals:

uncertainty.assign_parameter_uncertainties(correlated=True)
A = uncertainty.get_kinetic_covariance_matrix()
for i in range(A.shape[0]):
    source_i = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    for j in range(i):
        if A[i, j] != 0:
            try:
                if source_i['Rate Rules'][0] == 'H_Abstraction':
                    print(i, j)

                    source_j = uncertainty.reaction_sources_dict[uncertainty.reaction_list[j]]
                    print(source_i['Rate Rules'][1]['training'])
                    print()
                    print(source_j['Rate Rules'][1]['training'])
                    print()
                    print()
                    print()
            except KeyError:
                pass
# -

uncertainty.reaction_sources_dict[uncertainty.reaction_list[55]]



uncertainty.all_kinetic_sources['Rate Rules'].keys()

uncertainty.database.kinetics.families['H_Abstraction'].auto_generated

rxn_map = uncertainty.database.kinetics.families['Disproportionation'].get_reaction_matches(
            thermo_database=uncertainty.database.thermo,
            remove_degeneracy=True,
            get_reverse=True,
            exact_matches_only=False,
            fix_labels=True)





len(rxn_map['Root'])



uncertainty.reaction_sources_dict[uncertainty.reaction_list[44]]



uncertainty.database.kinetics.families['Disproportionation'].get_training_depository().entries[16]



uncertainty.database.kinetics.families['H_Abstraction'].get_training_depository().entries[16]



uncertainty.database.kinetics.families['H_Abstraction'].rules.entries['C/H3/Cd\H_Cd\H2;C_rad/H2/Cs\H3']

uncertainty.database.kinetics.families['H_Abstraction'].rules.entries['C/H3/Cd\H_Cd\H2;C_rad/H2/Cs\H3'][0].data

# +
for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Rate Rules' in source.keys() and source['Rate Rules'][0] == 'H_Abstraction':
        training_reactions = [t[1] for t in source['Rate Rules'][1]['training']]
#         print(training_reactions)
        for t in training_reactions:
            for j in range(len(uncertainty.reaction_list)):
                if uncertainty.reaction_list[j].is_isomorphic(t.item):
                    print(i, j, source)
                    print()
            
#         print(source)
#         print()
# -

uncertainty.reaction_sources_dict[uncertainty.reaction_list[18]]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[42]]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[44]]





uncertainty.reaction_sources_dict[uncertainty.reaction_list[16]]

uncertainty.kinetic_input_uncertainties[16]

uncertainty.kinetic_input_uncertainties[18]



uncertainty.reaction_sources_dict[uncertainty.reaction_list[16]]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[18]]

uncertainty.kinetic_input_uncertainties[18]

uncertainty.reaction_list[14]





source['Rate Rules']

uncertainty.kinetic_input_uncertainties[11]



uncertainty.database.kinetics.families['H_Abstraction'].rules.entries['C3H6 + C2H5 <=> C2H6 + C3H5']

uncertainty.k

for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Rate Rules' in source.keys():
        if source['Rate Rules'][0] == 'H_Abstraction':
            print(source['Rate Rules'])
#             for entry in source['Rate Rules'][1]['rules']:
#                 print(entry)
#             for train in source['Rate Rules'][1]['training']:
#                 print(train)
            
            print()



np.max(np.abs(A - B))

(A - B).where([np.abs(A - B) > 18])





for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        if A[i,j] != B[i,j]:
            print(i, j, A[i,j], B[i, j])



A[-1,-1]

B[-1, -1]



uncertainty.kinetic_input_uncertainties













uncertainty.thermo_input_uncertainties[-1]

uncertainty.thermo_input_uncertainties[-2]

uncertainty.get_thermo_covariance_matrix()[-2, -1]

.1*.2

uncertainty.thermo_input_uncertainties

A = np.zeros(3)

if A:
    print('yes')



uncertainty.thermo_input_uncertainties

uncertainty.species_sources_dict[uncertainty.species_list[2]]



np.diag(uncertainty.thermo_input_uncertainties)

uncertainty.get_thermo_covariance_matrix()

plt.matshow()


