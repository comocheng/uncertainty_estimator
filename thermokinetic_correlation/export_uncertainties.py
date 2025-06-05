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
# A script to export the reactions from an RMG mechanism and the database
# produces 3 files:

# make a dict entry for index

# species_dict_file = 'species_dict.pickle'
# reaction_dict_file = 'reaction_dict.pickle'

# correlated_uncertainty_file = 'correlated_uncertainty.pickle'

# -

# example script to show how to unpack information in pickled files
import pickle
import numpy as np

# +
import os
import copy
import itertools
import numpy as np
import scipy.stats
import rmgpy.data.thermo
import rmgpy.data.rmg
import rmgpy.chemkin
import rmgpy.exceptions

import rmgpy.tools.uncertainty
import matplotlib.pyplot as plt
# %matplotlib inline


import importlib
importlib.reload(rmgpy.tools.uncertainty)

# +
# load example 

gas_chemkin = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241020/chem_annotated-gas.inp'
gas_surface = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241020/chem_annotated-surface.inp'
sp_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241020/species_dictionary.txt'
test_species_list, test_reaction_list = rmgpy.chemkin.load_chemkin_file(gas_chemkin, sp_dict, surface_path=gas_surface)
# -

type(test_reaction_list[-1])

# Load the database
database = rmgpy.data.rmg.RMGDatabase()
database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = ['surfaceThermoPt111', 'primaryThermoLibrary'],
    reaction_libraries = ['Surface/CPOX_Pt/Deutschmann2006_adjusted'],
#     kinetics_families = ['Surface_Abstraction'],
    kinetics_families = ['surface'],
    kinetics_depositories = ['training'],
    depository = True,
)


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


# ## Manually add other reactions

# +
display(test_reaction_list[50])

manual_additions = [test_reaction_list[50]]

# -

test_reaction_list[50].products[1].thermo

# +
# Collect all of the training reactions

reaction_list = []
species_list = []  # build this up as you go to make sure all reactions use the same species object
for family in database.kinetics.families.keys():
    training_depository = database.kinetics.families[family].get_training_depository()
    for i, entry in enumerate(training_depository.entries):
        my_reaction = rmgpy.data.kinetics.family.TemplateReaction()
        my_reaction.reactants = training_depository.entries[entry].item.reactants
        my_reaction.products = training_depository.entries[entry].item.products
        my_reaction.family = family

        try:
            template_labels = database.kinetics.families[family].get_reaction_template_labels(my_reaction)
            template = database.kinetics.families[family].retrieve_template(template_labels)
            kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=my_reaction.degeneracy)[0]
            my_reaction.kinetics = kinetics
        except:
            continue
        reaction_list.append(my_reaction)
    
# add the other miscellaneous LIBRARY reactions
for i in range(len(manual_additions)):
    manual_additions[i].kinetics.comment += 'Library ' + manual_additions[i].family
    reaction_list.append(manual_additions[i])
print(len(reaction_list))

# -

'234'.endswith('4')

ord('9')


def increment_label(old_label):
    if len(old_label) >= 2:
        if '_' in old_label:
            tokens = old_label.split('_')
            try:
                old_num = int(tokens[-1])
                return old_label[:-len(tokens[-1])] + str(old_num + 1)
            except ValueError:
                pass
        elif '-' in old_label:
            tokens = old_label.split('-')
            try:
                old_num = int(tokens[-1])
                return old_label[:-len(tokens[-1])] + str(old_num + 1)
            except ValueError:
                pass
            
    return old_label + '_2'


increment_label('A-9')



# +
species_list = []


for i in range(len(reaction_list)):
     for j in range(len(reaction_list[i].reactants)):
        i_thing = get_i_thing(reaction_list[i].reactants[j], species_list)
        if i_thing < 0:
            reaction_list[i].reactants[j].thermo = database.thermo.get_thermo_data(reaction_list[i].reactants[j])
            new_sp = reaction_list[i].reactants[j]
            while new_sp.label in [sp.label for sp in species_list]:
                new_sp.label = increment_label(new_sp.label)
            species_list.append(new_sp)
        else:
            reaction_list[i].reactants[j] = species_list[i_thing]  # reuse species objects

for i in range(len(reaction_list)):
    for j in range(len(reaction_list[i].products)):
        i_thing = get_i_thing(reaction_list[i].products[j], species_list)
        if i_thing < 0:
            reaction_list[i].products[j].thermo = database.thermo.get_thermo_data(reaction_list[i].products[j])
            new_sp = reaction_list[i].products[j]
            while new_sp.label in [sp.label for sp in species_list]:
                new_sp.label = increment_label(new_sp.label)
            species_list.append(new_sp)
        else:
            reaction_list[i].products[j] = species_list[i_thing]  # reuse species objects
        
# -

len(reaction_list)



reaction_list[-1]

len(species_list)

len(species_list)

for i in range(len(species_list)):
    for j in range(i):
        if species_list[i].is_isomorphic(species_list[j]):
            print(i, j)
        if species_list[i].label == species_list[j].label:
            print(i, j)

species_list[39]

species_list[19]

# +
# for i in range(len(reaction_list)):
#     print(reaction_list[i].is_surface_reaction(), reaction_list[i])
# -



# +
# Trim the set to just surface reactions and the species reacting in them
# Per Joy's request, make sure to only include species involved in reactions
surface_only = True

reacting_species_set = set()
for i in range(len(reaction_list)):
    if surface_only and not reaction_list[i].is_surface_reaction():
        continue
    for sp in reaction_list[i].reactants + reaction_list[i].products:
        reacting_species_set.add(species_list.index(sp))
include_species_list = [species_list[i] for i in reacting_species_set]



include_reactions_set = set()
for i in range(len(reaction_list)):
    if surface_only and not reaction_list[i].is_surface_reaction():
        continue
    include_reactions_set.add(i)
include_reactions_list = [reaction_list[i] for i in include_reactions_set]


# -

[sp.label for sp in include_species_list]



len(include_species_list)

for i in include_reactions_list:
    print(i)

# # Check Sticking Coefficient Calculations

include_reactions_list[10].kinetics

# +
Tref = 1000.0
Hrxn = include_reactions_list[10].get_enthalpy_of_reaction(Tref)
print(Hrxn)
my_sticking = include_reactions_list[10].kinetics.to_arrhenius(Hrxn)  # THIS NAME IS WRONG-- it actually returns a StickingCoefficient type of kinetics
print(type(my_sticking))

new_reaction = copy.deepcopy(include_reactions_list[10])
new_reaction.kinetics = my_sticking
# -

include_reactions_list[10].kinetics.get_activation_energy(Hrxn)

include_reactions_list[10].kinetics.get_sticking_coefficient(300, Hrxn)

SDEN = 2.7200E-05  # mol/m^2
new_reaction.get_surface_rate_coefficient(1000, SDEN)



include_species_list[-2]

len(include_species_list)

include_species_list[i].number_of_surface_sites()

include_species_list[-1].number_of_surface_sites()

# # Fill out the species dictionary

# +
species_dict = {}
for i in range(len(include_species_list)):
    print(i, include_species_list[i])
    species_label = include_species_list[i].to_chemkin()
    
    if type(include_species_list[i].thermo) == rmgpy.thermo.thermodata.ThermoData:
        include_species_list[i].thermo = include_species_list[i].thermo.to_nasa(
            np.min(include_species_list[i].thermo.Tdata.value_si),
            np.max(include_species_list[i].thermo.Tdata.value_si),
            np.mean(include_species_list[i].thermo.Tdata.value_si)
        )
    
    species_entry = {
        'index': i,
        'label': species_label,
        'RMG_adjacency_list': include_species_list[i].to_adjacency_list(),
        'NASA0': {
            'coeffs': include_species_list[i].thermo.polynomials[0].coeffs,
            'Tmin': include_species_list[i].thermo.polynomials[0].Tmin.value_si,
            'Tmax': include_species_list[i].thermo.polynomials[0].Tmax.value_si
        },
        'NASA1': {
            'coeffs': include_species_list[i].thermo.polynomials[1].coeffs,
            'Tmin': include_species_list[i].thermo.polynomials[1].Tmin.value_si,
            'Tmax': include_species_list[i].thermo.polynomials[1].Tmax.value_si
        },
        'comment': include_species_list[i].thermo.comment,
        'is_surface_species': include_species_list[i].contains_surface_site(),
        'n_surface_sites': include_species_list[i].number_of_surface_sites(),
        'molecular_weight_kg': include_species_list[i].molecular_weight.value_si,
    }
    
    species_dict[species_label] = species_entry


# -

len(species_dict)



species_dict_file = 'species_dict.pickle'
with open(species_dict_file, 'wb') as f:
    pickle.dump(species_dict, f)

# read back in to check
with open(species_dict_file, 'rb') as f:
    my_data = pickle.load(f)


len(my_data)





# # Fill out the reaction dictionary

reaction_dict = {}
N = len(species_dict)
for i in range(len(include_reactions_list)):
    reaction_label = str(include_reactions_list[i])
    
    
    if type(include_reactions_list[i].kinetics) in [rmgpy.kinetics.surface.SurfaceArrheniusBEP, rmgpy.kinetics.surface.StickingCoefficientBEP]:
        kinetic_entry = {
            'A': include_reactions_list[i].kinetics.A.value,
            'n': include_reactions_list[i].kinetics.n.value,
            'alpha': include_reactions_list[i].kinetics.alpha.value,
            'E0': include_reactions_list[i].kinetics.E0.value,
            'A_units': include_reactions_list[i].kinetics.A.units,
            'E0_units': include_reactions_list[i].kinetics.E0.units
        }
    elif type(include_reactions_list[i].kinetics) in [rmgpy.kinetics.surface.SurfaceArrhenius, rmgpy.kinetics.surface.StickingCoefficient]:
        kinetic_entry = {
            'A': include_reactions_list[i].kinetics.A.value,
            'n': include_reactions_list[i].kinetics.n.value,
            'Ea': include_reactions_list[i].kinetics.Ea.value,
            'A_units': include_reactions_list[i].kinetics.A.units,
            'Ea_units': include_reactions_list[i].kinetics.Ea.units
        }
    
    reaction_entry = {
        'index': N + i,
        'parameterization': str(type(include_reactions_list[i].kinetics)),
        'kinetics': kinetic_entry,
        'comment': include_reactions_list[i].kinetics.comment,
        'reactants': [str(include_reactions_list[i].reactants[j].to_chemkin()) for j in range(len(include_reactions_list[i].reactants))],
        'products': [str(include_reactions_list[i].products[j].to_chemkin()) for j in range(len(include_reactions_list[i].products))],
    }
    
    reaction_dict[reaction_label] = reaction_entry


reaction_dict_file = 'reaction_dict.pickle'
with open(reaction_dict_file, 'wb') as f:
    pickle.dump(reaction_dict, f)

# read back in to check
with open(reaction_dict_file, 'rb') as f:
    my_data = pickle.load(f)


include_reactions_list[10].kinetics

my_data['OX_1 + OX_1 <=> X_3 + X_3 + O2(3)']



rmgpy.chemkin.save_species_dictionary('species_dictionary.txt', species_list)

rmgpy.chemkin.save_chemkin_file('my_reactions.inp', species_list, reaction_list)





# # Run uncertainty Analysis

uncertainty = rmgpy.tools.uncertainty.Uncertainty(species_list=species_list, reaction_list=reaction_list)

# +
# uncertainty = rmgpy.tools.uncertainty.Uncertainty()
# uncertainty.load_model('my_reactions.inp', 'species_dictionary.txt', transport_path=None, surface_path=None)
# -



len(uncertainty.species_list)

len(uncertainty.reaction_list)

uncertainty.database = database

uncertainty.extract_sources_from_model()

uncertainty.assign_parameter_uncertainties(correlated=True)

uncertainty.get_thermo_covariance_matrix()



# Show the thermo covariance matrix
plt.matshow(uncertainty.thermo_covariance_matrix)
plt.colorbar()

# Zoom in to see which parameters are correlated: the answer is not many -- onyl a handful of similar gorups
plt.matshow(uncertainty.thermo_covariance_matrix)
plt.colorbar()
plt.clim([0, 0.005])





# # Now do kinetics covariance

uncertainty.get_kinetic_covariance_matrix()
# Show the kinetics covariance matrix
plt.matshow(uncertainty.kinetic_covariance_matrix)
plt.colorbar()

# +
delta_ln_k = np.sqrt(np.max(uncertainty.kinetic_covariance_matrix))  # that's 1 std dev ln k

# which corresponds to kmax/k of
np.exp(delta_ln_k * np.sqrt(3))

# 6.7 orders of magnitude. Okay. For the worst case scenario, yeah...
np.log10(np.exp(delta_ln_k * np.sqrt(3)))
# -

# zoom in to highlight correlations
plt.matshow(uncertainty.kinetic_covariance_matrix)
plt.colorbar()
plt.clim([0, 0.11])

uncertainty.kinetic_covariance_matrix[:, -1]

uncertainty.overall_covariance_matrix[:, -1]



uncertainty.get_overall_covariance_matrix()
plt.matshow(uncertainty.overall_covariance_matrix)
plt.colorbar()

uncertainty.get_overall_covariance_matrix()
plt.matshow(uncertainty.overall_covariance_matrix)
plt.colorbar()
plt.clim([-0.01, 0.01])

np.save('thermo_cov.npy', uncertainty.thermo_covariance_matrix)
np.save('kinetic_cov.npy', uncertainty.kinetic_covariance_matrix)
np.save('overall_cov.npy', uncertainty.overall_covariance_matrix)



np.sqrt(uncertainty.kinetic_covariance_matrix[0,0])

for i in range(len(reaction_list)):
    print(f'{i}\t{np.sqrt(uncertainty.kinetic_covariance_matrix[i, i])}\t{reaction_list[i]}')

uncertainty.thermo_covariance_matrix[0,0] * 4184 * 4184

varH = np.float_power(28945.5, 2.0)

varH

type(uncertainty.reaction_list[10].kinetics)


