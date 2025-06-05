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
import sys
import pickle
import copy
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty
import rmgpy.kinetics.uncertainties

import random

import rmgpy.kinetics
import matplotlib.pyplot as plt
# %matplotlib inline
# -

# ## Settings for RMG-MAX-1

# + active=""
# Global RMG Settings:
#    database.directory   = /home/sevy/rmg/RMG-database/input (Default, relative to RMG-Py source code)
#    test_data.directory  = /home/sevy/rmg/RMG-Py/test/rmgpy/test_data (Default, relative to RMG-Py source code)
# RMG execution initiated at Mon Apr  1 09:21:32 2024
#
# #########################################################
# # RMG-Py - Reaction Mechanism Generator in Python       #
# # Version: 3.2.0                                        #
# # Authors: RMG Developers (rmg_dev@mit.edu)             #
# # P.I.s:   William H. Green (whgreen@mit.edu)           #
# #          Richard H. West (r.west@neu.edu)             #
# # Website: http://reactionmechanismgenerator.github.io/ #
# #########################################################
#
# The current git HEAD for RMG-Py is:
#     b'23473858ac31733ae0703d6ea692cdb5a83d7172'
#     b'Mon Apr 1 08:41:05 2024 -0400'
#
# The current git HEAD for RMG-database is:
#     b'cfb4910cfc64274981216acbfb7756aba6be0112'
#     b'Tue Dec 12 11:55:13 2023 -0500'
#
# Reading input file "/home/sevy/rmg/autoscience/with_lib/butane_top10_20240401/input.py"...
# # Data sources
#
# thermo_libs = [
#     'BurkeH2O2',
#     'primaryThermoLibrary',
#     'FFCM1(-)',
#     'CurranPentane',
#     'Klippenstein_Glarborg2016',
#     'thermo_DFT_CCSDTF12_BAC',
#     'DFT_QCI_thermo',
#     'CBS_QB3_1dHR',
# ]
#
# kinetic_libs = [
#     'FFCM1(-)',
#     'CurranPentane',
#     'combustion_core/version5',
#     'Klippenstein_Glarborg2016',
#     'C2H4+O_Klipp2017',
#     'BurkeH2O2inArHe',
#     'BurkeH2O2inN2',
# ]
# -

def plot_kinetics(rxns, labels=None):
    """Function for plotting reaction kinetics
    Takes in a list of RMG reactions (rmgpy.reaction.Reaction) or a single reaction
    """
    plt.xlabel('1000 / T (K^-1)')
    plt.ylabel('log10(k)')

    if type(rxns) != list:
        rxns = [rxns]

    T = np.linspace(300, 3000, 1001)
    for rxn in rxns:
        k = np.zeros(len(T))
        for i in range(0, len(T)):
            k[i] = rxn.get_rate_coefficient(T[i], 101325)
        plt.plot(1000.0 / T, np.log10(k))

    if labels:
        plt.legend(labels)
    plt.show()



# ## Run Gao Uncertainty

# +
# Load the model

# Must use annotated chemkin file
chemkin_file = 'RMG-MAX1/chem_annotated.inp'
dict_file = 'RMG-MAX1/species_dictionary.txt'

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file)


# +
# Run Gao estimation of input parameters (takes a long time to load database)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations')
uncertainty.load_model(chemkin_file, dict_file)


thermo_libs = [
    'BurkeH2O2',
    'primaryThermoLibrary',
    'FFCM1(-)',
    'CurranPentane',
    'Klippenstein_Glarborg2016',
    'thermo_DFT_CCSDTF12_BAC',
    'DFT_QCI_thermo',
    'CBS_QB3_1dHR',
]

kinetic_libs = [
    'FFCM1(-)',
    'CurranPentane',
    'combustion_core/version5',
    'Klippenstein_Glarborg2016',
    'C2H4+O_Klipp2017',
    'BurkeH2O2inArHe',
    'BurkeH2O2inN2',
]

uncertainty.load_database(
    thermo_libraries=thermo_libs,
    kinetics_families='default',
    reaction_libraries=kinetic_libs,
)
uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()

# +
# BM Fitting notebook

# reload with BM fitting thermo libraries

thermo_libraries = [
    'Klippenstein_Glarborg2016',
    'BurkeH2O2',
    'thermo_DFT_CCSDTF12_BAC', 
    'DFT_QCI_thermo',
    'primaryThermoLibrary',
    'primaryNS',
    'NitrogenCurran',
    'NOx2018',
    'FFCM1(-)',
    'SulfurLibrary',
    'SulfurGlarborgH2S',
    'SABIC_aromatics',
]

uncertainty.load_database(
    thermo_libraries=thermo_libraries,
    kinetics_families='default',
    reaction_libraries=[],
)
# -





# +
# with open('RMG-MAX1/rxn_sources.pickle', 'wb') as f:
#     pickle.dump(uncertainty.reaction_sources_dict, f)
    
other_dict = {}
for i in range(len(uncertainty.reaction_list)):
    other_dict[i] = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
with open('RMG-MAX1/rxn_sources.pickle', 'wb') as f:
    pickle.dump(other_dict, f)
    
    
with open('RMG-MAX1/rxn_uncertainties.pickle', 'wb') as f:
    pickle.dump(uncertainty.kinetic_input_uncertainties, f)
# -



for i in range(0, 200):
    if 'Rate Rules' in uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]].keys():
        print(i, uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]])
        print()

# +
# with open('RMG-MAX1/uncertainty.pickle', 'wb') as f:  # can technically do this, 
#but it saves a full copy of the .pickle file so actually not really
#     pickle.dump(uncertainty, f)

# +
# Go through and summarize kinetic sources?

for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Library' in source:
        pass
    if 'PDep' in source:
        pass
#         print('pdep')
    if 'Training' in source:
        pass # this is usually (always?) an exact match
#         print(source)
    if 'Rate Rules' in source:
        if source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['rules']:
            pass
#             print('both')
        if source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['template']:
            pass
#             print('both of these')
#         print(source)

# +
# Count exact
n_library = 0
n_pdep = 0
n_training = 0
n_rr = 0
n_rr_training_and_template = 0
n_rr_training_and_no_template = 0
n_rr_no_training_and_yes_template = 0
n_sidt = 0


# get maximum non-SIDT uncertainty estimate
max_u = 0
i_max_u = -1

sidt_us = []
sidt_vars = []
sidt_Ns = []
empirical_us = []

for i in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    if 'Library' in source:
        n_library += 1
    if 'PDep' in source:
        n_pdep += 1
    if 'Training' in source:
        n_training += 1
    if 'Rate Rules' in source:
        assert source['Rate Rules'][1]['template'] or source['Rate Rules'][1]['node']
        
        n_rr += 1
        if source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['rules']:
            print('both')
        if source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['template']:
            n_rr_training_and_template += 1
            
            # N_template matches # top-level nodes on tree, which is often but not always the number of reactants
#             N_template = len(uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['template'])
#             if N_template != len(uncertainty.reaction_list[i].reactants):
#                 print(i)
        if not source['Rate Rules'][1]['training'] and source['Rate Rules'][1]['template']:
            if uncertainty.kinetic_input_uncertainties[i] != 1.5:
                print(i)
            
            n_rr_no_training_and_yes_template += 1
        if source['Rate Rules'][1]['node']:
#             # find an example with few training reactions in the Disproportionation family
#             if (int(source['Rate Rules'][1]['node_n_train']) >= 3 and int(source['Rate Rules'][1]['node_n_train']) <=6 ) and \
#                 uncertainty.reaction_list[i].kinetics.Ea.value_si > 0:
#                 print(i)
            n_sidt += 1
            sidt_us.append(uncertainty.kinetic_input_uncertainties[i])
            sidt_vars.append(source['Rate Rules'][1]['node_std_dev'])
            sidt_Ns.append(source['Rate Rules'][1]['node_n_train'])
        else:
            empirical_us.append(uncertainty.kinetic_input_uncertainties[i])
            if uncertainty.kinetic_input_uncertainties[i] > max_u:
                max_u = uncertainty.kinetic_input_uncertainties[i]
                i_max_u = i
            
            
print(f'Library:\t\t{n_library}')
print(f'PDEP:\t\t\t{n_pdep}')
print(f'Training:\t\t{n_training}')
print(f'Rate Rules:\t\t{n_rr}')
print(f'n_rr_training_and_template:\t{n_rr_training_and_template}')
print(f'n_rr_training_and_no_template:\t{n_rr_training_and_no_template}')
print(f'n_rr_no_training_and_yes_template:\t{n_rr_no_training_and_yes_template}')
print(f'SIDT:\t{n_sidt}')
print(f'Manual RR Count:\t{n_sidt+n_rr_no_training_and_yes_template+n_rr_training_and_no_template+n_rr_training_and_template}')


# -

np.any(np.isnan(sidt_vars))

plt.scatter(sidt_Ns, sidt_vars, alpha=0.1)

len(sidt_us)

len(empirical_us)

# +
plt.hist(sidt_vars, 64, label='SIDT')
plt.hist(empirical_us, alpha=0.75, label='Empirical')
plt.legend()

print(np.median(sidt_vars))
print(np.median(empirical_us))
# -







max_u

i_max_u

# ### Library Example

i = 0
display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')
print(f'Library: {uncertainty.reaction_list[i].family}')


# ### PDEP Example

i = 138
display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')


# ### Training Example

# +
i = 83
display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')
print(f'Family: {uncertainty.reaction_list[i].family}')



# -

# ### Rate Rules 1: Templates and Training

# +
i = 135
i = 141
display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')
print(f'Family: {uncertainty.reaction_list[i].family}')
print()


M_train = len(uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['training'])
lnk = 0.5 + 1.0 + np.log10(M_train + 1) * 3.5
print(f'lnk: {lnk}')
# -

# ### Rate Rules 2: Templates and Rules

# +
i = 136
i = 240
i = 242
i = 187
display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')
print(f'Family: {uncertainty.reaction_list[i].family}')
print()


M_train = len(uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['training'])
lnk = 0.5 + 1.0 + np.log10(M_train + 1) * 3.5
print(f'lnk: {lnk}')
# -

# ### Rate Rules 3: SIDT Trees

# +
i = 80
display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')
print(f'Family: {uncertainty.reaction_list[i].family}')

N = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_n_train']
s = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_std_dev']
lnk = s + 1.0 + np.log10(N + 1) * 3.5
print(f'lnk: {lnk}')
# -



# +
# Need to get training reactions used in each node

# +
# for entry in uncertainty.database.kinetics.families[uncertainty.reaction_list[80].family].rules.entries:
#     if 'Root_Ext-1R!H-R_' in entry:
#         print(entry)
# -











uncertainty.kinetic_input_uncertainties[83]

uncertainty.reaction_list[138].kinetics

uncertainty.reaction_list[0].kinetics

uncertainty.reaction_sources_dict[uncertainty.reaction_list[0]]

uncertainty.reaction_list[0]

source['Library']

source['Rate Rules']



uncertainty.reaction_list[-2].kinetics

uncertainty.reaction_sources_dict[uncertainty.reaction_list[-2]]

uncertainty.kinetic_input_uncertainties[-2]

uncertainty.reaction_list[-2]



i = 95
node_name = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']
uncertainty.database.kinetics.families[uncertainty.reaction_list[i].family].rules.entries[node_name][0].long_desc
#                     std_dev_matches = re.search(r'Standard Deviation in ln\(k\): ([0-9]*.[0-9]*)', long_desc)
#                     std_dev = -1.0
#                     if std_dev_matches is not None:
#                         std_dev = float(std_dev_matches[1])

# # Rebuild tree to get training reactions in each node

# Create a giant dictionary with all of the reaction family information in it
auto_gen_families = {}
for family_name in uncertainty.database.kinetics.families.keys():
    if family_name == 'Intra_R_Add_Endocyclic' or family_name == 'Intra_R_Add_Exocyclic':
        continue
    if uncertainty.database.kinetics.families[family_name].auto_generated and family_name not in auto_gen_families.keys():
        auto_gen_families[family_name] = uncertainty.database.kinetics.families[family_name].rules.get_entries()
        auto_gen_families[f'{family_name}_labels'] = [entry.label for entry in uncertainty.database.kinetics.families[family_name].rules.get_entries()]
        auto_gen_families[f'{family_name}_rxn_map'] = uncertainty.database.kinetics.families[family_name].get_reaction_matches(
            thermo_database=uncertainty.database.thermo,
            remove_degeneracy=True,
            get_reverse=True,
            exact_matches_only=False,
            fix_labels=True)


# +
# i = 1011
i = 1549
# i = 1554
# i = 1559
# i = 1561
# 1563
# 1571
# 1575
# 1576
# 1581
# 1584
# 1591
# 1593
# i = 1599
# i = 1734
# i = 1935
# 2163
# 2301
# i = 2413

display(uncertainty.reaction_list[i])
print(f'Source Entry: {uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]}')
print(f'Uncertainty: {uncertainty.kinetic_input_uncertainties[i]}')
print(f'Family: {uncertainty.reaction_list[i].family}')

node = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']
N = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_n_train']
s = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_std_dev']
lnk = s + 1.0 + np.log10(N + 1) * 3.5
print(f'lnk: {lnk}')

print()
print(uncertainty.reaction_list[i].kinetics)
# -



# +
rxn_mod_arr = uncertainty.reaction_list[i]
rxn_bm = copy.deepcopy(rxn_mod_arr)
rxn_bm.kinetics = uncertainty.database.kinetics.families[uncertainty.reaction_list[i].family].rules.entries[node][0].data


plt.xlabel('1000 / T (K^-1)')
plt.ylabel('k')

labels=['Mod Arrhenius', 'Blowers-Masel']
T = np.linspace(300, 3000, 1001)
for rxn in [rxn_mod_arr, rxn_bm]:
    k = np.zeros(len(T))
    for j in range(0, len(T)):
        if type(rxn.kinetics) == rmgpy.kinetics.arrhenius.ArrheniusBM:
#             k[j] = rxn.get_rate_coefficient(T[j], rxn.get_enthalpy_of_reaction(T[j]))
            k[j] = rxn.get_rate_coefficient(T[j], rxn.get_enthalpy_of_reaction(1000))
        else:
            k[j] = rxn.get_rate_coefficient(T[j])
    plt.plot(1000.0 / T, k)
plt.yscale('log')
if labels:
    plt.legend(labels)
plt.show()




# plot_kinetics([rxn_mod_arr, rxn_bm], labels=['Mod Arrhenius', 'Blowers-Masel'])
# -

plot_kinetics(auto_gen_families['Disproportionation_rxn_map'][node], labels=[1, 2, 3, 4, 5, 6])

auto_gen_families['Disproportionation_rxn_map'][node][0].kinetics

list(auto_gen_families['Disproportionation_rxn_map'].keys())

# +
# Given a node, I want to plot all the kinetics, the average kinetics, and the std_dev
family = 'Disproportionation'

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
node = list(auto_gen_families[f'{family}_rxn_map'].keys())[44]

# # Find nodes that matter
# for a in range(len(list(auto_gen_families[f'{family}_rxn_map'].keys()))):
#     node = list(auto_gen_families[f'{family}_rxn_map'].keys())[a]
#     reaction_indices = []
#     for z in range(len(uncertainty.reaction_list)):
#         source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[z]]
#         if 'Rate Rules' in source:
#             if uncertainty.reaction_list[z].family == family:
#                 if source['Rate Rules'][1]['node'] and source['Rate Rules'][1]['node'] == node:
#                     # plot this reaction
#                     reaction_indices.append(z)
#     if len(reaction_indices) > 0:
#         print(a, len(reaction_indices))



print(f'Node: {node}')
plt.xlabel('1000 / T (K^-1)')
plt.ylabel('k')


T = np.linspace(300, 3000, 1001)
for z in range(len(auto_gen_families[f'{family}_rxn_map'][node])):
    rxn = auto_gen_families[f'{family}_rxn_map'][node][z]
    k = np.zeros(len(T))
    for j in range(0, len(T)):
        assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
        k[j] = rxn.get_rate_coefficient(T[j])
    plt.plot(1000.0 / T, k, label=z, color=colors[0], alpha=0.5)
    print(f'{rxn.kinetics.Ea}\t\t{rxn.kinetics.n}')
    
    
# Also plot the average of the kinetics? Except that isn't even relevant to the BM tree?
# It's an approximation. Good enough for just looking at the trees

# Maybe I plot the reactions that are estimated using that node...
reaction_indices = []
for z in range(len(uncertainty.reaction_list)):
    source = uncertainty.reaction_sources_dict[uncertainty.reaction_list[z]]
    if 'Rate Rules' in source:
        if uncertainty.reaction_list[z].family == family:
            if source['Rate Rules'][1]['node'] and source['Rate Rules'][1]['node'] == node:
                # plot this reaction
                reaction_indices.append(z)
                
                
                rxn = uncertainty.reaction_list[z]
                k = np.zeros(len(T))
                for j in range(0, len(T)):
                    assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
                    k[j] = rxn.get_rate_coefficient(T[j])
                plt.plot(1000.0 / T, k, label=z, color='black', alpha=0.1)
                
#                 break
#                 print(z)
print(f'Reactions Estimated with node: {len(reaction_indices)}')            
    
plt.yscale('log')

# plt.legend()
plt.show()




# -

node

uncertainty.database.kinetics.families[uncertainty.reaction_list[i].family].rules.entries[node][0].data

rxns = auto_gen_families[f'{family}_rxn_map'][node]
recipe = uncertainty.database.kinetics.families[uncertainty.reaction_list[i].family].forward_recipe


random.seed(400)


random.shuffle(rxns)
display(rxns[0])

kin = rmgpy.kinetics.arrhenius.ArrheniusBM().fit_to_reactions(rxns, recipe=recipe.actions)

kin





kin





# how many nodes?
len(auto_gen_families['Disproportionation_rxn_map'])

# how many training reactions?
len(auto_gen_families['Disproportionation_rxn_map']['Root'])

# +
# Plot a particular SIDT estimate's reactions and uncertainty



# +
# Want to see if I can reproduce a node's rule using a recently regenerated Disproportionation Tree

nodes = list(uncertainty.database.kinetics.families['Disproportionation'].rules.entries)

node = nodes[0]
node = 'Root_Ext-1R!H-R_4R->O'

rxns = auto_gen_families['Disproportionation_rxn_map'][node]


recipe = uncertainty.database.kinetics.families['Disproportionation'].forward_recipe
kin = rmgpy.kinetics.arrhenius.ArrheniusBM().fit_to_reactions(rxns, recipe=recipe.actions)
# -

kin

uncertainty.database.kinetics.families['Disproportionation'].rules.entries[node][0].data

len(rxns)

(6.666 - 5.8119)/ 6.666



# +
# Try to run the leave-one-out code
recipe = uncertainty.database.kinetics.families['Disproportionation'].forward_recipe
rxns = auto_gen_families['Disproportionation_rxn_map'][node]

dlnks = np.array([
    np.log(
        rmgpy.kinetics.arrhenius.ArrheniusBM().fit_to_reactions(rxns[list(set(range(len(rxns))) - {i})], recipe=recipe)
        .to_arrhenius(rxn.get_enthalpy_of_reaction(Tref))
        .get_rate_coefficient(T=Tref) / rxn.get_rate_coefficient(T=Tref)
    ) for i, rxn in enumerate(rxns)
])

# -



# +
recipe = uncertainty.database.kinetics.families['Disproportionation'].forward_recipe
rxns = auto_gen_families['Disproportionation_rxn_map'][node]
rxns = np.array(rxns)

label = node

data_mean = np.mean(np.log([r.kinetics.get_rate_coefficient(Tref) for r in rxns]))
Tref = 1000.0
n = len(rxns)

dlnks = np.array([
    np.log(
        rmgpy.kinetics.arrhenius.ArrheniusBM().fit_to_reactions(rxns[list(set(range(len(rxns))) - {i})], recipe=recipe.actions)
        .to_arrhenius(rxn.get_enthalpy_of_reaction(Tref))
        .get_rate_coefficient(T=Tref) / rxn.get_rate_coefficient(T=Tref)
    ) for i, rxn in enumerate(rxns)
])


varis = (np.array([rmgpy.kinetics.uncertainties.rank_accuracy_map[rxn.rank].value_si for rxn in rxns]) / (2.0 * 8.314 * Tref)) ** 2
# weighted average calculations
ws = 1.0 / varis
V1 = ws.sum()
V2 = (ws ** 2).sum()
mu = np.dot(ws, dlnks) / V1
s = np.sqrt(np.dot(ws, (dlnks - mu) ** 2) / (V1 - V2 / V1))

kin_uncertainty = rmgpy.kinetics.uncertainties.RateUncertainty(mu=mu, var=s ** 2, N=n, Tref=Tref, data_mean=data_mean, correlation=label)


print(kin_uncertainty)
# -

uncertainty.database.kinetics.families['Disproportionation'].rules.entries[node][0].data

2.6655742044913753 / .398

kin_uncertainty.get_expected_log_uncertainty() / .398

np.sqrt(2.6655742044913753) 

# +
rxns = np.array(rxns)

[rxns[list(set(range(len(rxns))) - {i})] for i, rxn in enumerate(rxns)]
# -




