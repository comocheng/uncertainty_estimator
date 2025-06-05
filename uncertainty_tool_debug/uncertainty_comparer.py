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
# script to compare different methods of uncertainty estimation

# +
import os
import sys
import pickle
import copy
import numpy as np
import rmgpy.chemkin
import rmgpy
import rmgpy.tools.uncertainty
import rmgpy.kinetics.uncertainties

import random

import rmgpy.kinetics
import matplotlib.pyplot as plt
# %matplotlib inline
# -

def reactions_in_same_direction(reactionA, reactionB):
    reactantsA = [x.smiles for x in reactionA.reactants]
    reactantsB = [x.smiles for x in reactionB.reactants]
        
    return reactantsA[0] in reactantsB


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



def get_node_std(rxns, family):
    if len(rxns) == 1:
        print('NO SIDT EST')
        return np.nan
#         return 1.329
    recipe = database.kinetics.families[family].forward_recipe

    rxns = np.array(rxns)

    label = ''
    Tref = 1000.0
    data_mean = np.mean(np.log([r.kinetics.get_rate_coefficient(Tref) for r in rxns]))

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

    std_dev = kin_uncertainty.get_expected_log_uncertainty() / .398

    return std_dev

# ## Load database

# +
database = rmgpy.data.rmg.RMGDatabase()

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

database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = thermo_libraries,
    transport_libraries = [],
    reaction_libraries = [],
    seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
    kinetics_families = 'all',
    kinetics_depositories = ['training'],
    #frequenciesLibraries = self.statmechLibraries,
    depository = False, # Don't bother loading the depository information, as we don't use it
)

# +
# load aramco
aramco_chemkin_file = '/home/moon/autoscience/aramco/chem_annotated.inp'
aramco_dict_file = '/home/moon/autoscience/aramco/species_dictionary.txt'

species_listA, reaction_listA = rmgpy.chemkin.load_chemkin_file(aramco_chemkin_file, aramco_dict_file)

# -

def get_aramco_rxn_i(rxn):
    for i in range(len(reaction_listA)):
        if rxn.is_isomorphic(reaction_listA[i]):
            return i
    return -1


# +
# Load the model

# Must use annotated chemkin file
# chemkin_file = 'RMG-MAX1/chem_annotated.inp'
# dict_file = 'RMG-MAX1/species_dictionary.txt'

# chemkin_file = 'RMG-min-7/chem_annotated.inp'
# dict_file = 'RMG-min-7/species_dictionary.txt'

chemkin_file = 'ethane_original_db/chemkin/chem_annotated.inp'
dict_file = 'ethane_original_db/chemkin/species_dictionary.txt'

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file)


# Run Gao estimation of input parameters (takes a long time to load database)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations')
uncertainty.load_model(chemkin_file, dict_file)


thermo_libs = [
    'primaryThermoLibrary',
]

kinetic_libs = [
]

uncertainty.load_database(
    thermo_libraries=thermo_libs,
    kinetics_families='default',
    reaction_libraries=kinetic_libs,
)
uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()
# -

# Create a giant dictionary with all of the reaction family information in it
auto_gen_families = {}
for family_name in database.kinetics.families.keys():
    if family_name == 'Intra_R_Add_Endocyclic' or family_name == 'Intra_R_Add_Exocyclic':
        continue
    if database.kinetics.families[family_name].auto_generated and family_name not in auto_gen_families.keys():
        auto_gen_families[family_name] = database.kinetics.families[family_name].rules.get_entries()
        auto_gen_families[f'{family_name}_labels'] = [entry.label for entry in database.kinetics.families[family_name].rules.get_entries()]
        auto_gen_families[f'{family_name}_rxn_map'] = database.kinetics.families[family_name].get_reaction_matches(
            thermo_database=database.thermo,
            remove_degeneracy=True,
            get_reverse=True,
            exact_matches_only=False,
            fix_labels=True)


# # For each reaction, plot the kinetics and uncertainty

# +
# Plot BM node variance uncertainty, Gao uncertainty, 
# aramco equation?
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# for i in [0, 10, 20]:
for i in range(len(uncertainty.reaction_list)):
    # Plot kinetics
    Tref = 1000.0
    T = np.linspace(300, 3000, 1001)
    P = 101325
    k = np.zeros(len(T))
    plt.clf()
    for j in range(0, len(T)):
        k[j] = uncertainty.reaction_list[i].get_rate_coefficient(T[j])
    plt.plot(1000.0 / T, k, label=f'R{i}', color=colors[1], alpha=0.7)

    
    gao_sigma_lnk = uncertainty.kinetic_input_uncertainties[i]
    gao_sigma_k = np.exp(gao_sigma_lnk)


    # Plot node std dev
    plt.fill_between(1000.0 / T, k, k * gao_sigma_k, alpha=0.5, color=colors[1], edgecolor=None, label='Gao std dev')
    plt.fill_between(1000.0 / T, k / gao_sigma_k, k, alpha=0.5, color=colors[1], edgecolor=None)

    # plot Matt's uncertainty if it's a BM tree node...
    if 'Rate Rules' in uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]].keys() and \
            uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']:
        family = uncertainty.reaction_list[i].family
        if family != 'Intra_R_Add_Endocyclic' and family != 'Intra_R_Add_Exocyclic':
            node = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']
            rxns = auto_gen_families[f'{family}_rxn_map'][node]
#             sigma_lnk = get_node_std(rxns, family)
            if len(rxns) == 1:
                print('NO SIDT EST')
            sigma_lnk = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_std_dev']
            
            if np.isnan(sigma_lnk):
                continue
            sigma_k = np.exp(sigma_lnk)
            plt.fill_between(1000.0 / T, k, k * sigma_k, alpha=0.5, color=colors[2], edgecolor=None, label='Matt SIDT std dev')
            plt.fill_between(1000.0 / T, k / sigma_k, k, alpha=0.5, color=colors[2], edgecolor=None)
        
        
            # also plot training reactions
            for z in range(len(auto_gen_families[f'{family}_rxn_map'][node])):
                rxn = auto_gen_families[f'{family}_rxn_map'][node][z]
                k = np.zeros(len(T))
                for j in range(0, len(T)):
                    assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
                    k[j] = rxn.get_rate_coefficient(T[j])
                if z == 0:
                    plt.plot(1000.0 / T, k, label='Training Reactions', color=colors[0], alpha=0.5)
                else:
                    plt.plot(1000.0 / T, k, label='_nolegend_', color=colors[0], alpha=0.5)
                plt.yscale('log')
        
    print(i)
    display(uncertainty.reaction_list[i])
    print(uncertainty.kinetic_input_uncertainties[i])
    print(uncertainty.reaction_list[i].kinetics)
        
    # look for Aramco version
    a = get_aramco_rxn_i(uncertainty.reaction_list[i])
    P = 101325
    if a >= 0:
        # reverse the reaction if need be:
        rxnA = reaction_listA[a]
        if not reactions_in_same_direction(uncertainty.reaction_list[i], rxnA):
            print('REVERSING')
            rxnA = reaction_listA[a].generate_reverse_rate_coefficient()
                                           
        
        k = np.zeros(len(T))
        for j in range(0, len(T)):
            k[j] = rxnA.get_rate_coefficient(T[j], P)
        plt.plot(1000.0 / T, k, label=f'Aramco', color='black')
    
    
    
    plt.title(f'U log_10 k = {np.round(np.log10(gao_sigma_k), 3)}')
    plt.xlabel ('1000 K / T')
    plt.ylabel('k')
    plt.yscale('log')
    plt.legend()
    plt.show()
    print()
    print()
#     print(uncertainty.reaction_list[i])


# -

len(auto_gen_families[f'{uncertainty.reaction_list[21].family}_rxn_map'][node])



# # Plot reaction kinetics using all parent nodes

# +
# Plot BM node variance uncertainty, Gao uncertainty, 
# aramco equation?
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# for i in range(len(uncertainty.reaction_list)):
for i in [20]:

    # Plot kinetics
    Tref = 1000.0
    T = np.linspace(300, 3000, 1001)
    P = 101325
    k = np.zeros(len(T))
    plt.clf()
    for j in range(0, len(T)):
        k[j] = uncertainty.reaction_list[i].get_rate_coefficient(T[j])
    plt.plot(1000.0 / T, k, label=f'R{i}', color=colors[1], alpha=0.7)

    
    gao_sigma_lnk = uncertainty.kinetic_input_uncertainties[i]
    gao_sigma_k = np.exp(gao_sigma_lnk)


    # Plot node std dev
    plt.fill_between(1000.0 / T, k, k * gao_sigma_k, alpha=0.5, color=colors[1], edgecolor=None, label='Gao std dev')
    plt.fill_between(1000.0 / T, k / gao_sigma_k, k, alpha=0.5, color=colors[1], edgecolor=None)

    # plot Matt's uncertainty if it's a BM tree node...
    if 'Rate Rules' in uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]].keys() and \
            uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']:
        family = uncertainty.reaction_list[i].family
        if family != 'Intra_R_Add_Endocyclic' and family != 'Intra_R_Add_Exocyclic':
            node = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']
            rxns = auto_gen_families[f'{family}_rxn_map'][node]
#             sigma_lnk = get_node_std(rxns, family)
            if len(rxns) == 1:
                print('NO SIDT EST')
            sigma_lnk = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_std_dev']
            
            if np.isnan(sigma_lnk):
                continue
            sigma_k = np.exp(sigma_lnk)
            plt.fill_between(1000.0 / T, k, k * sigma_k, alpha=0.5, color=colors[2], edgecolor=None, label='Matt SIDT std dev')
            plt.fill_between(1000.0 / T, k / sigma_k, k, alpha=0.5, color=colors[2], edgecolor=None)
        
        
            # also plot training reactions
            for z in range(len(auto_gen_families[f'{family}_rxn_map'][node])):
                rxn = auto_gen_families[f'{family}_rxn_map'][node][z]
                k = np.zeros(len(T))
                for j in range(0, len(T)):
                    assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
                    k[j] = rxn.get_rate_coefficient(T[j])
                if z == 0:
                    plt.plot(1000.0 / T, k, label='Training Reactions', color=colors[0], alpha=0.5)
                else:
                    plt.plot(1000.0 / T, k, label='_nolegend_', color=colors[0], alpha=0.5)
                plt.yscale('log')
        
    print(i)
    display(uncertainty.reaction_list[i])
    print(uncertainty.kinetic_input_uncertainties[i])
    print(uncertainty.reaction_list[i].kinetics)
        
    # look for Aramco version
    a = get_aramco_rxn_i(uncertainty.reaction_list[i])
    P = 101325
    if a >= 0:
        # reverse the reaction if need be:
        rxnA = reaction_listA[a]
        if not reactions_in_same_direction(uncertainty.reaction_list[i], rxnA):
            print('REVERSING')
            rxnA = reaction_listA[a].generate_reverse_rate_coefficient()
                                           
        
        k = np.zeros(len(T))
        for j in range(0, len(T)):
            k[j] = rxnA.get_rate_coefficient(T[j], P)
        plt.plot(1000.0 / T, k, label=f'Aramco', color='black')
    
    
    
#     plt.title(f'U log_10 k = {np.round(np.log10(gao_sigma_k), 3)}')
    plt.title(node)
    plt.xlabel ('1000 K / T')
    plt.ylabel('k')
    plt.yscale('log')
    plt.legend()
    plt.show()
    print()
    print()
    
    
    # plot parent nodes
    if 'Rate Rules' in uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]].keys() and \
            uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node']:
        family = uncertainty.reaction_list[i].family
        if family != 'Intra_R_Add_Endocyclic' and family != 'Intra_R_Add_Exocyclic':
#             parent_node = uncertainty.database.kinetics.families[family].groups.entries[node].parent.label
            parent_node = node
            while parent_node is not None:
                Tref = 1000.0
                kineticsBM = uncertainty.database.kinetics.families[family].rules.entries[parent_node][0].data
                kinetics = kineticsBM.to_arrhenius(uncertainty.reaction_list[i].get_enthalpy_of_reaction(Tref))
                
                T = np.linspace(300, 3000, 1001)
                P = 101325
                k = np.zeros(len(T))
                official_k = np.zeros(len(T))
                plt.clf()
                
                for j in range(0, len(T)):
                    k[j] = kinetics.get_rate_coefficient(T[j]) * uncertainty.reaction_list[i].degeneracy
#                     k[j] = kineticsBM.to_arrhenius(uncertainty.reaction_list[i].get_enthalpy_of_reaction(Tref)).get_rate_coefficient(T[j])
                plt.plot(1000.0 / T, k, label=f'R{i}-node', color='red', linestyle='dashed', alpha=1.0, zorder=10)
                
                for j in range(0, len(T)):
#                     k[j] = uncertainty.database.kinetics.families[family].rules.entries[parent_node][0].data
                    k[j] = uncertainty.reaction_list[i].get_rate_coefficient(T[j])
                    official_k[j] = uncertainty.reaction_list[i].get_rate_coefficient(T[j])
                plt.plot(1000.0 / T, k, label=f'R{i}', color=colors[1], alpha=1.0, zorder=10)
        

                gao_sigma_lnk = uncertainty.kinetic_input_uncertainties[i]
                gao_sigma_k = np.exp(gao_sigma_lnk)


                # Plot node std dev
                plt.fill_between(1000.0 / T, k, k * gao_sigma_k, alpha=0.5, color=colors[1], edgecolor=None, label='Gao std dev')
                plt.fill_between(1000.0 / T, k / gao_sigma_k, k, alpha=0.5, color=colors[1], edgecolor=None)
                
                
                rxns = auto_gen_families[f'{family}_rxn_map'][parent_node]
                for z in range(len(rxns)):
                    rxn = auto_gen_families[f'{family}_rxn_map'][parent_node][z]
                    k = np.zeros(len(T))
                    for j in range(0, len(T)):
                        assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
                        k[j] = rxn.get_rate_coefficient(T[j])
                    if z == 0:
                        plt.plot(1000.0 / T, k, label='Training Reactions', color=colors[0], alpha=0.25)
                    else:
                        plt.plot(1000.0 / T, k, label='_nolegend_', color=colors[0], alpha=0.25)
                    plt.yscale('log')
                    
                
                if len(rxns) == 1:
                    print('NO SIDT EST')
                sigma_lnk = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]['Rate Rules'][1]['node_std_dev']

                if np.isnan(sigma_lnk):
                    continue
                sigma_k = np.exp(sigma_lnk)
                plt.fill_between(1000.0 / T, official_k, official_k * sigma_k, alpha=0.75, color=colors[2], edgecolor=None, label='Matt SIDT std dev')
                plt.fill_between(1000.0 / T, official_k / sigma_k, official_k, alpha=0.75, color=colors[2], edgecolor=None)
        
                
                    
                a = get_aramco_rxn_i(uncertainty.reaction_list[i])
                P = 101325
                if a >= 0:
                    # reverse the reaction if need be:
                    rxnA = reaction_listA[a]
                    if not reactions_in_same_direction(uncertainty.reaction_list[i], rxnA):
                        print('REVERSING')
                        rxnA = reaction_listA[a].generate_reverse_rate_coefficient()


                    k = np.zeros(len(T))
                    for j in range(0, len(T)):
                        k[j] = rxnA.get_rate_coefficient(T[j], P)
                    plt.plot(1000.0 / T, k, label=f'Aramco', color='black')
                    
                    
                plt.title(parent_node)
                plt.show()
                
                parent_entry = uncertainty.database.kinetics.families[family].groups.entries[parent_node].parent
                if parent_entry:
                    parent_node = parent_entry.label
                else:
                    parent_node = None


                    

# -

# # debug difference in kinetics from BM vs chemkin

# for i in range(len(uncertainty.reaction_list)):
for i in [20]:
    if 'Root' in uncertainty.reaction_list[i].kinetics.comment:
        family = uncertainty.reaction_list[i].family
        exact, source = uncertainty.database.kinetics.families[family].extract_source_from_comments(uncertainty.reaction_list[i])
        if not exact:
            node = source[1]['node']
            print(i, uncertainty.reaction_list[i].kinetics.comment, node)
            print()

get_aramco_rxn_i(uncertainty.reaction_list[20])

uncertainty.reaction_list[20].family





uncertainty.database.kinetics.families[family].groups.entries[node].parent.label

parent_node

uncertainty.database.kinetics.families[family].rules.entries[node][0].data.get_rate_coefficient(uncertainty.reaction_list[i].get_enthalpy_of_reaction(1000))

uncertainty.reaction_list[i].get_enthalpy_of_reaction(1000)



uncertainty.reaction_list[21]

uncertainty.kinetic_input_uncertainties[21]

uncertainty.reaction_list[21].kinetics

uncertainty.reaction_list[21].family

get_node_std(auto_gen_families['Singlet_Carbene_Intra_Disproportionation_rxn_map']['Root'], 'Singlet_Carbene_Intra_Disproportionation')


uncertainty.database.kinetics.families['Singlet_Carbene_Intra_Disproportionation'].rules.entries['Root'][0].data

len(auto_gen_families['Singlet_Carbene_Intra_Disproportionation_rxn_map']['Root'])

uncertainty.reaction_list[19].kinetics

uncertainty.kinetic_input_uncertainties[19]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[19]]







get_aramco_rxn_i(uncertainty.reaction_list[0])

reaction_listA[188].get_rate_coefficient(1000, 101325)

type(reaction_list[1].kinetics)

uncertainty.reaction_sources_dict[uncertainty.reaction_list[0]]['Rate Rules'][1]['node']





reaction_list[0].kinetics

# +
# For each BM node, show reaction with uncertainty and training reaction...

# then also do this for catalysis



# -





# +
# pick some example reactions to estimate

for i in range(len(reaction_list)):
    if not hasattr(reaction_list[i], 'family'):
        continue
        
        
    z = 80
print(reaction_list[z].family)
display(reaction_list[z])
print(reaction_list[z].kinetics)

my_rxn = copy.deepcopy(reaction_list[z])

print(my_rxn.reactants[0].molecule[0].get_all_labeled_atoms())
print(my_rxn.reactants[1].molecule[0].get_all_labeled_atoms())
print()

database.kinetics.families[reaction_list[z].family].add_atom_labels_for_reaction(my_rxn)

print(my_rxn.reactants[0].molecule[0].get_all_labeled_atoms())
print(my_rxn.reactants[1].molecule[0].get_all_labeled_atoms())
print()

template_labels = database.kinetics.families[my_rxn.family].get_reaction_template_labels(my_rxn)
print(template_labels)
template = database.kinetics.families[my_rxn.family].retrieve_template(template_labels)
print(template)
print()

kinetics = database.kinetics.families[my_rxn.family].get_kinetics_for_template(template, degeneracy=my_rxn.degeneracy)[0]
print(kinetics)

plot_kinetics([reaction_list[z], my_rxn], labels=['Original', 'refitted'])

# -










