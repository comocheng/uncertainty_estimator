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
# Script to see how well your uncertainty scheme works for a given mechanism

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
import rmgpy.exceptions

import random

import rmgpy.kinetics
import matplotlib.pyplot as plt
# %matplotlib inline

sys.path.append('/home/moon/autoscience/reaction_calculator/database')
import database_fun


# -

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


def reactions_in_same_direction(reactionA, reactionB):
    reactantsA = [x.smiles for x in reactionA.reactants]
    reactantsB = [x.smiles for x in reactionB.reactants]
        
    return reactantsA[0] in reactantsB


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


# +
# load aramco
aramco_chemkin_file = '/home/moon/autoscience/aramco/chem_annotated.inp'
aramco_dict_file = '/home/moon/autoscience/aramco/species_dictionary.txt'

species_listA, reaction_listA = rmgpy.chemkin.load_chemkin_file(aramco_chemkin_file, aramco_dict_file)


# +
# load RMG-min-7
rmg_min_7_chemkin_file = '/home/moon/uncertainty_estimator/RMG-min-7/chem_annotated.inp'
rmg_min_7_dict_file = '/home/moon/uncertainty_estimator/RMG-min-7/species_dictionary.txt'

species_list7, reaction_list7 = rmgpy.chemkin.load_chemkin_file(rmg_min_7_chemkin_file, rmg_min_7_dict_file)

# -





# # Plot the reactions in the mechanism

# +
i = 280
# i = 124
# for i in disps[0:15]:

data_entry = {}

# collect data here and plot on next plot
SIDT_compare_idx = []
# RMG

# for i in range(500):
for i in range(50):
    my_rxn = copy.deepcopy(reaction_listA[i])
    plt.clf()
#     print(i)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    Tref = 1000.0
    T = np.linspace(300, 3000, 1001)
    P = 101325
    k = np.zeros(len(T))
    for j in range(0, len(T)):
        k[j] = my_rxn.get_rate_coefficient(T[j], P)
    plt.plot(1000.0 / T, k, label=f'Aramco', color='black')

    fam_rxn_list = database.kinetics.generate_reactions_from_families(
        reactants=my_rxn.reactants,
        products=my_rxn.products,
        only_families=None,
        resonance=True,
    )
    if not fam_rxn_list:
#         print('No matches')
        continue

    family = fam_rxn_list[0].family
    if family not in auto_gen_families.keys():
#         print(f'{family} is not autogenerated')
        continue
    
    
    SIDT_compare_idx.append(i)
    print(i)
    if not fam_rxn_list[0].is_forward:
        # skipping reverse for now?
        my_rxn.reactants = reaction_listA[i].products
        my_rxn.products = reaction_listA[i].reactants
        print('Reversing reaction')
        
        # plot aramco in reverse
        rev_kinA = reaction_listA[i].generate_reverse_rate_coefficient()
        for j in range(0, len(T)):
            k[j] = rev_kinA.get_rate_coefficient(T[j], P)
        plt.plot(1000.0 / T, k, label=f'Aramco-reverse', color='grey')
        
    my_rxn.degeneracy = fam_rxn_list[0].degeneracy
        
    for family in [r.family for r in fam_rxn_list]:
        database.kinetics.families[family].add_atom_labels_for_reaction(my_rxn)

        template_labels = database.kinetics.families[family].get_reaction_template_labels(my_rxn)
        template = database.kinetics.families[family].retrieve_template(template_labels)

        kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=my_rxn.degeneracy)[0]
    
    

    # also plot Matt's uncertainty, which means we need the families/nodes
    node = template[0].label
    rxns = auto_gen_families[f'{family}_rxn_map'][node]
    for z in range(len(rxns)):
        rxn = rxns[z]
        k = np.zeros(len(T))
        for j in range(0, len(T)):
            assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
            k[j] = rxn.get_rate_coefficient(T[j])
        if z == 0:
            plt.plot(1000.0 / T, k, label='Training Reactions', color=colors[0], alpha=0.5)
        else:
            plt.plot(1000.0 / T, k, label='_nolegend_', color=colors[0], alpha=0.5)
        plt.yscale('log')

    k = np.zeros(len(T))
    arrh_kinetics = kinetics.to_arrhenius(my_rxn.get_enthalpy_of_reaction(1000))
    for j in range(0, len(T)):
#         k[j] = kinetics.get_rate_coefficient(T[j], my_rxn.get_enthalpy_of_reaction(T[j]))
        k[j] = arrh_kinetics.get_rate_coefficient(T[j])
    plt.plot(1000.0 / T, k, label=f'RMG Estimate', color=colors[1], alpha=0.7)


    if len(rxns) == 1:
#         print('NO SIDT EST')
        pass
    else:
        r = rmgpy.reaction.Reaction()
        r.kinetics = kinetics

        sigma_lnk = get_node_std(rxns, family)

    #     sigma_lnk = database.kinetics.families[family].extract_source_from_comments(r)[1][1]['node_std_dev']
        sigma_k = np.exp(sigma_lnk)
        plt.fill_between(1000.0 / T, k, k * sigma_k, alpha=0.5, color=colors[2], edgecolor=None, label='Matt SIDT std dev')
        plt.fill_between(1000.0 / T, k / sigma_k, k, alpha=0.5, color=colors[2], edgecolor=None)    

    plt.title(node)
    plt.xlabel ('1000 K / T')
    plt.ylabel('k')
    plt.yscale('log')
    plt.legend()
    plt.show()
    # gao_sigma_lnk = uncertainty.kinetic_input_uncertainties[i]
    # gao_sigma_k = np.exp(gao_sigma_lnk)


    # # Plot node std dev
    # plt.fill_between(1000.0 / T, k, k * gao_sigma_k, alpha=0.5, color=colors[1], edgecolor=None, label='Gao std dev')
    # plt.fill_between(1000.0 / T, k / gao_sigma_k, k, alpha=0.5, color=colors[1], edgecolor=None)




# -







# # Debug Individual Node Estimate

# +
i = 0
print(i)
my_rxn = copy.deepcopy(reaction_listA[i])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
Tref = 1000.0
T = np.linspace(300, 3000, 1001)
P = 101325
k = np.zeros(len(T))
for j in range(0, len(T)):
    k[j] = my_rxn.get_rate_coefficient(T[j], P)
plt.plot(1000.0 / T, k, label=f'Aramco', color='black')

fam_rxn_list = database.kinetics.generate_reactions_from_families(
    reactants=my_rxn.reactants,
    products=my_rxn.products,
    only_families=None,
    resonance=True,
)
if not fam_rxn_list:
    print('No matches')

family = fam_rxn_list[0].family
if family not in auto_gen_families.keys():
    print(f'{family} is not autogenerated')


if not fam_rxn_list[0].is_forward:
    my_rxn.reactants = reaction_listA[i].products
    my_rxn.products = reaction_listA[i].reactants
    print('Reversing reaction')

my_rxn.degeneracy = fam_rxn_list[0].degeneracy
    
display(my_rxn)
for family in [r.family for r in fam_rxn_list]:
    database.kinetics.families[family].add_atom_labels_for_reaction(my_rxn)

    template_labels = database.kinetics.families[family].get_reaction_template_labels(my_rxn)
    template = database.kinetics.families[family].retrieve_template(template_labels)

    kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=my_rxn.degeneracy)[0]


# also plot Matt's uncertainty, which means we need the families/nodes
node = template[0].label
rxns = auto_gen_families[f'{family}_rxn_map'][node]
for z in range(len(rxns)):
    rxn = rxns[z]
    k = np.zeros(len(T))
    for j in range(0, len(T)):
        assert type(rxn.kinetics) == rmgpy.kinetics.arrhenius.Arrhenius
        k[j] = rxn.get_rate_coefficient(T[j])
    if z == 0:
        plt.plot(1000.0 / T, k, label='Training Reactions', color=colors[0], alpha=0.5)
    else:
        plt.plot(1000.0 / T, k, label='_nolegend_', color=colors[0], alpha=0.5)
    plt.yscale('log')
    
#     print('n_t=', rxn.kinetics.n)
#     print('A_t=', rxn.kinetics.A.value_si)
    print('Ea_t=', rxn.kinetics.Ea.value_si)
#     print('H_t=', rxn.get_enthalpy_of_reaction(1000))

k = np.zeros(len(T))
for j in range(0, len(T)):
    k[j] = kinetics.get_rate_coefficient(T[j], my_rxn.get_enthalpy_of_reaction(T[j]))
plt.plot(1000.0 / T, k, label=f'RMG Estimate', color=colors[1], alpha=0.7)


if len(rxns) == 1:
    print('NO SIDT EST')

else:
    r = rmgpy.reaction.Reaction()
    r.kinetics = kinetics

    sigma_lnk = get_node_std(rxns, family)
    print('my computed sigma is ', sigma_lnk)
    s_matt = r.kinetics.uncertainty.get_expected_log_uncertainty() / .398
    print('Matts sigma is ', s_matt)
    
    
#     print('n_est=', kinetics.n)
#     print('A_est=', kinetics.A.value_si)
    print('Ea_est=', kinetics.E0.value_si)
#     print('H_est=', my_rxn.get_enthalpy_of_reaction(1000))
    
    sigma_k = np.exp(sigma_lnk)
    plt.fill_between(1000.0 / T, k, k * sigma_k, alpha=0.5, color=colors[2], edgecolor=None, label='Matt SIDT std dev')
    plt.fill_between(1000.0 / T, k / sigma_k, k, alpha=0.5, color=colors[2], edgecolor=None)    

plt.title(node)
plt.xlabel ('1000 K / T')
plt.ylabel('k')
plt.yscale('log')
plt.legend()
plt.show()
# gao_sigma_lnk = uncertainty.kinetic_input_uncertainties[i]
# gao_sigma_k = np.exp(gao_sigma_lnk)


# # Plot node std dev
# plt.fill_between(1000.0 / T, k, k * gao_sigma_k, alpha=0.5, color=colors[1], edgecolor=None, label='Gao std dev')
# plt.fill_between(1000.0 / T, k / gao_sigma_k, k, alpha=0.5, color=colors[1], edgecolor=None)


# also plot the predicted kinetics for one of the training reactions
train0 = copy.deepcopy(rxns[0])
kin0 = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=rxns[0].degeneracy)[0]
k = np.zeros(len(T))
for j in range(0, len(T)):
    k[j] = kin0.get_rate_coefficient(T[j], train0.get_enthalpy_of_reaction(T[j]))
plt.plot(1000.0 / T, k, label=f'RMG Est Train 0', color=colors[5], alpha=0.7)


# -



# # Record the data 

# +
data_entries = []
data_entry = {}


# Tdatas = [300, 1000, 1500]
Tdatas = [1000]

for i in range(len(reaction_listA)):
    my_rxn = copy.deepcopy(reaction_listA[i])
    
    # get aramco or reverse aramco if applicable
    Tref = 1000.0
    P = 101325

    # get node estimate
    fam_rxn_list = database.kinetics.generate_reactions_from_families(
        reactants=my_rxn.reactants,
        products=my_rxn.products,
        only_families=None,
        resonance=True,
    )
    if not fam_rxn_list:
        continue

    family = fam_rxn_list[0].family
    if family not in auto_gen_families.keys():
        continue
    my_rxn.degeneracy = fam_rxn_list[0].degeneracy
    
    
    print(i)
    # get aramco kinetics in reverse if needed
    if not fam_rxn_list[0].is_forward:
        my_rxn.reactants = reaction_listA[i].products
        my_rxn.products = reaction_listA[i].reactants
        my_rxn.kinetics = reaction_listA[i].generate_reverse_rate_coefficient()
        
    # get node estimate
    for family in [r.family for r in fam_rxn_list]:
        database.kinetics.families[family].add_atom_labels_for_reaction(my_rxn)
        template_labels = database.kinetics.families[family].get_reaction_template_labels(my_rxn)
        template = database.kinetics.families[family].retrieve_template(template_labels)
        node = template[0].label
        kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=my_rxn.degeneracy)[0]
    arrh_kinetics = kinetics.to_arrhenius(my_rxn.get_enthalpy_of_reaction(Tref))

    # get the sigma from Matt's estimate
    rxns = auto_gen_families[f'{family}_rxn_map'][node]
    sigma_lnk = get_node_std(rxns, family)
    
    for T in Tdatas:
        # this will be a data entry
        data_entry = {
            'mech_idx': i,
            'mech_k': my_rxn.kinetics.get_rate_coefficient(T, P),
            'RMG_k': arrh_kinetics.get_rate_coefficient(T, P),
            'family': family,
            'node': node,
            'T': T,
            'reverse': not fam_rxn_list[0].is_forward, 
            'sigma_lnk': sigma_lnk
        }
        data_entries.append(data_entry)
    
    


# -

with open('my_data_entries', 'wb') as f:
    pickle.dump(data_entries, f)

# +
# Collect the sigma_lnks, RMG_k, mech_ks in to an array
RMG_ks = [x['RMG_k'] for x in data_entries]
mech_ks = [x['mech_k'] for x in data_entries]

plt.plot(np.log10(mech_ks), np.log10(mech_ks))
plt.scatter(np.log10(mech_ks), np.log10(RMG_ks), alpha=0.1)

plt.xlim([-10, 20])
plt.ylim([-10, 20])

ax = plt.gca()
ax.set_aspect('equal', 'box')

# +
# plot all the low sigma_k examples
RMG_ks = [x['RMG_k'] for x in data_entries if x['sigma_lnk'] < 1.329]
mech_ks = [x['mech_k'] for x in data_entries if x['sigma_lnk'] < 1.329]
sigmas = [np.log10(np.exp(x['sigma_lnk'])) for x in data_entries if x['sigma_lnk'] < 1.329]


plt.plot(np.log10(mech_ks), np.log10(mech_ks), color='black', linewidth=0.5, linestyle='dashed')
# plt.scatter(np.log10(mech_ks), np.log10(RMG_ks), alpha=0.1)


plt.errorbar(np.log10(mech_ks), np.log10(RMG_ks), xerr=sigmas, fmt='none', alpha=0.5)

plt.xlim([-10, 20])
plt.ylim([-10, 20])

ax = plt.gca()
ax.set_aspect('equal', 'box')

# +
# plot all the medium sigma_k examples
RMG_ks = [x['RMG_k'] for x in data_entries if (x['sigma_lnk'] >= 1.329 and x['sigma_lnk'] <= 10)]
mech_ks = [x['mech_k'] for x in data_entries if (x['sigma_lnk'] >= 1.329 and x['sigma_lnk'] <= 10)]
sigmas = [np.log10(np.exp(x['sigma_lnk'])) for x in data_entries if (x['sigma_lnk'] >= 1.329 and x['sigma_lnk'] <= 10)]


plt.plot(np.log10(mech_ks), np.log10(mech_ks), color='black', linewidth=0.5, linestyle='dashed')
# plt.scatter(np.log10(mech_ks), np.log10(RMG_ks), alpha=0.1)


plt.errorbar(np.log10(mech_ks), np.log10(RMG_ks), xerr=sigmas, fmt='none', alpha=0.5)

plt.xlim([-10, 20])
plt.ylim([-10, 20])

ax = plt.gca()
ax.set_aspect('equal', 'box')

# +
# plot all the high sigma_k examples
RMG_ks = [x['RMG_k'] for x in data_entries if x['sigma_lnk'] > 10]
mech_ks = [x['mech_k'] for x in data_entries if x['sigma_lnk'] > 10]
sigmas = [np.log10(np.exp(x['sigma_lnk'])) for x in data_entries if x['sigma_lnk'] > 10]


plt.plot(np.log10(mech_ks), np.log10(mech_ks), color='black', linewidth=0.5, linestyle='dashed')
# plt.scatter(np.log10(mech_ks), np.log10(RMG_ks), alpha=0.1)


plt.errorbar(np.log10(mech_ks), np.log10(RMG_ks), xerr=sigmas, fmt='none', alpha=0.5)

plt.xlim([-10, 20])
plt.ylim([-10, 20])

ax = plt.gca()
ax.set_aspect('equal', 'box')

# +
# Count up the times it's more than 95% out there
# -1.96 to 1.96 to

# Total

RMG_ks = [x['RMG_k'] for x in data_entries]
mech_ks = [x['mech_k'] for x in data_entries]
sigmas = [np.log10(np.exp(x['sigma_lnk'])) for x in data_entries]


# sort by the size of the sigma?
indices = np.arange(len(sigmas))





plt.plot(np.log10(mech_ks), np.log10(mech_ks), color='black', linewidth=0.5, linestyle='dashed')
# plt.scatter(np.log10(mech_ks), np.log10(RMG_ks), alpha=0.1)

outlier_count = 0
for i in range(len(RMG_ks)):
    if np.abs(np.log10(mech_ks[i]) - np.log10(RMG_ks[i])) > sigmas[i]:
        outlier_count += 1



plt.errorbar(np.log10(mech_ks), np.log10(RMG_ks), xerr=sigmas, fmt='none', alpha=0.5)

plt.xlim([-10, 20])
plt.ylim([-10, 20])

ax = plt.gca()
ax.set_aspect('equal', 'box')

print(f'{outlier_count}/{len(RMG_ks)}')
print(f'{np.round(outlier_count / len(RMG_ks), 5) * 100.0}%')


# -



# +
np.random.seed(400)
samples = np.random.normal(size=100000)

print(np.percentile(samples, 50))


# -

def get_confidence(s):
    # s is how many standard deviations
    return np.sum(np.abs(samples) < s) / len(samples)


# +
# Scale all results to a normalized gaussian curve

# Total
RMG_ks = np.array([x['RMG_k'] for x in data_entries])
mech_ks = np.array([x['mech_k'] for x in data_entries])
sigmas = np.array([np.log10(np.exp(x['sigma_lnk'])) for x in data_entries])

sigmas[sigmas == 0] = np.nan

# sigmas[np.isnan(sigmas)] = np.nanmax(sigmas)

err = [np.log10(RMG_ks[i]) - np.log10(mech_ks[i]) for i in range(len(RMG_ks))]
ax = plt.hist(samples[:len(err) - np.sum(np.isnan(err))], 100, range=[-10, 10], alpha=0.5)
ax = plt.hist(err, 100, range=[-10, 10], alpha=0.5)



# +
# Scale all results to a normalized gaussian curve
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# Total
RMG_ks = np.array([x['RMG_k'] for x in data_entries])
mech_ks = np.array([x['mech_k'] for x in data_entries])
sigmas = np.array([np.log10(np.exp(x['sigma_lnk'])) for x in data_entries])

sigmas[np.isnan(sigmas)] = np.inf
# sigmas[sigmas == 0] = np.nan

err = np.log10(mech_ks) - np.log10(RMG_ks)

indices = np.arange(len(sigmas))
sorted_order = [x for _, x in sorted(zip(sigmas, indices))][::-1]

plt.plot(sigmas[sorted_order], color=colors[0], label=r'Estimated 1$\sigma$')
plt.plot(-sigmas[sorted_order], color=colors[0])

plt.scatter(indices, err[sorted_order], marker='.', alpha=0.15, color='black', label='Reaction Error')

first_non_nan = np.argmax(np.invert(np.isinf(sigmas[sorted_order])))
plt.xlim([first_non_nan, len(indices)])

plt.xlabel('Data Index')
plt.ylabel('Error log10(k)')
plt.legend()
# -








