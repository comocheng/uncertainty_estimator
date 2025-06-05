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
# A notebook to diagnose problems in catalysis mechanisms

# +
import cantera as ct
import numpy as np
import os
import copy
import cycler

import rmgpy.chemkin

import matplotlib.pyplot as plt
# %matplotlib inline
# -

def get_i_thing(my_thing, thing_list):
    for j in range(len(thing_list)):
        if thing_list[j].is_isomorphic(my_thing):
            return j
    return -1


def reactions_in_same_direction(reactionA, reactionB):
    reactantsA = [x.smiles for x in reactionA.reactants]
    reactantsB = [x.smiles for x in reactionB.reactants]
        
    return reactantsA[0] in reactantsB


# +
# Plot the kinetics here

# match each reaction?
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



sevy_mech = './cpox_pt_20241101/pt55_og_lib_limited_fams'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
s_sp, s_rxn = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)


# -

for j in range(75, 95):
    display(s_rxn[j])
    print(s_rxn[j].family)
    for h in s_rxn[j].reactants + s_rxn[j].products:
        print(h.thermo.comment)
        
    thermo_comments = [x.thermo.comment for x in s_rxn[j].reactants + s_rxn[j].products]
    if any(['Binding energy corrected by LSR' in x for x in thermo_comments]):
        print('yes')
    print()
    print()

for j in range(35, 45):
    display(e_rxn[j])
    print(e_rxn[j].kinetics)
    for h in e_rxn[j].reactants + e_rxn[j].products:
        print(h.thermo.comment)
    print()
    print()

e_rxn[40].reactants[0].thermo.comment





dir(e_rxn[40].reactants[0])



e_rxn[41].kinetics

s_rxn[114].kinetics



s_rxn[41]

2.72e-9 / .0001

# for i in range(len(e_rxn)):
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# for i in range(5):
for i in range(len(e_rxn)):
    
    s = get_i_thing(e_rxn[i], s_rxn)
    
    
    if s < 0:
        print(f'reaction {i} not found')
        continue
        
    reverse = False
    if not reactions_in_same_direction(e_rxn[i], s_rxn[s]):
        reverse = True
        reverse_reaction = copy.deepcopy(s_rxn[s])
        reverse_reaction.products = s_rxn[s].reactants
        reverse_reaction.reactants = s_rxn[s].products
        reverse_reaction.kinetics = s_rxn[s].generate_reverse_rate_coefficient()
        print('REVERSED')
#     assert 
#     need to generate reverse reaction


    # plot kinetics
    P = 101325
    T = np.linspace(300, 1000, 1001)
    k_e = np.zeros(len(T))
    k_s = np.zeros(len(T))
    for j in range(0, len(T)):  # molm^2
        if reverse:
            k_s[j] = reverse_reaction.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
        else:
            k_s[j] = s_rxn[s].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
        k_e[j] = e_rxn[i].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
        


    # only plot reactions that deviate from emily's
    if np.abs(np.log10(k_s[0]) - np.log10(k_e[0])) < 1.0:
        continue

    
    display(i, s, e_rxn[i])
    plt.plot(1000.0 / T, k_e, label=f'Emily', color=colors[0])
    if reverse:
        plt.plot(1000.0 / T, k_s, label=f'Sevy (reverse)', color=colors[1], linestyle='dashed')
    else:
        plt.plot(1000.0 / T, k_s, label=f'Sevy', color=colors[1], linestyle='dashed')
    
    plt.legend()
    plt.xlabel('1000 K / T')
    plt.ylabel('k')
    plt.yscale('log')
    plt.show()

# +
# check out s_rxn[79]
s = 79
display(s_rxn[s])


reverse_reaction = copy.deepcopy(s_rxn[s])
reverse_reaction.products = s_rxn[s].reactants
reverse_reaction.reactants = s_rxn[s].products
reverse_reaction.kinetics = s_rxn[s].generate_reverse_rate_coefficient()

T = np.linspace(300, 1000, 1001)
k_s = np.zeros(len(T))
k_s_rev = np.zeros(len(T))
for j in range(0, len(T)):  # molm^2
    k_s[j] = s_rxn[s].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_s_rev[j] = reverse_reaction.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)

    
    
    
plt.plot(1000.0 / T, k_s, label=f'Sevy', color=colors[1], linestyle='dashed')
plt.plot(1000.0 / T, k_s_rev, label=f'Sevy (reverse)', color=colors[2], linestyle='dashed')
plt.legend()
plt.xlabel('1000 K / T')
plt.ylabel('k')
plt.yscale('log')
plt.show()
# -

display(e_rxn[60])
print(e_rxn[60].kinetics)

display(s_rxn[82].kinetics)



print(s_rxn[79])

# +
for i in range(len(e_rxn)):
    labels = [x.label for x in e_rxn[i].products + e_rxn[i].reactants]
    if 'CH4X' in labels:
        display(i, e_rxn[i])
        
        
for i in range(len(s_rxn)):
    labels = [x.label for x in s_rxn[i].products + s_rxn[i].reactants]
    if 'CH4X' in labels:
        display(i, s_rxn[i])

# +
# for sp in e_sp:
#     display(sp)
# -



Tref = 1000.0
    T = np.linspace(300, 3000, 1001)
    P = 101325
    k = np.zeros(len(T))
    for j in range(0, len(T)):
        k[j] = my_rxn.get_rate_coefficient(T[j], P)
    plt.plot(1000.0 / T, k, label=f'Aramco', color='black')





# # Load the Surface Mechanism
#

mech_yaml = './cpox_rh_emily/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241028/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241020/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

# +
# mech_yaml = './cpox_pt_20241020/chem_annotated_no79-gas.yaml'
# gas = ct.Solution(mech_yaml, 'gas')
# surf = ct.Interface(mech_yaml, 'surface1', [gas])
# print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
# print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')
# -

mech_yaml = './blondal_pt/blondal_pt.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')





# # Simulate it in a simple reactor

# +
# katrin's script
i_ar = gas.species_index('Ar')
i_ch4 = gas.species_index('CH4(2)')
i_o2 = gas.species_index('O2(3)')
i_co2 = gas.species_index('CO2(11)')
i_h2o = gas.species_index('H2O(10)')
i_h2 = gas.species_index('H2(4)')
i_co = gas.species_index('CO(13)')

# # emily's numbering
# i_ar = gas.species_index('Ar')
# i_ch4 = gas.species_index('CH4(2)')
# i_o2 = gas.species_index('O2(3)')
# i_co2 = gas.species_index('CO2(4)')
# i_h2o = gas.species_index('H2O(5)')
# i_h2 = gas.species_index('H2(6)')
# i_co = gas.species_index('CO(7)')

for idx in [i_ar, i_ch4, i_o2, i_co2, i_h2o, i_h2, i_co]:
    assert idx >= 0
# -

gas.species_names

# +
#######################################################################
# Input Parameters
#######################################################################
CH4_0 = 0.29577464788732394
O2_0 = 0.14788732394366197
Ar_0 = 1.0 - CH4_0 - O2_0


REACTOR_VOLUME = 1.0  # m^3
REACTOR_TEMPERATURE = 273.15  # K
REACTOR_PRESSURE = ct.one_atm
MAX_SIMULATION_TIME = 1.0
CONCENTRATIONS = {
    'CH4(2)': CH4_0,
    'O2(3)': O2_0,
    'Ar': Ar_0,
}


# initialize T and P
gas.TPX = REACTOR_TEMPERATURE, REACTOR_PRESSURE, CONCENTRATIONS
surf.TP = REACTOR_TEMPERATURE, REACTOR_PRESSURE

volume_flow = 1.0

# Catalyst settings (really need to double check these)
catalyst_weight = 4.24e-3
cat_site_per_wt = 5*61.67*1e-6*1e3 # [mol/kg] 1e-6mol/micromole, 1000g/kg
site_density = (
    surf.site_density * 1000
)  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
cat_area = (catalyst_weight * cat_site_per_wt) / site_density  # [m^3]
surf.coverages = "X(1):1.0"


gas_reactor = ct.IdealGasReactor(gas)
gas_reactor.volume = REACTOR_VOLUME
surface_reactor = ct.ReactorSurface(surf, gas_reactor, A=cat_area)

# set up mass flow controllers
inlet = ct.Reservoir(gas)
exhaust = ct.Reservoir(gas)

N = 1001
t_end = 1.0e-5
t = np.linspace(0, t_end, N)
delta_t = t[1] - t[0]

X_cov = np.zeros(len(t))
# CH_cov = np.zeros(len(t))
# O_cov = np.zeros(len(t))
# CO2_cov = np.zeros(len(t))

history = []
gas_history = []

for i in range(0, len(t)):
    X_cov[i] = surf.coverages[surf.species_names.index('X(1)')]
#     CO_cov[i] = surf.coverages[surf.species_names.index('OCX(14)')]
#     O_cov[i] = surf.coverages[surf.species_names.index('OX(8)')]
#     CO2_cov[i] = surf.coverages[surf.species_names.index('CO2X(13)')]
    history.append(surf.coverages)
    gas_history.append(surf.adjacent['gas'].X)
    surf.advance_coverages(delta_t)
    
    
# plt.plot(t, X_cov)
# plt.xlabel('time (s)')

# +
# get ranking of top 10 surface coverages
indices = np.arange(len(history[-1]))
sorted_order = [x for _, x in sorted(zip(history[-1], indices))][::-1]


for i in range(7):
    j = sorted_order[i]
    cov = [x[j] for x in history]
    plt.plot(t, cov, label=surf.species_names[j])
plt.legend()
plt.yscale('log')
# -



# +
# get ranking of top 10 surface coverages
indices = np.arange(len(gas_history[-1]))
sorted_order = [x for _, x in sorted(zip(gas_history[-1], indices))][::-1]


for i in range(7):
    j = sorted_order[i]
    cov = [x[j] for x in gas_history]
    plt.plot(t, cov, label=surf.adjacent['gas'].species_names[j])
plt.legend()
plt.yscale('log')
# -



surf.adjacent['gas'].X

surf.coverages







