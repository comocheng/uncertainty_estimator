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
# notebook to go through thermo estimation in RMG

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


# +

def plot_thermos(thermos, labels=None):
    if type(thermos) != list:
        thermos = [thermos]
    if labels is None:
        labels = ['' for t in thermos]
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(12, 3)
    fig.tight_layout()
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('H (kJ / mol)')
    ax[0].set_title('Enthalpy vs. Temperature')
    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('S (kJ / mol K)')
    ax[1].set_title('Entropy vs. Temperature')
    ax[2].set_xlabel('Temperature (K)')
    ax[2].set_ylabel('Cp (kJ / mol K)')
    ax[2].set_title('Heat Capacity vs. Temperature')
    T = np.linspace(300, 1500, 1001)
    for thermo in thermos:
        H = np.zeros(len(T))
        S = np.zeros(len(T))
        Cp = np.zeros(len(T))
        for i in range(0, len(T)):
            H[i] = thermo.get_enthalpy(T[i]) / 1000.0
            S[i] = thermo.get_entropy(T[i]) / 1000.0
            Cp[i] = thermo.get_heat_capacity(T[i]) / 1000.0
        ax[0].plot(T, H)
        ax[1].plot(T, S)
        ax[2].plot(T, Cp)
    ax[0].legend(labels)
    ax[1].legend(labels)
    ax[2].legend(labels)
    plt.subplots_adjust(wspace=0.25)
    plt.show()


# +
database = rmgpy.data.rmg.RMGDatabase()

thermo_libraries = [
    'surfaceThermoPt111',
    'primaryThermoLibrary',
#     'thermo_DFT_CCSDTF12_BAC',
#     'DFT_QCI_thermo'
]
reaction_libraries = [
    'Surface/CPOX_Pt/Deutschmann2006_adjusted',
    'BurkeH2O2inArHe'
]
database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = thermo_libraries,
    transport_libraries = [],
    reaction_libraries = reaction_libraries,
    seed_mechanisms = [],
    kinetics_families = ['default', 'surface'],
    kinetics_depositories = ['training'],
    depository = False, # Don't bother loading the depository information, as we don't use it
)
# -



# +
# database.thermo.libraries['thermo_DFT_CCSDTF12_BAC'].entries
# -

ethane = rmgpy.species.Species(smiles='C')
ethane.thermo = database.thermo.get_thermo_data(ethane)
display(ethane)
plot_thermos(ethane)
print(ethane.thermo)

sp2 = rmgpy.species.Species(smiles='[CH3]')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

sp2.molecule[0].symmetry_number

sp2 = rmgpy.species.Species(smiles='CC')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

sp2 = rmgpy.species.Species(smiles='CC(C)C')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

sp2 = rmgpy.species.Species(smiles='CCCCCCCC')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

sp2.molecule[0].symmetry_number



sp2 = rmgpy.species.Species(smiles='[CH2]C')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

sp2.thermo.label

sp2.

# +
Ts = np.linspace(300, 1500, 101)
Gs = np.zeros(len(Ts))
Hs = np.zeros(len(Ts))
for i in range(len(Ts)):
    Gs[i] = sp2.get_free_energy(Ts[i])
    Hs[i] = sp2.get_enthalpy(Ts[i]) / 1000.0
    

plt.plot(Ts, Hs*1.05-0.01*Ts, label='DFT')
plt.plot(Ts, Hs*1.05-0.03*Ts, label='DFTB')
plt.plot(Ts, Hs, color='black', label='True Value')
plt.xlabel('T (K)')
plt.ylabel('Group Enthalpy (kJ/mol)')
plt.legend()
ax = plt.gca()

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# -

sp2 = rmgpy.species.Species(smiles='C=CC=C')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

database.thermo.groups['group'].entries['Cds-Cds(Cds-Cds)H'].data

type(database.thermo.groups['group'].entries['C'].data

for item in database.thermo.groups['group'].entries.values():
    print(item.data)

database.thermo.get_thermo_data_from_groups

database.thermo.groups['group'].get_species(sp2)

sp2 = rmgpy.species.Species(smiles='[CH2]O')
sp2.thermo = database.thermo.get_thermo_data(sp2)
display(sp2)
plot_thermos(sp2)
print(sp2.thermo)

sp2.thermo.label

# +
OHX = rmgpy.species.Species().from_adjacency_list(
"""
1 X u0 p0 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}

"""
)
display(OHX)

OHX.thermo = database.thermo.get_thermo_data(OHX)

plot_thermos(OHX)
print(OHX.thermo)

# -



# +
# here we have different alphas
database.kinetics.families['Surface_Adsorption_Dissociative'].rules.entries['Adsorbate;VacantSite1;VacantSite2'][0].data


# -

database.kinetics.families['Surface_Adsorption_Dissociative'].groups.entries['H2O'].item


for entry in database.kinetics.families['Surface_Adsorption_Dissociative'].groups.entries:
    print(database.kinetics.families['Surface_Adsorption_Dissociative'].groups.entries[entry].data)

for entry in database.kinetics.families['Surface_Adsorption_Dissociative'].rules.entries:
    print(database.kinetics.families['Surface_Adsorption_Dissociative'].rules.entries[entry][0].data)











database.kinetics.families['Surface_Adsorption_vdW'].rules.entries['Adsorbate;VacantSite'][0].data

database.kinetics.families['Surface_Adsorption_Single'].rules.entries['Adsorbate;VacantSite'][0].data

for entry in database.kinetics.families['Surface_Adsorption_Single'].groups.entries:
    print(database.kinetics.families['Surface_Adsorption_Single'].groups.entries[entry].data)



OHX.thermo.label

# +
OHX = rmgpy.species.Species().from_adjacency_list(
"""
1 X u0 p0 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}

"""
)
display(OHX)

OHX.thermo = database.thermo.get_thermo_data(OHX)

plot_thermos(OHX)
print(OHX.thermo)

# -



print(database.thermo.libraries['surfaceThermoPt111'].entries['XOH'].item.to_adjacency_list())

# +
OHX = rmgpy.species.Species().from_adjacency_list(
"""
1 X u0 p0 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}

"""
)
display(OHX)

OHX.thermo = database.thermo.get_thermo_data(OHX)

plot_thermos(OHX)
print(OHX.thermo)


# +
CH4X = rmgpy.species.Species().from_adjacency_list(
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 X u0 p0 c0


"""
)
display(CH4X)

CH4X.thermo = database.thermo.get_thermo_data(CH4X)

plot_thermos(CH4X)
print(CH4X.thermo)


# +
CH3X = rmgpy.species.Species().from_adjacency_list(
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 X u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}



"""
)
display(CH3X)

CH3X.thermo = database.thermo.get_thermo_data(CH3X)

plot_thermos(CH3X)
print(CH3X.thermo)

# -

CH3X.molecule[0].get_desorbed_molecules()[0]



# +
sp1 = rmgpy.species.Species().from_adjacency_list(
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3 X u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
"""
)
display(sp1)

sp1.thermo = database.thermo.get_thermo_data(sp1)

plot_thermos(sp1)
print(sp1.thermo)



# -

database.thermo._add_adsorption_correction()

# +
my_mol = CH3X.molecule[0]

labeled_atoms = {'*': my_mol.atoms[4]}
node = database.thermo.groups['adsorptionPt111'].descend_tree(my_mol, labeled_atoms)

print(node)
# -

.des



adsorption_groups.descend_tree(molecule, labeled_atoms)

CH3X.molecule[0].atoms

rmgpy.species.Species(smiles='C[CH2]')

sp1.molecule[0]

database.thermo.groups['adsorptionPt111'].entries['C-*R3'].item

1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 X u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}







database

# # Look at correlations between G and k for BEP

# +
# # ! Reaction index: Chemkin #26; RMG #96
# # ! Template reaction: Surface_Adsorption_Dissociative
# # ! Flux pairs: H2O(5), OHX(31); X(1), HX(21); X(1), HX(21); 
# # ! Estimated using template [O-H;VacantSite1;VacantSite2] for rate rule [H2O;VacantSite1;VacantSite2]
# # ! Euclidian distance = 1.0
# # ! Multiplied by reaction path degeneracy 2.0
# # ! family: Surface_Adsorption_Dissociative
# X(1)+X(1)+H2O(5)<=>HX(21)+OHX(31)                   2.000e-01 0.000     34.487   
#     STICK

# +
# make the reaction
X = rmgpy.species.Species().from_adjacency_list("1 X u0 p0 c0")
H2O = rmgpy.species.Species(smiles='O')
HX = rmgpy.species.Species().from_adjacency_list(
    """
    1 H u0 p0 c0 {2,S}
    2 X u0 p0 c0 {1,S}
    """
)
OHX = rmgpy.species.Species().from_adjacency_list(
    """
    1 H u0 p0 c0 {2,S}
    2 O u0 p2 c0 {1,S} {3,S}
    3 X u0 p0 c0 {2,S}
    """
)

X.thermo = database.thermo.get_thermo_data(X)
H2O.thermo = database.thermo.get_thermo_data(H2O)
HX.thermo = database.thermo.get_thermo_data(HX)
OHX.thermo = database.thermo.get_thermo_data(OHX)


my_reaction = rmgpy.reaction.Reaction()
my_reaction.reactants = [X, X, H2O]
my_reaction.products = [HX, OHX]

display(my_reaction)

# +
# Get family match
fam_rxn_list = database.kinetics.generate_reactions_from_families(
    reactants=my_reaction.reactants,
    products=my_reaction.products,
    only_families=None,
    resonance=True,
)
family = fam_rxn_list[0].family
my_reaction.degeneracy = fam_rxn_list[0].degeneracy


# Estimate the kinetics
database.kinetics.families[family].add_atom_labels_for_reaction(my_reaction)
template_labels = database.kinetics.families[family].get_reaction_template_labels(my_reaction)
template = database.kinetics.families[family].retrieve_template(template_labels)
# node = template[0].label
kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=my_reaction.degeneracy)[0]


# -

print(kinetics)

T = 1000.0
kinetics.to_arrhenius(my_reaction.get_enthalpy_of_reaction(T))

for entry in database.kinetics.families['Surface_Adsorption_Dissociative'].rules.entries:
    print(database.kinetics.families['Surface_Adsorption_Dissociative'].rules.entries[entry][0].data)








