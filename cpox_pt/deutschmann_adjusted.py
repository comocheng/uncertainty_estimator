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
import cantera as ct
import numpy as np
import os
import copy

import rmgpy.chemkin

import matplotlib.pyplot as plt
# %matplotlib inline



# import os
# import sys
# import pickle
# import copy
# import numpy as np
# import rmgpy.chemkin
# import rmgpy
# import rmgpy.tools.uncertainty
# import rmgpy.kinetics.uncertainties
# import rmgpy.exceptions

# import random

# import rmgpy.kinetics
# import matplotlib.pyplot as plt
# # %matplotlib inline

# sys.path.append('/home/moon/autoscience/reaction_calculator/database')
# import database_fun
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
# -

display(e_rxn[41])
print(e_rxn[41].kinetics)
print(e_rxn[41].generate_reverse_rate_coefficient())

e_rxn[41].get_enthalpy_of_reaction(1000)



# +
adjusted_rev = copy.deepcopy(s_rxn[114])
adjusted_rev.products = s_rxn[114].reactants
adjusted_rev.reactants = s_rxn[114].products
adjusted_rev.kinetics = s_rxn[114].generate_reverse_rate_coefficient()


print(e_rxn[41].get_enthalpy_of_reaction(1000), e_rxn[41].kinetics.Ea.value_si)
print(adjusted_rev.get_enthalpy_of_reaction(1000), adjusted_rev.kinetics.Ea.value_si)

# +
orig_rev = copy.deepcopy(e_rxn[41])
orig_rev.products = e_rxn[41].reactants
orig_rev.reactants = e_rxn[41].products
orig_rev.kinetics = e_rxn[41].generate_reverse_rate_coefficient()


print(orig_rev.get_enthalpy_of_reaction(1000), orig_rev.kinetics.Ea.value_si)
print(s_rxn[114].get_enthalpy_of_reaction(1000), s_rxn[114].kinetics.Ea.value_si)
# -



display(s_rxn[114])
print(s_rxn[114].kinetics)
print(s_rxn[114].generate_reverse_rate_coefficient())

# +
# https://www-sciencedirect-com.ezproxy.neu.edu/science/article/pii/S0920586198003046#BIB6

dupont_pd = rmgpy.kinetics.surface.StickingCoefficient(
    A = 0.4,
    n = 0.0,
    Ea = (0.0, 'kJ/mol')
)

pd_rxn = copy.deepcopy(e_rxn[41])
pd_rxn.kinetics = dupont_pd


# +
# https://www.sciencedirect.com/science/article/pii/S0920586109000285

deutschmann_rh = rmgpy.kinetics.surface.StickingCoefficient(
    A = 0.01,
    n = 0.0,
    Ea = (0.0, 'kJ/mol')
)

rh_rxn = copy.deepcopy(e_rxn[41])
rh_rxn.kinetics = deutschmann_rh

# -



# +
# Plot Emily's reaction and my reverse
display(e_rxn[41])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
   

# plot kinetics
P = 101325
T = np.linspace(300, 1000, 1001)
k_e = np.zeros(len(T))
k_s = np.zeros(len(T))
k_pd = np.zeros(len(T))
k_rh = np.zeros(len(T))
for j in range(0, len(T)):
    k_s[j] = adjusted_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_e[j] = e_rxn[41].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_pd[j] = pd_rxn.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_rh[j] = rh_rxn.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)

plt.plot(1000.0 / T, k_e, label=f'Emily', color=colors[0])
plt.plot(1000.0 / T, k_s, label=f'Sevy (generated from reverse)', color=colors[1], linestyle='dashed')
plt.plot(1000.0 / T, k_pd, label=f'Pd (Dupont)', color=colors[2])
plt.plot(1000.0 / T, k_rh, label=f'Rh (Deutschmann)', color=colors[3])

plt.legend()
plt.xlabel('1000 K / T')
plt.ylabel('k')
plt.yscale('log')
plt.show()
# -

adjusted_rev.kinetics.A.value_si

e_rxn[41].kinetics.A.value_si



# +
# https://www.sciencedirect.com/science/article/pii/S0920586109000285

deutschmann_rh_rev = rmgpy.kinetics.surface.SurfaceArrhenius(
    A = (1.3e22, 'cm^2/(mol*s)'),
    n = 0.0,
    Ea = (355.20, 'kJ/mol')
)

rh_rxn_rev = copy.deepcopy(orig_rev)
rh_rxn_rev.kinetics = deutschmann_rh_rev
# -

orig_rev.kinetics.A

# +
# Plot my reaction and Emily's reverse
display(s_rxn[114])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
   

# plot kinetics
P = 101325
T = np.linspace(300, 1000, 1001)
k_e = np.zeros(len(T))
k_s = np.zeros(len(T))
k_rh = np.zeros(len(T))
for j in range(0, len(T)):
    k_s[j] = s_rxn[114].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_e[j] = orig_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_rh[j] = rh_rxn_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)


plt.plot(1000.0 / T, k_e, label=f'Emily (reverse)', color=colors[0])
plt.plot(1000.0 / T, k_s, label=f'Sevy', color=colors[1], linestyle='dashed')
plt.plot(1000.0 / T, k_rh, label=f'Rh (Deutschmann)', color=colors[3])

plt.legend()
plt.xlabel('1000 K / T')
plt.ylabel('k')
plt.yscale('log')
plt.show()
# -

adjusted_rev.get_enthalpy_of_reaction(1000)

adjusted_rev.kinetics.Ea.value_si

e_rxn[41].kinetics.Ea.value_si

e_rxn[41].reactants[0].get_enthalpy(T)

print(f'Reactants: {np.sum([r.get_enthalpy(T) for r in e_rxn[41].reactants])}')
print(f'Products: {np.sum([p.get_enthalpy(T) for p in e_rxn[41].products])}')

print(f'Reactants: {np.sum([r.get_enthalpy(T) for r in orig_rev.reactants])}')
print(f'Products: {np.sum([p.get_enthalpy(T) for p in orig_rev.products])}')

print(f'Reactants: {np.sum([r.get_enthalpy(T) for r in adjusted_rev.reactants])}')
print(f'Products: {np.sum([p.get_enthalpy(T) for p in adjusted_rev.products])}')

print(f'Reactants: {np.sum([r.get_enthalpy(T) for r in s_rxn[114].reactants])}')
print(f'Products: {np.sum([p.get_enthalpy(T) for p in s_rxn[114].products])}')



# plot the enthalpy
T = 700
display(e_rxn[41])
plt.plot([1, 2, 3], [0, e_rxn[41].kinetics.Ea.value_si, e_rxn[41].get_enthalpy_of_reaction(T)], marker='_', color=colors[0])
plt.plot([1, 2, 3], [0, adjusted_rev.kinetics.Ea.value_si, adjusted_rev.get_enthalpy_of_reaction(T)], marker='_', color=colors[1])


# plot the enthalpy
T = 700
display(s_rxn[114])
plt.plot([1, 2, 3], [0, orig_rev.kinetics.Ea.value_si, orig_rev.get_enthalpy_of_reaction(T)], marker='_', color=colors[0], label='emily')
plt.plot([1, 2, 3], [0, s_rxn[114].kinetics.Ea.value_si, s_rxn[114].get_enthalpy_of_reaction(T)], marker='_', color=colors[1], label='sevy')
plt.plot([1, 2, 3], [0, rh_rxn_rev.kinetics.Ea.value_si, rh_rxn_rev.get_enthalpy_of_reaction(T)], marker='_', color=colors[3], label='Rh deutschmann')
plt.legend()
plt.ylabel('H (kJ/mol)')
plt.xlabel('Reaction coordinate')

rh_rxn_rev.get_enthalpy_of_reaction(T)

rh_rxn_rev.kinetics.Ea.value_si

s_rxn[114]

s_rxn[114].kinetics.Ea.value_si



s_rxn[114].kinetics.Ea.value_si

s_rxn[114].get_enthalpy_of_reaction(1000)

# # what if I raise the barrier to the endothermicity?

# +
my_rxn = copy.deepcopy(s_rxn[114])
display(my_rxn)
my_kinetics = rmgpy.kinetics.surface.SurfaceArrhenius(
    A = (3.7e21, 'cm^2/(mol*s)'),
    n = 0.0,
    Ea = (412.6, 'kJ/mol') # Ea = (278.7, 'kJ/mol')  # raise barrier to enthalpy of reaction
)

my_rxn.kinetics = my_kinetics


my_rxn_rev = 

# +
# Plot my reaction and Emily's reverse
display(s_rxn[114])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
   

# plot kinetics
P = 101325
T = np.linspace(300, 1000, 1001)
k_e = np.zeros(len(T))
k_s = np.zeros(len(T))
k_rh = np.zeros(len(T))
k_mod = np.zeros(len(T))
for j in range(0, len(T)):
    k_s[j] = s_rxn[114].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_e[j] = orig_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_rh[j] = rh_rxn_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_mod[j] = my_rxn.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)

plt.plot(1000.0 / T, k_e, label=f'Emily (reverse)', color=colors[0])
plt.plot(1000.0 / T, k_s, label=f'Sevy', color=colors[1], linestyle='dashed')
plt.plot(1000.0 / T, k_rh, label=f'Rh (Deutschmann)', color=colors[3])
plt.plot(1000.0 / T, k_mod, label=f'Sevy (raised barrier)', color=colors[4])

plt.legend()
plt.xlabel('1000 K / T')
plt.ylabel('k')
plt.yscale('log')
plt.show()
# -

my_rxn_rev = copy.deepcopy(my_rxn)
my_rxn_rev.products = my_rxn.reactants
my_rxn_rev.reactants = my_rxn.products
my_rxn_rev.kinetics = my_rxn.generate_reverse_rate_coefficient()


# +
# Plot Emily's reaction and my reverse
display(e_rxn[41])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
   

# plot kinetics
P = 101325
T = np.linspace(300, 1000, 1001)
k_e = np.zeros(len(T))
k_s = np.zeros(len(T))
k_pd = np.zeros(len(T))
k_rh = np.zeros(len(T))
k_mod_rev = np.zeros(len(T))
for j in range(0, len(T)):
    k_s[j] = adjusted_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_e[j] = e_rxn[41].get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_pd[j] = pd_rxn.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_rh[j] = rh_rxn.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)
    k_mod_rev[j] = my_rxn_rev.get_rate_coefficient(T[j], P, surface_site_density=2.72E-5)

plt.plot(1000.0 / T, k_e, label=f'Emily', color=colors[0])
plt.plot(1000.0 / T, k_s, label=f'Sevy (reverse)', color=colors[1], linestyle='dashed')
plt.plot(1000.0 / T, k_pd, label=f'Pd (Dupont)', color=colors[2])
plt.plot(1000.0 / T, k_rh, label=f'Rh (Deutschmann)', color=colors[3])
plt.plot(1000.0 / T, k_mod_rev, label=f'Sevy (reverse) (raised barrier)', color=colors[4])

plt.legend()
plt.xlabel('1000 K / T')
plt.ylabel('k')
plt.yscale('log')
plt.show()

# +
my_rxn = copy.deepcopy(s_rxn[114])

display(my_rxn)
# -

my_rxn.products[0]

my_rxn.products[0].thermo

my_rxn.reactants[1].thermo.E0

for i in range(len(s_rxn)):
    try:
        s_rxn[i].fix_barrier()
        print('it worked!', i)
    except AttributeError:
        pass

my_rxn.products[0].get_thermo_data()

my_rxn.kinetics



[spec.get_thermo_data().E0.value_si if spec.get_thermo_data().E0 is not None else spec.get_thermo_data().to_wilhoit().E0.value_si for spec in my_rxn.products]



my_rxn.fix_barrier_height()









# +
# what if i try estimating this with families?



database = rmgpy.data.rmg.RMGDatabase()


thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo']
reactionLibraries = ['Surface/CPOX_Pt/Deutschmann2006_adjusted', 'BurkeH2O2inArHe']


database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = thermoLibraries,
    transport_libraries = [],
    reaction_libraries = reactionLibraries,
    seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
    kinetics_families = ['surface', 'default'],
    kinetics_depositories = ['training'],
    #frequenciesLibraries = self.statmechLibraries,
    depository = False, # Don't bother loading the depository information, as we don't use it
)
# -







print(e_rxn[41].reactants[2].to_adjacency_list())

display(e_rxn[41])
my_rxn = copy.deepcopy(e_rxn[41])
fam_rxn_list = database.kinetics.generate_reactions_from_families(
    reactants=my_rxn.reactants,
    products=my_rxn.products,
    only_families=None,
    resonance=True,
)

fam_rxn_list





database.kinetics.families[family].add_atom_labels_for_reaction(my_rxn)

database.kinetics.families







from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
import rmgpy.constants as constants
from rmgpy.reaction import Reaction
from rmgpy.species import Species, TransitionState
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech.rotation import NonlinearRotor
from rmgpy.statmech.torsion import HinderedRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.thermo import Wilhoit, ThermoData, NASA, NASAPolynomial
from rmgpy.kinetics import (
    Arrhenius,
    ArrheniusEP,
    MultiArrhenius,
    PDepArrhenius,
    MultiPDepArrhenius,
    ThirdBody,
    Troe,
    Lindemann,
    Chebyshev,
    SurfaceArrhenius,
    StickingCoefficient,
)


# +

class TestReaction:
    """
    Contains unit tests of the Reaction class.
    """

    def setup_class(self):
        """
        A method that is called prior to each unit test in this class.
        """
        ethylene = Species(
            label="C2H4",
            conformer=Conformer(
                E0=(44.7127, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(28.0313, "amu"),
                    ),
                    NonlinearRotor(
                        inertia=(
                            [3.41526, 16.6498, 20.065],
                            "amu*angstrom^2",
                        ),
                        symmetry=4,
                    ),
                    HarmonicOscillator(
                        frequencies=(
                            [
                                828.397,
                                970.652,
                                977.223,
                                1052.93,
                                1233.55,
                                1367.56,
                                1465.09,
                                1672.25,
                                3098.46,
                                3111.7,
                                3165.79,
                                3193.54,
                            ],
                            "cm^-1",
                        ),
                    ),
                ],
                spin_multiplicity=1,
                optical_isomers=1,
            ),
        )

        hydrogen = Species(
            label="H",
            conformer=Conformer(
                E0=(211.794, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(1.00783, "amu"),
                    ),
                ],
                spin_multiplicity=2,
                optical_isomers=1,
            ),
        )

        ethyl = Species(
            label="C2H5",
            conformer=Conformer(
                E0=(111.603, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(29.0391, "amu"),
                    ),
                    NonlinearRotor(
                        inertia=(
                            [4.8709, 22.2353, 23.9925],
                            "amu*angstrom^2",
                        ),
                        symmetry=1,
                    ),
                    HarmonicOscillator(
                        frequencies=(
                            [
                                482.224,
                                791.876,
                                974.355,
                                1051.48,
                                1183.21,
                                1361.36,
                                1448.65,
                                1455.07,
                                1465.48,
                                2688.22,
                                2954.51,
                                3033.39,
                                3101.54,
                                3204.73,
                            ],
                            "cm^-1",
                        ),
                    ),
                    HinderedRotor(
                        inertia=(1.11481, "amu*angstrom^2"),
                        symmetry=6,
                        barrier=(0.244029, "kJ/mol"),
                        semiclassical=None,
                    ),
                ],
                spin_multiplicity=2,
                optical_isomers=1,
            ),
        )

        TS = TransitionState(
            label="TS",
            conformer=Conformer(
                E0=(266.694, "kJ/mol"),
                modes=[
                    IdealGasTranslation(
                        mass=(29.0391, "amu"),
                    ),
                    NonlinearRotor(
                        inertia=(
                            [6.78512, 22.1437, 22.2114],
                            "amu*angstrom^2",
                        ),
                        symmetry=1,
                    ),
                    HarmonicOscillator(
                        frequencies=(
                            [
                                412.75,
                                415.206,
                                821.495,
                                924.44,
                                982.714,
                                1024.16,
                                1224.21,
                                1326.36,
                                1455.06,
                                1600.35,
                                3101.46,
                                3110.55,
                                3175.34,
                                3201.88,
                            ],
                            "cm^-1",
                        ),
                    ),
                ],
                spin_multiplicity=2,
                optical_isomers=1,
            ),
            frequency=(-750.232, "cm^-1"),
        )

        self.reaction = Reaction(
            reactants=[hydrogen, ethylene],
            products=[ethyl],
            kinetics=Arrhenius(
                A=(501366000.0, "cm^3/(mol*s)"),
                n=1.637,
                Ea=(4.32508, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(300, "K"),
                Tmax=(2500, "K"),
            ),
            transition_state=TS,
            degeneracy=2,
        )
        self.reaction.kinetics.comment = """
        Multiplied by reaction path degeneracy 2.0
        """

        # CC(=O)O[O]
        acetylperoxy = Species(
            label="acetylperoxy",
            molecule=[Molecule(smiles="CC(=O)O[O]")],
            thermo=Wilhoit(
                Cp0=(4.0 * constants.R, "J/(mol*K)"),
                CpInf=(21.0 * constants.R, "J/(mol*K)"),
                a0=-3.95,
                a1=9.26,
                a2=-15.6,
                a3=8.55,
                B=(500.0, "K"),
                H0=(-6.151e04, "J/mol"),
                S0=(-790.2, "J/(mol*K)"),
            ),
        )

        # C[C]=O
        acetyl = Species(
            label="acetyl",
            molecule=[Molecule(smiles="C[C]=O")],
            thermo=Wilhoit(
                Cp0=(4.0 * constants.R, "J/(mol*K)"),
                CpInf=(15.5 * constants.R, "J/(mol*K)"),
                a0=0.2541,
                a1=-0.4712,
                a2=-4.434,
                a3=2.25,
                B=(500.0, "K"),
                H0=(-1.439e05, "J/mol"),
                S0=(-524.6, "J/(mol*K)"),
            ),
        )

        # [O][O]
        oxygen = Species(
            label="oxygen",
            molecule=[Molecule(smiles="[O][O]")],
            thermo=Wilhoit(
                Cp0=(3.5 * constants.R, "J/(mol*K)"),
                CpInf=(4.5 * constants.R, "J/(mol*K)"),
                a0=-0.9324,
                a1=26.18,
                a2=-70.47,
                a3=44.12,
                B=(500.0, "K"),
                H0=(1.453e04, "J/mol"),
                S0=(-12.19, "J/(mol*K)"),
            ),
        )

        self.reaction2 = Reaction(
            reactants=[acetyl, oxygen],
            products=[acetylperoxy],
            kinetics=Arrhenius(
                A=(2.65e12, "cm^3/(mol*s)"),
                n=0.0,
                Ea=(0.0, "kJ/mol"),
                T0=(1, "K"),
                Tmin=(300, "K"),
                Tmax=(2000, "K"),
            ),
        )

        oxygen_atom = Species().from_smiles("[O]")
        so2 = Species().from_smiles("O=S=O")
        so3 = Species().from_smiles("O=S(=O)=O")

        self.reaction3 = Reaction(
            reactants=[oxygen_atom, so2],
            products=[so3],
            kinetics=Arrhenius(A=(3.7e11, "cm^3/(mol*s)"), n=0, Ea=(1689, "cal/mol"), T0=(1, "K")),
        )

        H2 = Species().from_smiles("[H][H]")
        PO3 = Species().from_smiles("[O]P(=O)=O")
        HOPO2 = Species().from_smiles("OP(=O)=O")
        H_atom = Species().from_smiles("[H]")

        self.reaction4 = Reaction(
            reactants=[H2, PO3],
            products=[HOPO2, H_atom],
            kinetics=Arrhenius(A=(2.4e7, "cm^3/(mol*s)"), n=1.38, Ea=(15.38, "kcal/mol"), T0=(1, "K")),
        )
        self.reaction4_pairs = [(PO3, HOPO2), (H2, H_atom)]


# -

my_test = TestReaction()

my_test.setup_class()

my_test.reaction2.fix_barrier_height()

my_test.reaction2.reactants[0].thermo.E0

my_test.reaction2.reactants[0].thermo

my_rxn.reactants[0].thermo.get_entropy(298)

my_rxn.reactants[0].get_thermo_data().to_wilhoit()

my_rxn.reactants[0].thermo.to_wilhoit()



my_rxn.products[2].thermo.to_wilhoit()

my_sp = rmgpy.species.Species(
)
