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
# script to make perturbed yamls
# -

import os
import sys
import time
import cantera as ct
import numpy as np
import pandas as pd
import concurrent.futures
import rmgpy.chemkin
import subprocess





# +

# get the table index from input for easy parallelization
DELTA_J_MOL = 418.4  # J/mol, but equals 0.1 kcal/mol
R = 8.3144598  # gas constant in J/mol
DELTA = 0.01
chemkin = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/chem_annotated-gas.inp'
surface = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/chem_annotated-surface.inp'

working_dir = os.path.join(os.path.dirname(chemkin))


# +
def perturb_species(species):
    # takes in an RMG species object
    # change the enthalpy offset
    increase = None
    for poly in species.thermo.polynomials:
        new_coeffs = poly.coeffs
        if not increase:
            # Only define the increase in enthalpy once or you'll end up with numerical gaps in continuity
            # increase = DELTA * new_coeffs[5]
            increase = DELTA_J_MOL / R
        new_coeffs[5] += increase
        poly.coeffs = new_coeffs


def perturb_reaction(rxn):
    # takes in an RMG reaction object
    # delta is the ln(k) amount to perturb the A factor
    # delta is a multiplicative factor- units don't matter, yay!
    # does not deepycopy because there's some issues with rmgpy.reactions copying
    if type(rxn.kinetics) == rmgpy.kinetics.chebyshev.Chebyshev:
        rxn.kinetics.coeffs.value_si[0][0] += np.log10(1.0 + DELTA)
    elif type(rxn.kinetics) in [rmgpy.kinetics.falloff.Troe, rmgpy.kinetics.falloff.ThirdBody, rmgpy.kinetics.falloff.Lindemann]:
        if hasattr(rxn.kinetics, 'arrheniusHigh'):
            rxn.kinetics.arrheniusHigh.A.value *= np.exp(DELTA)
        if hasattr(rxn.kinetics, 'arrheniusLow'):
            rxn.kinetics.arrheniusLow.A.value *= np.exp(DELTA)
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.MultiArrhenius:
        for j in range(len(rxn.kinetics.arrhenius)):
            rxn.kinetics.arrhenius[j].A.value *= np.exp(DELTA)
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.PDepArrhenius:
        for j in range(len(rxn.kinetics.arrhenius)):
            if type(rxn.kinetics.arrhenius[j]) == rmgpy.kinetics.arrhenius.Arrhenius:
                rxn.kinetics.arrhenius[j].A.value *= np.exp(DELTA)
            elif type(rxn.kinetics.arrhenius[j]) == rmgpy.kinetics.arrhenius.MultiArrhenius:
                for k in range(len(rxn.kinetics.arrhenius[j].arrhenius)):
                    rxn.kinetics.arrhenius[j].arrhenius[k].A.value *= np.exp(DELTA)
            else:
                raise ValueError(f'weird kinetics {str(rxn.kinetics)}')
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.MultiPDepArrhenius:
        for i in range(len(rxn.kinetics.arrhenius)):
            for j in range(len(rxn.kinetics.arrhenius[i].arrhenius)):
                if type(rxn.kinetics.arrhenius[i].arrhenius[j]) == rmgpy.kinetics.arrhenius.Arrhenius:
                    rxn.kinetics.arrhenius[i].arrhenius[j].A.value *= np.exp(DELTA)
                elif type(rxn.kinetics.arrhenius[i].arrhenius[j]) == rmgpy.kinetics.arrhenius.MultiArrhenius:
                    for k in range(len(rxn.kinetics.arrhenius[i].arrhenius[j].arrhenius)):
                        rxn.kinetics.arrhenius[i].arrhenius[j].arrhenius[k].A.value *= np.exp(DELTA)
                else:
                    raise ValueError(f'weird kinetics {str(rxn.kinetics)}')

    else:  # Arrhenius
        rxn.kinetics.A.value *= np.exp(DELTA)



# +

transport = os.path.join(working_dir, 'tran.dat')
species_dict = os.path.join(working_dir, 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin, dictionary_path=species_dict, transport_path=transport, use_chemkin_names=True)
species_list_surface, reaction_list_surface = rmgpy.chemkin.load_chemkin_file(surface, dictionary_path=species_dict, transport_path=transport, use_chemkin_names=True)

print(f'Loaded {len(species_list)} species, {len(reaction_list)} reactions')
base_yaml_path = os.path.join(working_dir, 'base.yaml')
perturbed_chemkin = os.path.join(working_dir, 'perturbed.inp')
perturbed_chemkin_surface = os.path.join(working_dir, 'perturbed_surface.inp')
perturbed_yaml_path = os.path.join(working_dir, 'perturbed.yaml')
# -

len(species_list) + len(species_list_surface)

len(reaction_list) + len(reaction_list_surface)

# +
# write base cantera
subprocess.run(['ck2yaml', f'--input={chemkin}', f'--transport={transport}', f'--surface={surface}', f'--output={base_yaml_path}'])

for i in range(len(species_list)):
    perturb_species(species_list[i])
for i in range(len(species_list_surface)):
    perturb_species(species_list_surface[i])

for i in range(len(reaction_list)):
    perturb_reaction(reaction_list[i])
for i in range(len(reaction_list_surface)):
    perturb_reaction(reaction_list_surface[i])

# save the results
rmgpy.chemkin.save_chemkin_file(perturbed_chemkin, species_list, reaction_list, verbose=True, check_for_duplicates=True)
rmgpy.chemkin.save_chemkin_surface_file(perturbed_chemkin_surface, species_list_surface, reaction_list_surface, verbose=True, check_for_duplicates=True)



subprocess.run(['ck2yaml', f'--input={perturbed_chemkin}', f'--transport={transport}', f'--surface={perturbed_chemkin_surface}', f'--output={perturbed_yaml_path}'])


# # load the 2 ctis
# base_gas = ct.Solution(base_yaml_path)
# perturbed_gas = ct.Solution(perturbed_yaml_path)
# -



# +

def perturb_species_amount(species, AMOUNT_KCAL_MOL):
    AMOUNT_J_MOL = 4184 * AMOUNT_KCAL_MOL
    # takes in an RMG species object
    # change the enthalpy offset
    increase = None
    for poly in species.thermo.polynomials:
        new_coeffs = poly.coeffs
        if not increase:
            # Only define the increase in enthalpy once or you'll end up with numerical gaps in continuity
            # increase = DELTA * new_coeffs[5] 
            increase = AMOUNT_J_MOL / R
        new_coeffs[5] += increase
        poly.coeffs = new_coeffs


def perturb_reaction_amount(rxn, AMOUNT):
    # takes in an RMG reaction object
    # delta is the ln(k) amount to perturb the A factor
    # delta is a multiplicative factor- units don't matter, yay!
    # does not deepycopy because there's some issues with rmgpy.reactions copying
    if type(rxn.kinetics) == rmgpy.kinetics.chebyshev.Chebyshev:
        rxn.kinetics.coeffs.value_si[0][0] += np.log10(1.0 + AMOUNT)
    elif type(rxn.kinetics) in [rmgpy.kinetics.falloff.Troe, rmgpy.kinetics.falloff.ThirdBody, rmgpy.kinetics.falloff.Lindemann]:
        if hasattr(rxn.kinetics, 'arrheniusHigh'):
            rxn.kinetics.arrheniusHigh.A.value *= np.exp(AMOUNT)
        if hasattr(rxn.kinetics, 'arrheniusLow'):
            rxn.kinetics.arrheniusLow.A.value *= np.exp(AMOUNT)
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.MultiArrhenius:
        for j in range(len(rxn.kinetics.arrhenius)):
            rxn.kinetics.arrhenius[j].A.value *= np.exp(AMOUNT)
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.PDepArrhenius:
        for j in range(len(rxn.kinetics.arrhenius)):
            if type(rxn.kinetics.arrhenius[j]) == rmgpy.kinetics.arrhenius.Arrhenius:
                rxn.kinetics.arrhenius[j].A.value *= np.exp(AMOUNT)
            elif type(rxn.kinetics.arrhenius[j]) == rmgpy.kinetics.arrhenius.MultiArrhenius:
                for k in range(len(rxn.kinetics.arrhenius[j].arrhenius)):
                    rxn.kinetics.arrhenius[j].arrhenius[k].A.value *= np.exp(AMOUNT)
            else:
                raise ValueError(f'weird kinetics {str(rxn.kinetics)}')
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.MultiPDepArrhenius:
        for i in range(len(rxn.kinetics.arrhenius)):
            for j in range(len(rxn.kinetics.arrhenius[i].arrhenius)):
                if type(rxn.kinetics.arrhenius[i].arrhenius[j]) == rmgpy.kinetics.arrhenius.Arrhenius:
                    rxn.kinetics.arrhenius[i].arrhenius[j].A.value *= np.exp(AMOUNT)
                elif type(rxn.kinetics.arrhenius[i].arrhenius[j]) == rmgpy.kinetics.arrhenius.MultiArrhenius:
                    for k in range(len(rxn.kinetics.arrhenius[i].arrhenius[j].arrhenius)):
                        rxn.kinetics.arrhenius[i].arrhenius[j].arrhenius[k].A.value *= np.exp(AMOUNT)
                else:
                    raise ValueError(f'weird kinetics {str(rxn.kinetics)}')

    else:  # Arrhenius
        rxn.kinetics.A.value *= np.exp(AMOUNT)

# -



# +
# Make 100 perturbation files
covariance_matrix = np.array([[2.25      , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 2.25      , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 2.25      , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 2.25      , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 2.25      ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        2.25      , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 2.25      , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 2.25      , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 2.25      , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.11803399]])

np.random.seed(400)

perturbation_mat = np.random.multivariate_normal(np.zeros(10), covariance_matrix, size=1000, check_valid='warn', tol=1e-8)

# -

perturbation_mat[0, :]



# +
mech_dir = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/perturbed_mechs'

perturbed_chemkin_temp = os.path.join(mech_dir, 'perturbed_chemkin_temp.inp')
perturbed_chemkin_surface_temp = os.path.join(mech_dir, 'perturbed_chemkin_surface_temp.inp')

# Must use annotated chemkin file
chemkin_file = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/chem_annotated-gas.inp'
surface_file = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/chem_annotated-surface.inp'
dict_file = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/species_dictionary.txt'

sp_indices = [7, 8, 3, 36, 5, 34, 4, 30, 27]
rxn_index = 105
for i in range(100, 200):
    
    output_file = os.path.join(mech_dir, f'chem_{i:04}.yaml')
    
    total_species_list, total_reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file, surface_path=surface_file)
    
    
    for a, j in enumerate(sp_indices):
        amount = perturbation_mat[i, a]
        perturb_species_amount(total_species_list[j], amount)
    
    perturb_reaction_amount(total_reaction_list[105], perturbation_mat[i, 9])
    
    
    species_list = [x for x in total_species_list if not x.contains_surface_site()]
    species_list_surface = [x for x in total_species_list if x.contains_surface_site()]
    
    reaction_list = [x for x in total_reaction_list if not x.is_surface_reaction()]
    reaction_list_surface = [x for x in total_reaction_list if x.is_surface_reaction()]
    
    
    
    # save the results
    rmgpy.chemkin.save_chemkin_file(perturbed_chemkin_temp, species_list, reaction_list, verbose=True, check_for_duplicates=True)
    rmgpy.chemkin.save_chemkin_surface_file(perturbed_chemkin_surface_temp, species_list_surface, reaction_list_surface, verbose=True, check_for_duplicates=True)



    subprocess.run(['ck2yaml', f'--input={perturbed_chemkin_temp}', f'--transport={transport}', f'--surface={perturbed_chemkin_surface_temp}', f'--output={output_file}'])

    
# # load the 2 ctis
# base_gas = ct.Solution(base_yaml_path)
# perturbed_gas = ct.Solution(perturbed_yaml_path)
# -

species_list[0].contains_surface_site




