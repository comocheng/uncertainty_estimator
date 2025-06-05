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

# example script to show how to unpack information in pickled files
import pickle
import numpy as np

# +
species_dict_file = 'species_dict.pickle'
reaction_dict_file = 'reaction_dict.pickle'
species_index_file = 'species_index_list.pickle'
reaction_index_file = 'reaction_index_list.pickle'

kinetics_uncertainty_file = 'kinetics_correlation.pickle'
thermo_uncertainty_file = 'thermo_correlation.pickle'


# -

def get_thermo_from_NASA(NASA0, NASA1, T):
    # compute thermo properties from nasa polynomials
    # NASA0 is the lower temperature range and NASA1 is the higher
    # expecting NASA polynomials in the following dictionary format:
#     {'coeffs': array([ 3.53732118e+00, -1.21570202e-03,  5.31615358e-06, -4.89440364e-09,
#          1.45843807e-12, -1.03858843e+03,  4.68368633e+00]),
#      'Tmin': (100,'K'),
#      'Tmax': (1074.56,'K')}
    
    assert T >= NASA0['Tmin']
    assert T <= NASA1['Tmax']
    
    a_low = NASA0['coeffs']
    a_high = NASA1['coeffs']
    
    if T < NASA0['Tmax']:
        cp = a_low[0] + a_low[1] * T + a_low[2] * T**2.0 + a_low[3] * T**3.0 + a_low[4] * T**4.0
        h = a_low[0] * T + a_low[1] / 2.0 * T**2.0 + a_low[2] / 3.0 * T**3.0 + a_low[3] / 4.0 * T**4.0 + a_low[4] / 5.0 * T**5.0 + a_low[5]
        s = a_low[0] * np.log(T) + a_low[1] * T + a_low[2] / 2.0 * T**2.0 + a_low[3] / 3.0 * T**3.0 + a_low[4] / 4.0 * T**4.0 + a_low[6]
    else:
        cp = a_high[0] + a_high[1] * T + a_high[2] * T**2.0 + a_high[3] * T**3.0 + a_high[4] * T**4.0
        h = a_high[0] * T + a_high[1] / 2.0 * T**2.0 + a_high[2] / 3.0 * T**3.0 + a_high[3] / 4.0 * T**4.0 + a_high[4] / 5.0 * T**5.0 + a_high[5]
        s = a_high[0] * np.log(T) + a_high[1] * T + a_high[2] / 2.0 * T**2.0 + a_high[3] / 3.0 * T**3.0 + a_high[4] / 4.0 * T**4.0 + a_high[6]

    R = 8.314472
    cp *= R
    h *= R
    s *= R

    return cp, h, s


# ## Look at the species in the mechanism

# +
# Open the index list - describes the numbering for each species in the mechanism
with open(species_index_file, 'rb') as f:
    species_index_list = pickle.load(f)

# Load the thermo data for each species (just the NASA polynomial info)
with open(species_dict_file, 'rb') as f:
    species_dict = pickle.load(f)
    
# Load the thermo uncertainty correlation matrix
with open(thermo_uncertainty_file, 'rb') as f:
    thermo_uncertainty = pickle.load(f)
# -

# print the species and reaction details
for i in range(len(species_index_list)):
    print('Species index: {}'.format(i))
    print('Species Equation: {}'.format(species_index_list[i]))
    print('Thermo Uncertainty: {}'.format(thermo_uncertainty[i, i]))
    print('Species details: {}'.format(species_dict[species_index_list[i]]))
    
    # get the thermo at 1000 K
    T = 1000
    NASA0 = species_dict[species_index_list[i]]['NASA0']
    NASA1 = species_dict[species_index_list[i]]['NASA1']
    
    Cp, H, S = get_thermo_from_NASA(NASA0, NASA1, T)
    print(f'H @1000K (J/mol): {H}'.format(species_dict[species_index_list[i]]))
    
    print('')



get_thermo_from_NASA(species_dict['OX(8)']['NASA0'], species_dict['OX(8)']['NASA1'], 1000)





species_dict['OCX(14)']['NASA0']['coeffs']

species_dict['OCX(14)']







# ## Look at the reactions in the mechanism

# +
# get the numbering order for the reactions
with open(reaction_index_file, 'rb')as f:
    reaction_index_list = pickle.load(f)
    
# read the reaction dictionary - this has kinetics
with open(reaction_dict_file, 'rb') as f:
    reaction_dict = pickle.load(f)
    
# read kinetics uncertainty matrix pickle file
with open(kinetics_uncertainty_file, 'rb') as f:
    kinetics_uncertainty = pickle.load(f)
# -

for i in range(len(reaction_index_list)):
    print('Reaction index: {}'.format(i))
    print('Reaction Equation: {}'.format(reaction_index_list[i]))
    #print('Reaction Type: {}'.format(reaction_dict[reaction_index_list[i]]['type']))
    print('Kinetics Uncertainty: {}'.format(kinetics_uncertainty[i, i]))
    print('Reaction details: {}'.format(reaction_dict[reaction_index_list[i]]))
    print('')

reaction_dict['OX(8) + OX(8) <=> vacantX(4) + vacantX(4) + O2(3)']

# +
# # confirm that the thermo matches the RMG version
# import rmgpy.chemkin
# species, reaction = rmgpy.chemkin.load_chemkin_file('chem_annotated-gas.inp', dictionary_path='species_dictionary.txt')
# h_error = 0
# cp_error = 0
# s_error = 0
# i = 2
# for T in np.linspace(300, 3000, 20):
    
#     NASA0 = species_dict[species_index_list[i]]['NASA0']
#     NASA1 = species_dict[species_index_list[i]]['NASA1']
#     Cp, H, S = get_thermo_from_NASA(NASA0, NASA1, T)
    
#     h_error += np.float_power(species[i].thermo.get_enthalpy(T) - H, 2.0)
#     cp_error += np.float_power(species[i].thermo.get_heat_capacity(T) - Cp, 2.0)
#     s_error += np.float_power(species[i].thermo.get_entropy(T) - S, 2.0)
# -


