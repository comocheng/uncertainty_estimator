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

# Updated script to show how to unpack information in pickled files 1/20/2025
import pickle
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline

# +
# read in the files
# species_dict_file = 'species_dict.pickle'
# reaction_dict_file = 'reaction_dict.pickle'

# thermo_uncertainty_file = 'thermo_cov.npy'
# kinetics_uncertainty_file = 'kinetic_cov.npy'
# overall_uncertainty_file = 'overall_cov.npy'

species_dict_file = '../saved_data_files_v2.2/species_dict.pickle'
reaction_dict_file = '../saved_data_files_v2.2/reaction_dict.pickle'

thermo_uncertainty_file = '../saved_data_files_v2.2/thermo_cov.npy'
kinetics_uncertainty_file = '../saved_data_files_v2.2/kinetic_cov.npy'
overall_uncertainty_file = '../saved_data_files_v2.2/overall_cov.npy'


with open(species_dict_file, 'rb') as f:
    species_info_dict = pickle.load(f)
with open(reaction_dict_file, 'rb') as f:
    reaction_info_dict = pickle.load(f)

thermo_cov = np.load(thermo_uncertainty_file)
kinetics_cov = np.load(kinetics_uncertainty_file)
overall_cov = np.load(overall_uncertainty_file)
# -

len(thermo_cov)

len(species_info_dict)


def get_thermo_from_NASA(NASA0, NASA1, T):
    # compute thermo properties from nasa polynomials  units are Joules and mols
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

print(f'Index\t{"Species Label" + " " * (20 - len("Species Label"))}\tVar(H) (kcal/mol)^2')
for i, spec_key in enumerate(species_info_dict.keys()):
    print(f'{i}\t{spec_key + " " * (20 - len(spec_key))}\t{thermo_cov[i, i]}')

# +
# Example get thermo at 1000K
T = 1000
my_species_label = 'CH2X_1'
NASA0 = species_info_dict[my_species_label]['NASA0']
NASA1 = species_info_dict[my_species_label]['NASA1']

cp, h, s = get_thermo_from_NASA(NASA0, NASA1, T)
print(h)  # in J/mol
# -



# plot the thermo covariance
plt.matshow(thermo_cov)
plt.colorbar()
plt.title('Thermo Covariance (kcal/mol)^2')
# plt.clim([0, 0.01])  # uncomment this to zoom in and see off-diagonals better

# # Look at the Reactions in the File

print(f'Index\t{"Reaction Label" + " " * (40 - len("Reaction Label"))}\tVar(ln k)')
for i, reaction_key in enumerate(reaction_info_dict.keys()):
    print(f'{i}\t{reaction_key + " " * (40 - len(reaction_key))}\t{kinetics_cov[i, i]}')

reaction_info_dict[reaction_key]

# see the info for a given reaction
reaction_info_dict['H2COX_vdW + X_3 <=> H2CO_2X']

# plot the kinetics covariance
plt.matshow(kinetics_cov)
plt.colorbar()
plt.title('Kinetics (ln k) Covariance')
# plt.clim([-0.01, 0.01])  # uncomment this to zoom in and see off-diagonals better

# # plot the overall combined thermo/kinetics covariance matrix

# plot the overall covariance -- species in top left, reactions in bottom right
plt.matshow(overall_cov)
plt.colorbar()
plt.title('Combined Thermo/Kinetic Covariance')
# plt.clim([-0.01, 0.01])  # uncomment this to zoom in and see off-diagonals better




