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
# a notebook where i can work on developing the code that will export the uncertainty covariance maxtrices

import os
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty
import matplotlib.pyplot as plt

# %matplotlib inline

import importlib
importlib.reload(rmgpy.tools.uncertainty)
importlib.reload(rmgpy.data.thermo)


# -





# +
# pick an annotated chemkin file to analyze
chemkin_file = '/home/moon/uncertainty_estimator/uncertainty_tool_dev/ethane_limit_families/chemkin/chem_annotated.inp'
dict_file = '/home/moon/uncertainty_estimator/uncertainty_tool_dev/ethane_limit_families/chemkin/species_dictionary.txt'


chemkin_file = '/home/moon/rmg/RMG-Py/examples/rmg/minimal/chemkin/chem_annotated.inp'
dict_file = '/home/moon/rmg/RMG-Py/examples/rmg/minimal/chemkin/species_dictionary.txt'

# -

# make an uncertainty object
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations_eth')


# +
# load the reactions and species (directly from the file)
uncertainty.load_model(chemkin_file, dict_file)

# or indirectly by passing the species/reaction lists through memory
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(
    species_list=species_list,
    reaction_list=reaction_list,
    output_directory='uncertainty_calculations'
)

# -

def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


i_ethane = get_i_thing(rmgpy.species.Species(smiles='CC'), uncertainty.species_list)
i_H2 = get_i_thing(rmgpy.species.Species(smiles='[H][H]'), uncertainty.species_list)



# +
# Load the database - this takes a while because the averaging up apparently needs Julia
# make sure all of the libraries/families match what was used to generate the mechanism
# thermo_libraries = ['primaryThermoLibrary', 'BurkeH2O2']
# reaction_libraries = ['BurkeH2O2inN2']
# kinetics_families = [
#     'Disproportionation',
#     'H_Abstraction',
#     'intra_H_migration',
#     'R_Recombination',
#     'Intra_Disproportionation',
# ]


thermo_libraries = ['primaryThermoLibrary']
reaction_libraries = []
kinetics_families = 'default'

uncertainty.load_database(
    thermo_libraries=thermo_libraries,
    kinetics_families=kinetics_families,
    reaction_libraries=reaction_libraries,
)


# +
# Extract sources from model fills out reaction sources dict and the species thermo sources dict
uncertainty.extract_sources_from_model()

# can be accessed like
print(uncertainty.species_sources_dict[uncertainty.species_list[0]])
print()
print(uncertainty.reaction_sources_dict[uncertainty.reaction_list[0]])
# uncertainty.compile_all_sources appears to be completely useless, might want to remove it from RMG

# -

# compile sources helps summarize where all the kinetics/thermo came from
uncertainty.compile_all_sources()
uncertainty.all_kinetic_sources['Rate Rules']['H_Abstraction']

# +
correlated = True

uncertainty.assign_parameter_uncertainties(correlated=correlated)
# -





# +
# run sensitivity analysis

initial_mole_fractions={
    uncertainty.species_list[i_ethane]: 1.0,  # ethane
}
sensitive_species = [uncertainty.species_list[i_ethane], uncertainty.species_list[i_H2]]
# sensitive_species = []
T = 1500.0  # K
P = 100000.0  # Pa
end_times = np.linspace(1e-9, 5.0e-4, 21)
# end_times = [1e-3]

total_var0 = np.zeros(len(end_times))  # species concentrations
total_var1 = np.zeros(len(end_times))

for j, termination_time in enumerate(end_times):

    uncertainty.sensitivity_analysis(
        initial_mole_fractions,
        sensitive_species,
        T,
        P,
        termination_time,
        use_cantera=True,
    )
    
    results = uncertainty.local_analysis(
        sensitive_species,
        reaction_system_index=0,
        correlated=correlated,
        number=5,
        fileformat='.png'
    )

    total_var0[j] = results[sensitive_species[0]][0]
    total_var1[j] = results[sensitive_species[1]][0]



# +
# results[sensitive_species[0]]
# -

std0 = np.float_power(total_var0, 0.5)
std1 = np.float_power(total_var1, 0.5)






# +
# total variance is in var (ln C)

# multiplicative factor is e^var
# std dev = sqrt
mult0 = np.exp(std0)
mult1 = np.exp(std1)
# -





# +
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# +
# plot the simulation
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
csvfile_path = os.path.join(uncertainty.output_directory, 'solver',
                            'simulation_1_26.csv')
time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)
times = time.data
  
for data in data_list:
    if data.label.lower().strip() in ['volume', 'temperature', 'pressure']:
        continue
    if np.median(np.abs(data.data)) < 1e-4:
        continue
#     if 'dG' in data.label:
#         continue

    if data.label.split()[-1] == 'ethane(1)':
        plt.plot(times, data.data, label=data.label.split()[-1])
        closest_indices = []
        for i in range(len(end_times)):
            closest_indices.append(find_nearest_idx(times, end_times[i]))
        
        plt.vlines(times[closest_indices], ymin=np.divide(data.data[closest_indices], mult0), ymax=np.multiply(data.data[closest_indices], mult0), color=colors[0])
        plt.scatter(times[closest_indices], np.multiply(data.data[closest_indices], mult0), marker='_', color=colors[0])
        plt.scatter(times[closest_indices], np.divide(data.data[closest_indices], mult0), marker='_', color=colors[0])
#     elif data.label.split()[-1] == '[H][H](11)':
#         plt.plot(times, data.data, label=data.label.split()[-1])
#         closest_indices = []
#         for i in range(len(end_times)):
#             closest_indices.append(find_nearest_idx(times, end_times[i]))
        
#         plt.scatter(times[closest_indices], np.multiply(data.data[closest_indices], mult1), marker='_', color='black')
#         plt.scatter(times[closest_indices], np.divide(data.data[closest_indices], mult1), marker='_', color='black')
    else:
        plt.plot(times, data.data, label=data.label.split()[-1])

    plt.title('Reactions')
    plt.xlabel('time (s)')
    plt.ylabel('mol / m^3')
    plt.legend()
#     plt.legend(bbox_to_anchor=(1.1, 1.05))
#     plt.yscale('log')
#     plt.ylim([1e-9, 1e1])
plt.show()
# -

uncertainty.output_directory

uncertainty.get_thermo_covariance_matrix()

plt.matshow(uncertainty.thermo_covariance_matrix)
plt.colorbar()
plt.clim([0, 0.1])

uncertainty.database.thermo.groups['radical'].entries

# +
results = uncertainty.local_analysis(
    sensitive_species,
    reaction_system_index=0,
    correlated=False,
    number=9,
    fileformat='.png'
)

print(results[sensitive_species[0]])
# -

print(results[sensitive_species[1]])



# # Parse sensitivity
#

# +
csvfile_path = os.path.join(uncertainty.output_directory, 'solver',
                            'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                 sens_species.index))
time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)
times = time.data
  
for i, data in enumerate(data_list):
    print(data.index, data.label)
# -

uncertainty.reaction_list[59]



# +
results = uncertainty.local_analysis(
    sensitive_species,
    reaction_system_index=0,
    correlated=False,
    number=10,
    fileformat='.png'
)

print(results[sensitive_species[0]])
# total_variance = [resut, reaction_uncertainty, thermo_uncertainty

# +
sens_species = sensitive_species[0]
reaction_system_index = 0
csvfile_path = os.path.join(uncertainty.output_directory, 'solver',
                            'simulation_20{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                 sens_species.index))
time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)
times = time.data
  
for data in data_list:
    if 'dG' in data.label:
        continue
    plt.plot(times, data.data, label=data.label.split()[-1])

    plt.title('Reactions')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
    
#     for data in data_list4:
#         if 'dln[k' in data.label:
#             continue
# #         if 'dG' in data.label:
# #             continue
#         plt.plot(times4, data.data, label=data.label.split()[-1])

#     for data in data_list2:
#         if 'dln[k' in data.label:
#             continue
# #         if 'dG' in data.label:
# #             continue
#         plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
#     plt.title('Species')
#     plt.xlabel('time (s)')
#     plt.ylabel('sensitivity')
#     plt.legend(bbox_to_anchor=(1.1, 1.05))
#     plt.show()

# +
sens_species = sensitive_species[0]
reaction_system_index = 0
csvfile_path = os.path.join(uncertainty.output_directory, 'solver',
                            'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                 sens_species.index))
time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)
times = time.data
  
for data in data_list:
    print(data.label)
#     if 'dG' in data.label:
#         continue
#     plt.plot(times, data.data, label=data.label.split()[-1])

#     plt.title('Reactions')
#     plt.xlabel('time (s)')
#     plt.ylabel('sensitivity')
#     plt.legend(bbox_to_anchor=(1.1, 1.05))
#     plt.show()
    
# #     for data in data_list4:
# #         if 'dln[k' in data.label:
# #             continue
# # #         if 'dG' in data.label:
# # #             continue
# #         plt.plot(times4, data.data, label=data.label.split()[-1])

# #     for data in data_list2:
# #         if 'dln[k' in data.label:
# #             continue
# # #         if 'dG' in data.label:
# # #             continue
# #         plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
# #     plt.title('Species')
# #     plt.xlabel('time (s)')
# #     plt.ylabel('sensitivity')
# #     plt.legend(bbox_to_anchor=(1.1, 1.05))
# #     plt.show()
# -


