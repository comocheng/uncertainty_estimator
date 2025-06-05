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



# -



# pick an annotated chemkin file to analyze
chemkin_file = '/home/moon/uncertainty_estimator/cpox_pt/toy_surface2/chemkin/chem_annotated-gas.inp'
chemkin_file_surface = '/home/moon/uncertainty_estimator/cpox_pt/toy_surface2/chemkin/chem_annotated-surface.inp'
dict_file = '/home/moon/uncertainty_estimator/cpox_pt/toy_surface2/chemkin/species_dictionary.txt'


# +
# make an uncertainty object
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_surf2')
# load the reactions and species (directly from the file)
uncertainty.load_model(chemkin_file, dict_file)

# or indirectly by passing the species/reaction lists through memory
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file, surface_path=chemkin_file_surface)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(
    species_list=species_list,
    reaction_list=reaction_list,
    output_directory='uncertainty_calculations'
)


# +
# Load the database - this takes a while because the averaging up apparently needs Julia
# make sure all of the libraries/families match what was used to generate the mechanism
thermo_libraries = ['surfaceThermoPt111_reduced', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo']
reaction_libraries = ['Surface/CPOX_Pt/Deutschmann2006_adjusted','BurkeH2O2inArHe']
kinetics_families = [
    'default',
    'Surface_Abstraction',
    'Surface_Abstraction_Single_vdW',
    'Surface_Abstraction_vdW',
    #'Surface_Addition_Single_vdW',
    #'Surface_Adsorption_Abstraction_vdW',
    'Surface_Adsorption_Bidentate',
    'Surface_Adsorption_Dissociative',
    'Surface_Adsorption_Dissociative_Double',
    'Surface_Adsorption_Double',
    'Surface_Adsorption_Single',
    'Surface_Adsorption_vdW',
    'Surface_Bidentate_Dissociation',
    'Surface_Dissociation',
    'Surface_Dissociation_Beta',
    'Surface_Dissociation_Double_vdW',
    'Surface_Dissociation_vdW',
    'Surface_Dissociation_to_Bidentate',
    #'Surface_DoubleBond_to_Bidentate',
    #'Surface_Dual_Adsorption_vdW',
    'Surface_EleyRideal_Addition_Multiple_Bond',
    'Surface_Migration',
    'Surface_vdW_to_Bidentate'
]

uncertainty.load_database(
    thermo_libraries=thermo_libraries,
    kinetics_families=kinetics_families,
    reaction_libraries=reaction_libraries,
)

# -

uncertainty.extract_sources_from_model()
uncertainty.compile_all_sources()

uncertainty.assign_parameter_uncertainties(correlated=False)


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


i_methane = get_i_thing(rmgpy.species.Species(smiles='C'), uncertainty.species_list)
i_O2 = get_i_thing(rmgpy.species.Species(smiles='[O][O]'), uncertainty.species_list)

uncertainty.species_list[22]

uncertainty.species_list[23]

np.logspace(-10, -4, 9)

# +
# run sensitivity analysis
correlated = False
initial_mole_fractions={
    uncertainty.species_list[i_methane]: 0.5,  # methane
    uncertainty.species_list[i_O2]: 0.5,  # O2
}

initial_surf_mol_fractions={
    uncertainty.species_list[22]: 1.0  # Pt
}

sensitive_species = [uncertainty.species_list[22], uncertainty.species_list[23]]  # X, #HX
# sensitive_species = []
T = 1500.0  # K
P = 100000.0  # Pa
end_times = np.logspace(-10, -4, 9)
# end_times = [1e-3]

total_var0 = np.zeros(len(end_times))  # species concentrations
total_var1 = np.zeros(len(end_times))

for j, termination_time in enumerate(end_times):

    uncertainty.sensitivity_analysis(
        initial_mole_fractions,
        sensitive_species,
        T,
        P,
        termination_time=termination_time,
        initial_surface_coverages=initial_surf_mol_fractions,
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


# -

std0 = np.float_power(total_var0, 0.5)
std1 = np.float_power(total_var1, 0.5)
mult0 = np.exp(std0)
mult1 = np.exp(std1)


# +
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# -

# %matplotlib inline

times

# +
# plot the simulation
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
csvfile_path = os.path.join(uncertainty.output_directory, 'solver',
                            'simulation_1_51.csv')
time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)
times = time.data
  
for data in data_list:
    if data.label.lower().strip() in ['volume', 'temperature', 'pressure']:
        continue
    if np.median(np.abs(data.data)) < 1e-4:
        if data.label != 'HX(21)':
            continue
    if not 'X' in data.label and not 'Pt' in data.label:
        continue
#     if 'dG' in data.label:
#         continue

    if data.label.split()[-1] == 'X(1)':
        plt.plot(times, data.data, label=data.label.split()[-1])
        closest_indices = []
        for i in range(len(end_times)):
            closest_indices.append(find_nearest_idx(times, end_times[i]))
        
        plt.vlines(times[closest_indices], ymin=np.divide(data.data[closest_indices], mult0), ymax=np.multiply(data.data[closest_indices], mult0), color=colors[0])
        plt.scatter(times[closest_indices], np.multiply(data.data[closest_indices], mult0), marker='_', color=colors[0])
        plt.scatter(times[closest_indices], np.divide(data.data[closest_indices], mult0), marker='_', color=colors[0])
#     if data.label.split()[-1] == 'ethane(1)':
#         plt.plot(times, data.data, label=data.label.split()[-1])
#         closest_indices = []
#         for i in range(len(end_times)):
#             closest_indices.append(find_nearest_idx(times, end_times[i]))
        
#         plt.vlines(times[closest_indices], ymin=np.divide(data.data[closest_indices], mult0), ymax=np.multiply(data.data[closest_indices], mult0), color=colors[0])
#         plt.scatter(times[closest_indices], np.multiply(data.data[closest_indices], mult0), marker='_', color=colors[0])
#         plt.scatter(times[closest_indices], np.divide(data.data[closest_indices], mult0), marker='_', color=colors[0])
#     elif data.label.split()[-1] == '[H][H](11)':
#         plt.plot(times, data.data, label=data.label.split()[-1])
#         closest_indices = []
#         for i in range(len(end_times)):
#             closest_indices.append(find_nearest_idx(times, end_times[i]))
        
#         plt.scatter(times[closest_indices], np.multiply(data.data[closest_indices], mult1), marker='_', color='black')
#         plt.scatter(times[closest_indices], np.divide(data.data[closest_indices], mult1), marker='_', color='black')
    else:
        plt.plot(times, data.data, label=data.label.split()[-1])

    plt.title('Surface Coverages')
    plt.xlabel('time (s)')
    plt.ylabel('Coverage')
    plt.legend()
#     plt.legend(bbox_to_anchor=(1.1, 1.05))
#     plt.yscale('log')
#     plt.ylim([1e-9, 1e1])
    plt.xscale('log')
#     plt.xlim([1e-14, 1e-3])
plt.show()
# -


