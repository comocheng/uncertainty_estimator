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
# script to test out changes in the uncertainty tool's sensitivity function

# +
import os
import rmgpy.tools.uncertainty
import rmgpy.chemkin

import numpy as np

import matplotlib.pyplot as plt
# %matplotlib inline

import importlib
importlib.reload(rmgpy.tools.uncertainty)


import rmgpy.rmg.settings
# -



# +
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='test_sens_ct')

sevy_mech = '/home/moon/rmg/RMG-Py/examples/rmg/superminimal/chemkin'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict)

uncertainty.species_list = species_list
uncertainty.reaction_list = reaction_list


# +
T = 1500.0
P = 100000.0
# initial_mole_fractions='H2(1):0.67, O2(2):0.33'
initial_mole_fractions={
    [x for x in uncertainty.species_list if x.to_chemkin() == 'H2(1)'][0]: 0.67,
    [x for x in uncertainty.species_list if x.to_chemkin() == 'O2(2)'][0]: 0.33, 
}


reaction_system_index = 0
sensitive_species = [
    [x for x in uncertainty.species_list if x.to_chemkin() == 'H2(1)'][0],
    [x for x in uncertainty.species_list if x.to_chemkin() == 'O2(2)'][0]
]
termination_time = 1e-0

uncertainty.sensitivity_analysis(
    initial_mole_fractions,
    sensitive_species,
    T,
    P,
    termination_time,
    use_cantera=True,
)

# +
# repeat for non-cantera comparison
T = 1500.0
P = 100000.0
uncertainty2 = rmgpy.tools.uncertainty.Uncertainty(output_directory='test_sens_rmg')

sevy_mech = '/home/moon/rmg/RMG-Py/examples/rmg/superminimal/chemkin'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict)

uncertainty2.species_list = species_list
uncertainty2.reaction_list = reaction_list


initial_mole_fractions2={
    [x for x in uncertainty2.species_list if x.to_chemkin() == 'H2(1)'][0]: 0.67,
    [x for x in uncertainty2.species_list if x.to_chemkin() == 'O2(2)'][0]: 0.33, 
}


sensitive_species2 = [
    [x for x in uncertainty2.species_list if x.to_chemkin() == 'H2(1)'][0],
    [x for x in uncertainty2.species_list if x.to_chemkin() == 'O2(2)'][0]
]

uncertainty2.sensitivity_analysis(
    initial_mole_fractions2,
    sensitive_species2,
    T,
    P,
    termination_time,
    use_cantera=False,
)

# +
# plot the results of the sensitivity (for comparison against the non-cantera way of doing it)


reaction_system_index = 0
for sens_species in sensitive_species:
    csvfile_path1 = os.path.join(uncertainty.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time1, data_list1 = rmgpy.tools.plot.parse_csv_data(csvfile_path1)
    
    csvfile_path2 = os.path.join(uncertainty2.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time2, data_list2 = rmgpy.tools.plot.parse_csv_data(csvfile_path2)
    
    
    times1 = time1.data
    times2 = time2.data
    for data in data_list1:
        if 'dG' in data.label:
            continue
        plt.plot(times1, data.data, label=data.label.split()[-1])

    for data in data_list2:
        if 'dG' in data.label:
            continue
        plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
    plt.title('sens rxns')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
    
    for data in data_list1:
        if 'dln[k' in data.label:
            continue
        plt.plot(times1, data.data, label=data.label.split()[-1])

    for data in data_list2:
        if 'dln[k' in data.label:
            continue
        plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
    plt.title('sens specs')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
# -



# +
# forget sensitivities, just plot the results of the simulations


simfile1 = os.path.join(uncertainty.output_directory, 'solver', f'simulation_1_{len(uncertainty.species_list):d}.csv')
simfile2 = os.path.join(uncertainty2.output_directory, 'solver', f'simulation_1_{len(uncertainty2.species_list):d}.csv')
            
time1, data_list1 = rmgpy.tools.plot.parse_csv_data(simfile1)
time2, data_list2 = rmgpy.tools.plot.parse_csv_data(simfile2)


times1 = time1.data
times2 = time2.data
THRESHOLD = 1e-15

for data in data_list1:
    if data.label.strip().lower() == 'volume':
        v1 = data.data

for data in data_list1:
    if data.label.strip().lower() in ['volume', 'temperature', 'pressure']:
        continue
    if np.sum(data.data) < THRESHOLD:
        continue
#     plt.plot(times1, data.data, label=data.label.split()[-1])
    plt.plot(times1, np.multiply(data.data, v1), label=data.label.split()[-1])

for data in data_list2:
    if data.label.strip().lower() in ['volume', 'temperature', 'pressure']:
        continue
    if np.sum(data.data) < THRESHOLD:
        continue
    plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
plt.title('Simulation Results')
plt.xlabel('time (s)')
plt.ylabel('Concentration (moles)')
# plt.xlim([0, 0.5])
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()


# -



# +
# try a surface mechanism

uncertainty3 = rmgpy.tools.uncertainty.Uncertainty(output_directory='test_sens_ct_surf')

sevy_mech = '/home/moon/nitridation/fe110_20241206/'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_chemkin_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_chemkin_surface)

uncertainty3.species_list = species_list
uncertainty3.reaction_list = reaction_list


# -



# +
# NH3 mech
initial_mole_fractions = {
    [x for x in uncertainty3.species_list if x.to_chemkin() == 'O2(3)'][0]: 0.5,  # O2
    [x for x in uncertainty3.species_list if x.to_chemkin() == 'NH3(2)'][0]: 0.5,  # NH3
}
initial_surface_coverages = {
    [x for x in uncertainty3.species_list if x.to_chemkin() == 'X(1)'][0]: 1.0,  # X
}

sensitive_species = [x for x in uncertainty3.species_list if x.to_chemkin() == 'X(1)'] + \
    [x for x in uncertainty3.species_list if x.to_chemkin() == 'HX(8)']

T = 700
P = 100000.0

termination_time = 1e-3

uncertainty3.sensitivity_analysis(
    initial_mole_fractions,
    sensitive_species,
    T,
    P,
    termination_time,
    use_cantera=True,
    sensitivity_threshold=1e-1,
    initial_surface_coverages=initial_surface_coverages
)

# +
# plot the results of the sensitivity (for comparison against the non-cantera way of doing it)

reaction_system_index = 0
for sens_species in sensitive_species:
    csvfile_path3 = os.path.join(uncertainty3.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time3, data_list3 = rmgpy.tools.plot.parse_csv_data(csvfile_path3)

    times3 = time3.data
    for data in data_list3:
        if 'dG' in data.label:
            continue
        plt.plot(times3, data.data, label=data.label.split()[-1])

    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
# -

# # test out manual species sensitivity on gas-phase

# +
uncertainty4 = rmgpy.tools.uncertainty.Uncertainty(output_directory='manual_sens')

sevy_mech = '/home/moon/rmg/RMG-Py/examples/rmg/superminimal/chemkin'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict)

uncertainty4.species_list = species_list
uncertainty4.reaction_list = reaction_list


T = 1500.0
P = 100000.0
termination_time = 1e-0

sensitive_species = [
    [x for x in uncertainty.species_list if x.to_chemkin() == 'H2(1)'][0],
    [x for x in uncertainty.species_list if x.to_chemkin() == 'O2(2)'][0]
]



initial_mole_fractions4={
    [x for x in uncertainty4.species_list if x.to_chemkin() == 'H2(1)'][0]: 0.67,
    [x for x in uncertainty4.species_list if x.to_chemkin() == 'O2(2)'][0]: 0.33, 
}


sensitive_species4 = [
    [x for x in uncertainty4.species_list if x.to_chemkin() == 'H2(1)'][0],
    [x for x in uncertainty4.species_list if x.to_chemkin() == 'O2(2)'][0]
]

uncertainty4.sensitivity_analysis(
    initial_mole_fractions4,
    sensitive_species4,
    T,
    P,
    termination_time,
    use_cantera=True,
#     sensitivity_threshold=1e-7,
    manual_sens=True
)


# plot the results of the sensitivity (for comparison against the non-cantera way of doing it)


reaction_system_index = 0
for sens_species in sensitive_species4:
    csvfile_path4 = os.path.join(uncertainty4.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time4, data_list4 = rmgpy.tools.plot.parse_csv_data(csvfile_path4)
    
    csvfile_path2 = os.path.join(uncertainty2.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time2, data_list2 = rmgpy.tools.plot.parse_csv_data(csvfile_path2)
    
    
    times4 = time4.data
    times2 = time2.data
    for data in data_list4:
        if 'dG' in data.label:
            continue
        plt.plot(times4, data.data, label=data.label.split()[-1])

    for data in data_list2:
        if 'dG' in data.label:
            continue
        plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
    plt.title('Reactions')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
    
    for data in data_list4:
        if 'dln[k' in data.label:
            continue
#         if 'dG' in data.label:
#             continue
        plt.plot(times4, data.data, label=data.label.split()[-1])

    for data in data_list2:
        if 'dln[k' in data.label:
            continue
#         if 'dG' in data.label:
#             continue
        plt.plot(times2, data.data, label=data.label.split()[-1], linestyle='dashed')
    plt.title('Species')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
# -

# # Test out species sensitivity on surface mechanism

# +
uncertainty5 = rmgpy.tools.uncertainty.Uncertainty(output_directory='manual_sens_surf')

sevy_mech = '/home/moon/nitridation/fe110_20241206/'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_chemkin_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_chemkin_surface)

uncertainty5.species_list = species_list
uncertainty5.reaction_list = reaction_list

# NH3 mech
initial_mole_fractions = {
    [x for x in uncertainty5.species_list if x.to_chemkin() == 'O2(3)'][0]: 0.5,  # O2
    [x for x in uncertainty5.species_list if x.to_chemkin() == 'NH3(2)'][0]: 0.5,  # NH3
}
initial_surface_coverages = {
    [x for x in uncertainty5.species_list if x.to_chemkin() == 'X(1)'][0]: 1.0,  # X
}

sensitive_species5 = [x for x in uncertainty5.species_list if x.to_chemkin() == 'X(1)'] + \
    [x for x in uncertainty5.species_list if x.to_chemkin() == 'HX(8)']

T = 700
P = 100000.0

termination_time = 1e-3

uncertainty5.sensitivity_analysis(
    initial_mole_fractions,
    sensitive_species5,
    T,
    P,
    termination_time,
    use_cantera=True,
#     sensitivity_threshold=1e-1,
    initial_surface_coverages=initial_surface_coverages,
    manual_sens=True
)



# plot the results of the sensitivity (for comparison against the non-cantera way of doing it)

THRESHOLD = 1e-15
reaction_system_index = 0
for sens_species in sensitive_species5:
    csvfile_path5 = os.path.join(uncertainty5.output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time5, data_list5 = rmgpy.tools.plot.parse_csv_data(csvfile_path5)
    

    
    times5 = time5.data
    for data in data_list5:
#         if 'dln[k' in data.label:
#             continue
        if 'dG' in data.label:
            continue
        if np.sum(np.abs(data.data)) < THRESHOLD:
            continue
        plt.plot(times5, data.data, label=data.label.split()[-1])
    plt.title('reactions')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
    
    
    for data in data_list5:
        if 'dln[k' in data.label:
            continue
        if np.sum(np.abs(data.data)) < THRESHOLD:
            continue
#         if 'dG' in data.label:
#             continue

        data.data[data.data == -np.inf] = 0.0

        plt.plot(times5, data.data, label=data.label.split()[-1])
    plt.title('species')
    plt.xlabel('time (s)')
    plt.ylabel('sensitivity')
#     plt.ylim([-1.0, 1.0])
#     plt.xlim([-1e-11, 1e-9])
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
# -



# +
# forget sensitivities, just plot the results of the simulations


simfile5 = os.path.join(uncertainty5.output_directory, 'solver', f'simulation_1_{len(uncertainty5.species_list):d}.csv')
time5, data_list5 = rmgpy.tools.plot.parse_csv_data(simfile5)

times5 = time5.data

THRESHOLD = 1e-15

for data in data_list5:
    if data.label.strip().lower() == 'volume':
        v5 = data.data

for data in data_list5:
    if data.label.strip().lower() in ['volume', 'temperature', 'pressure']:
        continue
    if np.sum(data.data) < THRESHOLD:
        continue
    plt.plot(times5, np.multiply(data.data, v5), label=data.label.split()[-1])
plt.title('Simulation Results')
plt.xlabel('time (s)')
plt.ylabel('Concentration')
# plt.xlim([0, 0.5])
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()

for data in data_list5:
    if data.label.strip().lower() in ['volume', 'temperature', 'pressure']:
        continue
    if np.sum(data.data) < THRESHOLD:
        continue
    if 'X' not in data.label.split()[-1] and 'Pt' not in data.label.split()[-1]:
        continue
    plt.plot(times5, data.data / 2.7200000000000002e-08, label=data.label.split()[-1])
plt.title('Simulation Results')
plt.xlabel('time (s)')
plt.ylabel('Concentration')
# plt.xlim([0, 0.5])
plt.xlim([0, 1e-7])
plt.legend()
# plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()


# -

data_list5[21].data

np.log(np.nan)



[x[len(reaction_list)+2] for x in au]




