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
import rmgpy.species

import glob

import cycler

import matplotlib.patches
import matplotlib.pyplot as plt
# %matplotlib inline
# -

mech_yaml = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241112/og_lib/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

# +

# Emily's numbering
i_ar = gas.species_index('Ar')
i_ch4 = gas.species_index('CH4(2)')
i_o2 = gas.species_index('O2(3)')
i_co2 = gas.species_index('CO2(4)')
i_h2o = gas.species_index('H2O(5)')
i_h2 = gas.species_index('H2(6)')
i_co = gas.species_index('CO(7)')

for idx in [i_ar, i_ch4, i_o2, i_co2, i_h2o, i_h2, i_co]:
    assert idx >= 0
    
on_catalyst = 1000  # catalyst length 10mm, but it doesn't say where.  let's guess at 1 cm?
off_catalyst = 2000
dt = 1.0
# -

gas_names = gas.species_names
surf_names = surf.species_names

gas_names

dist_array = np.load('../distance_array.npy')
T_array = np.load('../T_array.npy')

# +
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
def plotZoom(gas_out, first=False):
    my_alpha = 0.2
    if first:
        ax = plt.gca()
        catalyst_rect = matplotlib.patches.Rectangle((10, 0), 10, 6, facecolor="gainsboro")
        ax.add_patch(catalyst_rect)
    
#     fig, axs = plt.subplots(1, 1)
#     axs.set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))

#     for i in range(len(gas_out[0, :])):
    for z, i in enumerate([gas_names.index('CH4(2)'), gas_names.index('O2(3)'), gas_names.index('H2(6)'), gas_names.index('CO(7)')]):
        species_name = gas_names[i]
        if species_name.endswith(')'):
            if species_name[-3] == '(':
                species_name = species_name[0:-3]
            else:
                species_name = species_name[0:-4]
        

        if gas_out[:, i].max() > 5.e-3:
            #             print(gas_names[i])
            if first:
                plt.plot(dist_array, gas_out[:, i], label=species_name, alpha=my_alpha, color=colors[z])
            else:
                plt.plot(dist_array, gas_out[:, i], label='_no_label', alpha=my_alpha, color=colors[z])


    axs = plt.gca()
#     axs.set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))
    
    if first:
        axs.plot([dist_array[on_catalyst], dist_array[on_catalyst]], [0, 0.2], linestyle='--', color='xkcd:grey')
        axs.plot([dist_array[off_catalyst], dist_array[off_catalyst]], [0, 0.2], linestyle='--', color='xkcd:grey')
        axs.annotate("catalyst", fontsize=18, xy=(14, 0.07), va='bottom', ha='left')

        leg = axs.legend(loc='lower right', fontsize=18)
        for lh in leg.legendHandles: 
            lh.set_alpha(1)
            lh.set_linewidth(4)
#         axs.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=4)
        axs.set_ylim(0., 0.12)
        axs.set_xlim(8, 25)
        axs.set_xlabel('Distance (mm)', fontsize=22)
        axs.set_ylabel('Flow mol/min', fontsize=22)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        fig = plt.gcf()
        fig.set_figheight(6)
        fig.set_figwidth(12)


# -

results_dir = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/results'
gas_files = glob.glob(os.path.join(results_dir, '*_gas_out.npy'))
surf_files = glob.glob(os.path.join(results_dir, '*_surf_out.npy'))
assert len(gas_files) == len(surf_files)





dir(leg)

# +
gas_out = np.load(gas_files[0])
plotZoom(gas_out, first=True)

for i in range(1, len(gas_files)):
    gas_out = np.load(gas_files[i])
    plotZoom(gas_out)
# -

results_dir = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/perturbed_mechs/results'
gas_files = glob.glob(os.path.join(results_dir, '*_gas_out.npy'))
surf_files = glob.glob(os.path.join(results_dir, '*_surf_out.npy'))
assert len(gas_files) == len(surf_files)

# +
gas_out = np.load(gas_files[0])
plotZoom(gas_out, first=True)

for i in range(1, len(gas_files)):
    gas_out = np.load(gas_files[i])
    plotZoom(gas_out)
# -





# +
results_dir = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/results'
my_file = os.path.join(results_dir, 'gas_rxn_0024_gas_out.npy')
my_file2 = os.path.join(results_dir, 'gas_rxn_0024_surf_out.npy')

gas_out = np.load(my_file)
surf_out = np.load(my_file2)
plotZoom(gas_out)

my_file = os.path.join(results_dir, 'gas_rxn_0025_gas_out.npy')
my_file2 = os.path.join(results_dir, 'gas_rxn_0025_surf_out.npy')

gas_out = np.load(my_file)
surf_out = np.load(my_file2)
plotZoom(gas_out)
# -

plt.plot(gas_out[:, 3])
plt.plot(gas_out[:, 8])
plt.xlim([800, 2500])

# # Get the 10 most sensitive reactions/species

results_dir = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/results'
gas_files = glob.glob(os.path.join(results_dir, '*_gas_out.npy'))
surf_files = glob.glob(os.path.join(results_dir, '*_surf_out.npy'))
assert len(gas_files) == len(surf_files)

gas_files.sort()

'/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/results/gas_spec_0000_gas_out.npy' in gas_files

gas.species()

gas_files[0:10]



# +
# change Ne properties

baseline_gas = np.load('/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/results/gas_spec_0001_gas_out.npy')

baseline_CH4 = baseline_gas[2000, gas_names.index('CH4(2)')]

# +
delta = np.zeros(len(gas_files))

for i in range(len(gas_files)):
    gas_out = np.load(gas_files[i])
    delta[i] = gas_out[2000, gas_names.index('CH4(2)')] - baseline_CH4



# +
abs_delta = np.abs(delta)


indices = np.arange(len(delta))
sorted_order = [x for _, x in sorted(zip(abs_delta, indices))][::-1]
# -

for i in range(10):
    print(gas_files[sorted_order[i]])



# +
sens_spec = [
    rmgpy.species.Species(smiles='[H][H]'),
    rmgpy.species.Species(smiles='[C-]#[O+]'),
    rmgpy.species.Species(smiles='C'),
    rmgpy.species.Species().from_adjacency_list("""1 O u0 p2 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,S} {5,D}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {2,D}"""),
    rmgpy.species.Species(smiles='O=C=O'),
    rmgpy.species.Species().from_adjacency_list("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0"""),
    rmgpy.species.Species(smiles='[O][O]'),
    rmgpy.species.Species().from_adjacency_list("""1 C u0 p0 c0 {2,S} {3,T}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,T}"""),
    rmgpy.species.Species().from_adjacency_list("""1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}"""),

]

for sp in sens_spec:
    display(sp)
# -



gas.species()[7]

gas.species()[8]

gas.species()[3]

surf.species()[14]

gas.species()[5]

surf.species()[12]

gas.species()[4]

surf.species()[8]

surf.reactions()[59]

surf.species()[5]







baseline_gas.shape

dist_array[2000]


