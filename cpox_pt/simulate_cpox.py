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

import cycler

import matplotlib.pyplot as plt
# %matplotlib inline
# -

# load the mechanism
mech_yaml = './cpox_pt_20241020/chem_annotated-gas.yaml'
mech_yaml = './cpox_rh_emily/chem_annotated-gas.yaml'
gas, _ = ct.import_phases(mech_yaml, ["gas", "surface1"])
surf = ct.Interface(mech_yaml, 'surface1')
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')


mech_yaml = './cpox_rh_emily/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241020/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241020/chem_annotated_no79-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241028/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

#mod mech
mech_yaml = './my_mech-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

#mod mech
mech_yaml = './my_mech2-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

#mod mech
mech_yaml = './emily_thermo-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

#mod mech
mech_yaml = './emily_thermo_intersect-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

# #copy of emily's mech
mech_yaml = './emily_copy-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

#emily's reactions, my species
mech_yaml = './emily_kinetics-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241028/chem_annotated_O2-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241031/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241101/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241101/pt55_og_lib_limited_fams/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241101/pt55_og_lib_limited_fams/chem_annotated-gas_O2.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241108/adjusted_lib/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241108/adjusted_lib55/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241108/og_lib/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './blondal_pt/blondal_pt.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_pt_20241107/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241112/vlachos_pt111/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241112/mincat/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241112/c1/cantera/chem_annotated.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')

mech_yaml = './cpox_rh_20241112/handpicked/chem_annotated-gas.yaml'
gas = ct.Solution(mech_yaml, 'gas')
surf = ct.Interface(mech_yaml, 'surface1', [gas])
print(f'This mechanism contains {gas.n_species} gas species and {surf.n_species} surface species')
print(f'This mechanism contains {gas.n_reactions} gas reactions and {surf.n_reactions} surface reactions')



# +
# # katrin's script
# i_ar = gas.species_index('Ar')
# i_ch4 = gas.species_index('CH4(2)')
# i_o2 = gas.species_index('O2(3)')
# i_co2 = gas.species_index('CO2(11)')
# i_h2o = gas.species_index('H2O(10)')
# i_h2 = gas.species_index('H2(4)')
# i_co = gas.species_index('CO(13)')

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
# -





# +
#######################################################################
# Input Parameters
#######################################################################

# unit conversion factors to SI
mm = 0.001
cm = 0.01
ms = mm
minute = 60.0

t_in = 700  # K - in the paper, it was ~698.15K at the start of the cat surface and ~373.15 for the gas inlet temp
t_cat = t_in
length = 70 * mm  # Reactor length- m
diam = 16.5*mm  # Reactor diameter - in m
area = (diam/2.0)**2*np.pi  # Reactor cross section area (area of tube) in m^2
porosity = 0.81  # Monolith channel porosity, from Horn ref 17 sec 2.2.2
cat_area_per_vol = 16000  # I made this up, in m-1. 4500 is lowest that "work" for all base
flow_rate = 4.7  # slpm
flow_rate = flow_rate*.001/60  # m^3/s
tot_flow = 0.208  # from Horn 2007, constant inlet flow rate in mol/min, equivalent to 4.7 slpm
velocity = flow_rate/area  # m/s
# The PFR will be simulated by a chain of 'NReactors' stirred reactors.
NReactors = 7001
# NReactors = 101

on_catalyst = 1000  # catalyst length 10mm, but it doesn't say where.  let's guess at 1 cm?
off_catalyst = 2000
dt = 1.0

reactor_len = length/(NReactors-1)
rvol = area * reactor_len * porosity
# catalyst area in one reactor
cat_area = cat_area_per_vol * rvol
# -



# +

def plotSurf(a):
    gas_out, surf_out, gas_names, surf_names, dist_array, T_array = a

    fig, axs = plt.subplots(1, 2)
    axs[0].set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))

    for i in range(len(gas_out[0, :])):
        if i != i_ar:
            if gas_out[:, i].max() > 5.e-3:
                axs[0].plot(dist_array, gas_out[:, i], label=gas_names[i])
                species_name = gas_names[i]
                if species_name.endswith(')'):
                    if species_name[-3] == '(':
                        species_name = species_name[0:-3]
                    else:
                        species_name = species_name[0:-4]
                if species_name == "O2":
                    axs[0].annotate("O$_2$", fontsize=18,
                                    xy=(dist_array[2200], gas_out[:, i][2200] + gas_out[:, i][2200] / 100.0),
                                    va='bottom', ha='center')
                elif species_name == "CO2":
                    axs[0].annotate("CO$_2$", fontsize=18,
                                    xy=(dist_array[2200], gas_out[:, i][2200] + gas_out[:, i][2200] / 10.0), va='top',
                                    ha='center')
                elif species_name == "CO":
                    axs[0].annotate("CO", fontsize=18, xy=(dist_array[2200], gas_out[:, i][2200] + 0.001),
                                    va='bottom', ha='center')
                elif species_name == "CH2O":
                    axs[0].annotate("CH$_2$O", fontsize=18, xy=(dist_array[2200], gas_out[:, i][2200] + 0.001),
                                    va='bottom', ha='center')
                elif species_name == "CH4":
                    axs[0].annotate("CH$_4$", fontsize=18,
                                    xy=(dist_array[2200], gas_out[:, i][2200] + gas_out[:, i][2200] / 100.0),
                                    va='bottom', ha='center')
                elif species_name == "H2O":
                    axs[0].annotate("H$_2$O", fontsize=18,
                                    xy=(dist_array[2200], gas_out[:, i][2200] + gas_out[:, i][2200] / 40.0), va='top',
                                    ha='center')
                else:
                    axs[0].annotate(species_name, fontsize=18,
                                    xy=(dist_array[-1], gas_out[:, i][-1] + gas_out[:, i][-1] / 10.0), va='top',
                                    ha='center')
            else:
                axs[0].plot(0, 0)

    axs[1].set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))
    # Plot two temperatures (of gas-phase and surface vs only surface.)
    for i in range(len(surf_out[0, :])):
        if surf_out[:, i].max() > 5.e-3:
            axs[1].semilogy(dist_array, surf_out[:, i], label=surf_names[i])
    axs[0].set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))

    axs[0].plot([dist_array[on_catalyst], dist_array[on_catalyst]], [0, 0.2], linestyle='--', color='xkcd:grey')
    axs[0].plot([dist_array[off_catalyst], dist_array[off_catalyst]], [0, 0.2], linestyle='--', color='xkcd:grey')
    axs[0].annotate("catalyst", fontsize=18, xy=(dist_array[on_catalyst], 0.175), va='bottom', ha='left')
    axs[1].plot([dist_array[on_catalyst], dist_array[on_catalyst]], [0, 1.2], linestyle='--', color='xkcd:grey')
    axs[1].plot([dist_array[off_catalyst], dist_array[off_catalyst]], [0, 1.2], linestyle='--', color='xkcd:grey')
    axs[1].annotate("catalyst", fontsize=18, xy=(dist_array[on_catalyst], 1.1), va='bottom', ha='left')

    for item in (
            axs[0].get_xticklabels() + axs[0].get_yticklabels() + axs[1].get_xticklabels() + axs[1].get_yticklabels()):
        item.set_fontsize(18)

    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=2)
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=4)
    axs[0].set_ylim(0., 0.1)
    axs[1].set_ylim(1e-10, 1.2)
    axs[0].set_xlim(5, 25)
    axs[1].set_xlim(9, 21)
    axs[0].set_xlabel('Distance (mm)', fontsize=22)
    axs[1].set_xlabel('Distance (mm)', fontsize=22)
    axs[0].set_ylabel('flow/ mol/min', fontsize=22)
    axs[1].set_ylabel('Site fraction', fontsize=22)
    fig.delaxes(axs[0])  # THIS DELETES THE EXTRA SUBPLOT!

    fig.set_figheight(6)
    fig.set_figwidth(18)

    for n in range(len(gas_names)):
        if gas_names[n] == 'CH4(2)':
            c_in = gas_out[0][n]
        if gas_names[n] == 'O2(3)':
            o_in = gas_out[0][n]
    ratio = c_in / (o_in * 2)
    ratio = round(ratio, 1)

    out_dir = 'figures'
    os.path.exists(out_dir) or os.makedirs(out_dir)
    fig.clf()


# +

def plotZoom(a):
    gas_out, surf_out, gas_names, surf_names, dist_array, T_array = a

    fig, axs = plt.subplots(1, 2)
    axs[0].set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))

    for i in range(len(gas_out[0, :])):
        if i != i_ar:
            if gas_out[:, i].max() > 5.e-3:
                #             print(gas_names[i])
                axs[0].plot(dist_array, gas_out[:, i], label=gas_names[i])
                species_name = gas_names[i]
                if species_name.endswith(')'):
                    if species_name[-3] == '(':
                        species_name = species_name[0:-3]
                    else:
                        species_name = species_name[0:-4]
#                 if species_name == "O2":
#                     axs[0].annotate("O$_2$", fontsize=18, color='y',
#                                     xy=(dist_array[1100], gas_out[:, i][1100] + gas_out[:, i][1100] / 100.0),
#                                     va='bottom', ha='center')
#                 elif species_name == "CO2":
#                     axs[0].annotate("CO$_2$", fontsize=18, color='c',
#                                     xy=(dist_array[2400], gas_out[:, i][2400] + gas_out[:, i][2400] / 10.0), va='bottom',
#                                     ha='center')
#                 elif species_name == "CO":
#                     axs[0].annotate("CO", fontsize=18, color='g', xy=(dist_array[2100], gas_out[:, i][2100] + 0.001),
#                                     va='bottom', ha='center')
#                 elif species_name == "H2":
#                     axs[0].annotate("H$_2$", fontsize=18, color='k', xy=(dist_array[2200], gas_out[:, i][2200] - 0.001),
#                                     va='top', ha='center')
#                 elif species_name == "CH4":
#                     axs[0].annotate("CH$_4$", fontsize=18, color='b',
#                                     xy=(dist_array[1100], gas_out[:, i][1100] + gas_out[:, i][1100] / 100.0),
#                                     va='bottom', ha='center')
#                 elif species_name == "H2O":
#                     axs[0].annotate("H$_2$O", fontsize=18, color='r',
#                                     xy=(dist_array[2100], gas_out[:, i][2100] + gas_out[:, i][2100] / 40.0 + 0.001), va='bottom',
#                                     ha='center')
#                 else:
#                     axs[0].annotate(species_name, fontsize=18,
#                                     xy=(dist_array[-1], gas_out[:, i][-1] + gas_out[:, i][-1] / 10.0), va='top',
#                                     ha='center')
            else:
                axs[0].plot(0, 0)

    axs[1].set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))
    ax2 = axs[0].twinx()
    ax2.plot(dist_array, T_array, label='temperature', color='r', linestyle=':')
    axs[0].set_prop_cycle(cycler.cycler('color', ['m', 'g', 'b', 'y', 'c', 'r', 'k', 'g']))

    axs[0].plot([dist_array[on_catalyst], dist_array[on_catalyst]], [0, 0.2], linestyle='--', color='xkcd:grey')
    axs[0].plot([dist_array[off_catalyst], dist_array[off_catalyst]], [0, 0.2], linestyle='--', color='xkcd:grey')
    axs[0].annotate("catalyst", fontsize=18, xy=(dist_array[on_catalyst], 0.175), va='bottom', ha='left')
    axs[1].plot([dist_array[on_catalyst], dist_array[on_catalyst]], [600.0, 2000], linestyle='--', color='xkcd:grey')
    axs[1].plot([dist_array[off_catalyst], dist_array[off_catalyst]], [600.0, 2000], linestyle='--', color='xkcd:grey')
    axs[1].annotate("catalyst", fontsize=18, xy=(dist_array[on_catalyst], 1800), va='bottom', ha='left')

    for item in (
            axs[0].get_xticklabels() + axs[0].get_yticklabels() + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(18)

    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=4)
    axs[0].set_ylim(0., 0.1)
    axs[1].set_ylim(600.0, 2000)
    axs[0].set_xlim(8, 25)
    axs[1].set_xlim(8, 25)
    axs[0].set_xlabel('Distance (mm)', fontsize=22)
    axs[1].set_xlabel('Distance (mm)', fontsize=22)  # axs[0,1].set_xlabel('time (s)'); axs[1,1].set_xlabel('time (s)')
    axs[0].set_ylabel('flow/ mol/min', fontsize=22)
    ax2.set_ylabel('Temperature (K)', fontsize=22)
    ax2.set_ylim(600, 2000)
    ax2.set_xlim(8, 25)
    fig.delaxes(axs[1])  # THIS DELETES THE EXTRA SUBPLOT!
    fig.set_figheight(6)
    fig.set_figwidth(24)

    for n in range(len(gas_names)):
        if gas_names[n] == 'CH4(2)':
            c_in = gas_out[0][n]
        if gas_names[n] == 'O2(3)':
            o_in = gas_out[0][n]
    ratio = c_in / (o_in * 2)
    ratio = round(ratio, 1)

    out_dir = 'figures'
    os.path.exists(out_dir) or os.makedirs(out_dir)



# +

def monolithFull(gas, surf, temp, mol_in, verbose=False, sens=False):
    """
    Verbose prints out values as you go along
    Sens is for sensitivity, in the form [perturbation, reaction #]
    """
    ch4, o2, ar = mol_in
    ratio = ch4/(2*o2)
    ratio = round(ratio, 1)
    ch4 = str(ch4)
    o2 = str(o2)
    ar = str(ar)
    X = str('CH4(2):' + ch4 + ', O2(3):' + o2 + ', Ar:' + ar)
    if verbose:
        print(X)
    gas.TPX = 273.15, ct.one_atm, X  # need to initialize mass flow rate at STP
    mass_flow_rate = flow_rate * gas.density_mass
    gas.TPX = temp, ct.one_atm, X
    temp_cat = temp
    surf.TP = temp_cat, ct.one_atm
    surf.coverages = 'X(1):1.0'
    gas.set_multiplier(1.0)

    TDY = gas.TDY
    cov = surf.coverages

    if verbose is True:
        print('  distance(mm)   X_CH4        X_O2        X_H2       X_CO       X_H2O       X_CO2')

    # create a new reactor
    gas.TDY = TDY
    r = ct.IdealGasReactor(gas)
    r.volume = rvol

    # create a reservoir to represent the reactor immediately upstream. Note
    # that the gas object is set already to the state of the upstream reactor
    upstream = ct.Reservoir(gas, name='upstream')

    # create a reservoir for the reactor to exhaust into. The composition of
    # this reservoir is irrelevant.
    downstream = ct.Reservoir(gas, name='downstream')

    # Add the reacting surface to the reactor. The area is set to the desired
    # catalyst area in the reactor.
    rsurf = ct.ReactorSurface(surf, r, A=cat_area)

    # The mass flow rate into the reactor will be fixed by using a
    # MassFlowController object.
    # mass_flow_rate = velocity * gas.density_mass * area  # kg/s
    # mass_flow_rate = flow_rate * gas.density_mass
    m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)

    # We need an outlet to the downstream reservoir. This will determine the
    # pressure in the reactor. The value of K will only affect the transient
    # pressure difference.
#     v = ct.PressureController(r, downstream, primary=m, K=1e-5)
    v = ct.PressureController(r, downstream, master=m, K=1e-5)

    sim = ct.ReactorNet([r])
    sim.max_err_test_fails = 12

    # set relative and absolute tolerances on the simulation
    sim.rtol = 1.0e-10
    sim.atol = 1.0e-20

    gas_names = gas.species_names
    surf_names = surf.species_names
    gas_out = []
    surf_out = []
    dist_array = []
    T_array = []

    surf.set_multiplier(0.0)  # no surface reactions until the gauze
    for n in range(NReactors):
        # Set the state of the reservoir to match that of the previous reactor
        gas.TDY = r.thermo.TDY
        upstream.syncState()
        if n == on_catalyst:
            surf.set_multiplier(1.0)
            if sens is not False:
                surf.set_multiplier(1.0 + sens[0], sens[1])
        if n == off_catalyst:
            surf.set_multiplier(0.0)
        sim.reinitialize()
        sim.advance_to_steady_state(max_steps=100000)
        dist = n * reactor_len * 1.0e3  # distance in mm
        dist_array.append(dist)
        T_array.append(surf.T)
        # print "mass_flow_rate", mass_flow_rate,  v.mdot(sim.time), "kg/s"
        kmole_flow_rate = mass_flow_rate / gas.mean_molecular_weight  # kmol/s
        gas_out.append(1000 * 60 * kmole_flow_rate * gas.X.copy())  # molar flow rate in moles/minute
        surf_out.append(surf.X.copy())

        # make reaction diagrams
        out_dir = 'rxnpath'
        os.path.exists(out_dir) or os.makedirs(out_dir)
        elements = ['H', 'O']
        locations_of_interest = [1000, 1150, 1160, 1183, 1196, 1999]
        if sens is False:
            for l in locations_of_interest:
                if n == l:
                    location = str(n / 100)

                    diagram = ct.ReactionPathDiagram(surf, 'X')
                    diagram.title = 'rxn path'
                    diagram.label_threshold = 1e-9
                    dot_file = out_dir + '/rxnpath-' + str(ratio) + '-x-' + location + 'mm.dot'
                    img_file = out_dir + '/rxnpath-' + str(ratio) + '-x-' + location + 'mm.png'
                    img_path = os.path.join(out_dir, img_file)
                    diagram.write_dot(dot_file)
                    os.system('dot {0} -Tpng -o{1} -Gdpi=200'.format(dot_file, img_file))

                    for element in elements:
                        diagram = ct.ReactionPathDiagram(surf, element)
                        diagram.title = element + 'rxn path'
                        diagram.label_threshold = 1e-9
                        dot_file = out_dir + '/rxnpath-' + str(ratio) + '-surf-' + location + 'mm-' + element + '.dot'
                        img_file = out_dir + '/rxnpath-' + str(ratio) + '-surf-' + location + 'mm-' + element + '.png'
                        img_path = os.path.join(out_dir, img_file)
                        diagram.write_dot(dot_file)
                        os.system('dot {0} -Tpng -o{1} -Gdpi=200'.format(dot_file, img_file))
        else:
            pass

        if verbose is True:
            if not n % 100:
                print('  {0:10f}  {1:10f}  {2:10f}  {3:10f} {4:10f} {5:10f} {6:10f}'.format(dist, *gas[
                    'CH4(2)', 'O2(3)', 'H2(6)', 'CO(7)', 'H2O(5)', 'CO2(4)'].X * 1000 * 60 * kmole_flow_rate))

    gas_out = np.array(gas_out)
    surf_out = np.array(surf_out)
    gas_names = np.array(gas_names)
    surf_names = np.array(surf_names)
    data_out = gas_out, surf_out, gas_names, surf_names, dist_array, T_array
    return data_out


ratio = 1.0
fo2 = 1 / (2. * ratio + 1 + 79 / 21)
fch4 = 2 * fo2 * ratio
far = 79 * fo2 / 21
ratio_in = [fch4, fo2, far]

a = monolithFull(gas, surf, t_in, ratio_in, verbose=False)
gas_out, surf_out, gas_names, surf_names, dist_array, T_array = a
# -

mech_yaml

plotSurf(a)

plotZoom(a)

gas_out.shape

surf_out.shape



# i = 8
for i in range(0, 10):
    plt.plot(surf_out[:, i], label=surf.species_names[i])
plt.legend(loc='lower right')
plt.yscale('log')
plt.ylabel('coverage')
plt.show()



len(dist_array)

gas.n_species

gas_out.shape

plt.plot(gas_out[:, 5])


