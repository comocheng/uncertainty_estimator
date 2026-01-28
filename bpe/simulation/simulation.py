"""Module for running simulations of the Horn CPOX experiment"""

import cantera as ct
import numpy as np
import signal
import os
import copy
import scipy.interpolate
import logging

import pandas as pd

import matplotlib.pyplot as plt

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

# constants/objects to set up once at the beginning and then reuse for each simulation
# these are the same for all simulations and do not need to be reloaded each time
N_REACTORS = 1001
MIN_SIM_DIST = -0.004  # in meters
MAX_SIM_DIST = 0.015  # in meters
TOTAL_PFR_LEN = MAX_SIM_DIST - MIN_SIM_DIST
INDIVIDUAL_CSTR_LEN = TOTAL_PFR_LEN / N_REACTORS
dist_array = np.linspace(MIN_SIM_DIST, MAX_SIM_DIST, N_REACTORS)

TIMEOUT_SECONDS = 30  # seconds to wait before timing out a simulation

UNCERTAINTY_REPO = os.environ['UNCERTAINTY_REPO']
# set up the temperature profile as a function of distance along the reactor
# pt_data = os.path.join(os.path.dirname(__file__), '../../cpox_pt/horn_data/pt_profiles_smooth.csv')
pt_data = os.path.join(UNCERTAINTY_REPO, './cpox_pt/horn_data/pt_profiles_smooth.csv')
df = pd.read_csv(pt_data)
distances = (df['Distance (mm)'] - 10.0) / 1000.0  # ignore the 10mm of no/catalyst space
exp_Ts = df['Temperature (K)']
temperature_function = scipy.interpolate.interp1d(distances, exp_Ts, fill_value='extrapolate')
temperature_profile = temperature_function(dist_array)
use_temperature_profile = True

# set up the catalyst on/off indices
CAT_ON_INDEX = np.argmin(np.abs(dist_array - 0))
CAT_OFF_INDEX = np.argmin(np.abs(dist_array - 0.01))

# ----------------------- Parameters specific to Horn CPOX experiment -----------------------
REACTOR_DIAMETER = 0.0165
CROSS_SECTION_AREA = (REACTOR_DIAMETER / 2.0) ** 2.0 * np.pi
POROSITY = 0.81  # Monolith channel porosity, from Horn ref 17 sec 2.2.2
CAT_AREA_PER_VOL = 16000  # made-up
FLOW_RATE_SLPM = 4.7  # slpm
FLOW_RATE = FLOW_RATE_SLPM * 0.001 / 60  # m^3/s
velocity = FLOW_RATE / CROSS_SECTION_AREA  # m/s
residence_time = INDIVIDUAL_CSTR_LEN / velocity # unit in s
individual_cstr_vol = CROSS_SECTION_AREA * INDIVIDUAL_CSTR_LEN * POROSITY
individual_cstr_cat_area = CAT_AREA_PER_VOL * individual_cstr_vol
T_INLET = 700.0  # K  -- won't matter if using the temperature profile
x_CH4 = 0.296  # inlet mole fraction of methane
x_O2 = 0.147
x_Ar = 1.0 - x_CH4 - x_O2
P_INLET = ct.one_atm  # Pa

# Wall settings, but we'll probably ignore in favor of the temperature profile
wall_on = False
wall_area = INDIVIDUAL_CSTR_LEN * np.pi * REACTOR_DIAMETER
U_hx = 100.0 # heat transfer rate for the wall
env_temp = 300.0  # K, temperature of the environment for the wall boundary condition

# Plotting settings
FIG_HEIGHT = 6.0
FIG_WIDTH = 12.0


class TimeoutException(Exception):   # Custom exception class
    pass
def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


# TODO add handling to make sure all the colors are consistent across species

def get_i_thing(ref_composition, phase):
    """Helper function for getting the index of a species in a Cantera phase given its composition"""
    for i in range(phase.n_species):
        if phase.species()[i].composition == ref_composition:
            return i
    assert False, f"Could not find species with composition {ref_composition} in phase {phase.name}"

def increase_enthalpy(phase, species_index, increase_enthalpy_J_per_mol):
    """Helper function for increasing the enthalpy of a species in a Cantera phase by a specified amount (in J/mol)"""
    data_copy = copy.deepcopy(phase.species()[species_index].input_data)
    for i in range(len(data_copy['thermo']['data'])):
        data_copy['thermo']['data'][i][5] += increase_enthalpy_J_per_mol * 1000.0 / ct.gas_constant  # Cantera needs J / kmol
    new_sp = ct.Species().from_dict(data_copy)
    phase.modify_species(species_index, new_sp)


def run_simulation(
    mech_yaml=None,
    gas=None,
    surf=None,
    surf_thermo_perturb=None,
    surf_kinetics_perturb=None,
):
    """Run a simulation of the Horn CPOX experiment

    Parameters
    ----------
    mech_yaml : str
        Path to the Cantera mechanism file
    gas : cantera.Solution
        Cantera gas object
    surf : cantera.Interface
        Cantera surface object
    surf_thermo_perturb : dict
        key is index of the species in the surface phase, value is the amount to perturb the enthalpy in J/kmol
    surf_kinetics_perturb : dict
        key is index of the reaction in the surface phase, value is the multiplier to apply to the reaction rate
        
    Returns
    -------
    output gas concentration
         
    """

    # accept mechanism file or the gas and surf objects directly
    if not gas is None and not surf is None:
        assert isinstance(gas, ct.Solution), "gas must be a Cantera Solution object"
        assert isinstance(surf, ct.Interface), "surf must be a Cantera Interface object"
        logging.debug("Using provided gas and surface objects for simulation")
    elif mech_yaml is not None:
        gas = ct.Solution(mech_yaml)
        surf = ct.Interface(mech_yaml, "surface1", [gas])
        logging.debug(f"Loaded mechanism file: {mech_yaml}")
    else:
        raise ValueError("Either gas and surf objects or a mechanism file must be provided")
    

    # Apply perturbations to the surface thermodynamics
    if surf_thermo_perturb:
        for index in surf_thermo_perturb.keys():
            increase_enthalpy(surf, index, surf_thermo_perturb[index])
            logging.debug(f"Increased enthalpy of surface species {surf.species_names[index]} by {surf_thermo_perturb[index]} J/mol")
    
    # Get indices of key species
    i_Ar = get_i_thing({'Ar': 1.0}, gas)
    i_CH4 = get_i_thing({'C': 1.0, 'H': 4.0}, gas)
    i_O2 = get_i_thing({'O': 2.0}, gas)
    i_X = get_i_thing({'X': 1.0}, surf)
    
    # Initialize the reactor
    X = f'{gas.species_names[i_CH4]}: {x_CH4}, {gas.species_names[i_O2]}: {x_O2}, {gas.species_names[i_Ar]}: {x_Ar}'
    gas.TPX = 273.15, ct.one_atm, X  # need to initialize mass flow rate at STP
    mass_flow_rate = FLOW_RATE * gas.density_mass
    gas.TPX = T_INLET, P_INLET, X
    surf.TP = T_INLET, P_INLET

    initial_coverages = np.zeros_like(surf.coverages)
    initial_coverages[i_X] = 1.0
    surf.coverages = initial_coverages

    
    # create a new reactor
    if use_temperature_profile:
        r = ct.IdealGasReactor(gas, energy='off')
    else:
        r = ct.IdealGasReactor(gas, energy='on')
    r.volume = individual_cstr_vol
    upstream = ct.Reservoir(gas, name='upstream')
    downstream = ct.Reservoir(gas, name='downstream')
    rsurf = ct.ReactorSurface(surf, r, A=individual_cstr_cat_area)
    m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)
    if ct.__version__ in ['2.6.0']:
        v = ct.PressureController(r, downstream, master=m, K=1e-5)
    else:
        v = ct.PressureController(r, downstream, primary=m, K=1e-5)

    if wall_on and not use_temperature_profile:
        air = ct.Solution("air.yaml")
        air.TP = env_temp, ct.one_atm
        env = ct.Reservoir(air)

        # add a wall to the reactor
        w = ct.Wall(r, env, A=wall_area, U=U_hx)

    sim = ct.ReactorNet([r])
    sim.max_err_test_fails = 12

    # we already know the shape of the output arrays, so we should preallocate them for speed
    gas_out = np.zeros((N_REACTORS, gas.n_species)) + np.nan
    # Don't need to store these other outputs, but we could
    surf_out = np.zeros((N_REACTORS, surf.n_species)) + np.nan
    gas_rates = np.zeros((N_REACTORS, gas.n_reactions)) + np.nan
    surf_rates = np.zeros((N_REACTORS, surf.n_reactions)) + np.nan
    # T_array = np.zeros(N_REACTORS) + np.nan


    # start in the non-catalyst region
    surf.set_multiplier(0.0)
    
    for n in range(N_REACTORS):
        signal.alarm(TIMEOUT_SECONDS)  # Set the alarm for the timeout

        # Set the state of the reservoir to match that of the previous reactor
        if use_temperature_profile:
            gas.TDY = temperature_profile[n], r.thermo.DP[0], r.thermo.Y
            # surf.TD = temperature_profile[n], surf.TD[1]
            r.syncState()
            upstream.syncState()
        else:
            gas.TDY = r.thermo.TDY
            upstream.syncState()
        
        sim.reinitialize()
        if n == CAT_ON_INDEX:
            # turn the surface reactions on as we enter the catalyst region
            surf.set_multiplier(1.0)

            # apply kinetic perturbations to the surface reactions
            if surf_kinetics_perturb:
                for index in surf_kinetics_perturb.keys():
                    surf.set_multiplier(surf_kinetics_perturb[index], index)
                    logging.debug(f"Set multiplier of surface reaction {surf.reactions()[index].equation} to {surf_kinetics_perturb[index]}")

        elif n == CAT_OFF_INDEX:
            # turn the surface reactions off as we exit the catalyst region
            surf.set_multiplier(0.0)

        try:
            sim.initial_time = 0.0
            sim.advance(sim.time + 1e4 * residence_time)
            # add timeout handling here if the simulation takes too long
        except (ct.CanteraError, ct._utils.CanteraError, TimeoutException) as e:
            signal.alarm(0)
        #except (ct.CanteraError, ct._utils.CanteraError) as e:
            logging.error(f"Cantera error at reactor {n}: {e}")
            return gas_out, surf_out, gas_rates, surf_rates

        
        kmole_flow_rate = mass_flow_rate / gas.mean_molecular_weight  # kmol/s
        gas_out[n, :] = 1000 * 60 * kmole_flow_rate * gas.X.copy()  # molar flow rate in moles/minute
        surf_out[n, :] = surf.X.copy()
        gas_rates[n, :] = gas.net_rates_of_progress
        surf_rates[n, :] = surf.net_rates_of_progress
        # T_array[n] = surf.T
        
        # Reset the alarm
        signal.alarm(0)

    return gas_out, surf_out, gas_rates, surf_rates


def plot_gas_results(gas_out, gas, outfile=None, T_array=None):
    """Helper function for plotting the results of the gas phase simulation"""

    if isinstance(gas, str):
        gas = ct.Solution(gas)

    # Get indices of key species
    i_CH4 = get_i_thing({'C': 1.0, 'H': 4.0}, gas)
    i_O2 = get_i_thing({'O': 2.0}, gas)
    i_H2 = get_i_thing({'H': 2.0}, gas)
    i_H2O = get_i_thing({'H': 2.0, 'O': 1.0}, gas)
    i_CO = get_i_thing({'C': 1.0, 'O': 1.0}, gas)
    i_CO2 = get_i_thing({'C': 1.0, 'O': 2.0}, gas)

    # Plot the result
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    linewidth2 = 2.0

    fig = plt.figure()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)

    # Simulation Results
    plt.plot(dist_array, gas_out[:, i_O2], label='O2', color=colors[0], linewidth=linewidth2)
    plt.plot(dist_array, gas_out[:, i_CH4], label='CH4', color=colors[1], linewidth=linewidth2)
    plt.plot(dist_array, gas_out[:, i_CO], label='CO', color=colors[4], linewidth=linewidth2)
    plt.plot(dist_array, gas_out[:, i_CO2], label='CO2', color=colors[5], linewidth=linewidth2)
    plt.plot(dist_array, gas_out[:, i_H2], label='H2', color=colors[2], linewidth=linewidth2)
    plt.plot(dist_array, gas_out[:, i_H2O], label='H2O', color=colors[3], linewidth=linewidth2)

    # Experimental data
    plt.plot((df['Distance (mm)'].values - 10.0) / 1000.0, df['O2 (mol/min)'].values, linestyle='dashed', label='EXP O2', color=colors[0])
    plt.plot((df['Distance (mm)'].values - 10.0) / 1000.0, df['CH4 (mol/min)'].values, linestyle='dashed', label='EXP CH4', color=colors[1])
    plt.plot((df['Distance (mm)'].values - 10.0) / 1000.0, df['H2 (mol/min)'].values, linestyle='dashed', label='EXP H2', color=colors[2])
    plt.plot((df['Distance (mm)'].values - 10.0) / 1000.0, df['H2O (mol/min)'].values, linestyle='dashed', label='EXP H2O', color=colors[3])
    plt.plot((df['Distance (mm)'].values - 10.0) / 1000.0, df['CO (mol/min)'].values, linestyle='dashed', label='EXP CO', color=colors[4])
    plt.plot((df['Distance (mm)'].values - 10.0) / 1000.0, df['CO2 (mol/min)'].values, linestyle='dashed', label='EXP CO2', color=colors[5])
    
    ax1 = plt.gca()
    ylim = ax1.get_ylim()
    # Mark the catalyst extent
    plt.plot([dist_array[CAT_ON_INDEX], dist_array[CAT_ON_INDEX]], [0, 1.0], linestyle='--', color='xkcd:grey')
    plt.plot([dist_array[CAT_OFF_INDEX], dist_array[CAT_OFF_INDEX]], [0, 1.0], linestyle='--', color='xkcd:grey')
    plt.ylim(ylim)
    plt.title('Gas Phase Molar Flow Rates Along Reactor')

    ax2 = ax1.twinx()
    if T_array is not None:
        ax2.plot(dist_array, T_array, label='Temperature', color='k', linestyle=':', linewidth=linewidth2)
        ax2.set_ylabel('Temperature (K)')
        ax2.yaxis.get_major_formatter().set_useOffset(False)
    else:
        ax2.plot(dist_array, temperature_profile, label='Temperature', color='k', linestyle=':', linewidth=linewidth2)
        
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('Flow (mol/min)')
    ax1.legend(bbox_to_anchor=(1.15, 0.5))
    # ax1.set_ylim((-0.005686988947011412, 0.1))
    # ax1.set_xlim((-0.0004950495049504951, 0.010396039603960397))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def plot_top_n(result, phase, n=10, outfile=None, species=True, title=None, xlim=None, ylim=None):
    """Helper function for plotting the top n species/reaction rates in a Cantera phase"""

    max_value = np.max(np.abs(result), axis=0)
    indices = np.arange(result.shape[1])
    sorted_order = [x for _, x in sorted(zip(max_value, indices))][::-1]

    # infer surface or gas from presence of X in names
    surface = 'X' in phase.species_names[0]

    fig = plt.figure()
    fig.set_figheight(FIG_HEIGHT)
    fig.set_figwidth(FIG_WIDTH)
    for i in range(n):
        j = sorted_order[i]
        label = phase.species_names[j] if species else phase.reactions()[j].equation
        plt.plot(dist_array, result[:, j], label=label)

    if title is not None:
        plt.title(title)
    
    plt.xlabel('Distance (m)')
    ylabel = 'Reaction Rate'
    if species:
        ylabel = 'Surface Coverage' if surface else 'Flow (mol/min)'
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1.15, 0.5))

    # draw the catalyst extent
    plt.plot([dist_array[CAT_ON_INDEX], dist_array[CAT_ON_INDEX]], plt.ylim(), linestyle='--', color='xkcd:grey')
    plt.plot([dist_array[CAT_OFF_INDEX], dist_array[CAT_OFF_INDEX]], plt.ylim(), linestyle='--', color='xkcd:grey')

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def plot_top_n_gas_species(result, phase, n=10, outfile=None, xlim=None, ylim=None):
    """Helper function for plotting the top n gas phase species from the PFR run"""
    if isinstance(phase, str):
        phase = ct.Solution(phase)
    plot_top_n(result, phase, n=n, outfile=outfile, species=True, title='Gas Species', xlim=xlim, ylim=ylim)

def plot_top_n_surface_species(result, phase, n=10, outfile=None, xlim=None, ylim=None):
    """Helper function for plotting the top n surface species from the PFR run"""
    if isinstance(phase, str):
        phase = ct.Interface(phase, 'surface1')
    plot_top_n(result, phase, n=n, outfile=outfile, species=True, title='Surface Species', xlim=xlim, ylim=ylim)

def plot_top_n_gas_reactions(result, phase, n=10, outfile=None, xlim=None, ylim=None):
    """Helper function for plotting the top n gas phase reactions from the PFR run"""
    if isinstance(phase, str):
        phase = ct.Solution(phase)
    plot_top_n(result, phase, n=n, outfile=outfile, species=False, title='Gas Reactions', xlim=xlim, ylim=ylim)

def plot_top_n_surface_reactions(result, phase, n=10, outfile=None, xlim=None, ylim=None):
    """Helper function for plotting the top n surface reactions from the PFR run"""
    if isinstance(phase, str):
        phase = ct.Interface(phase, 'surface1')
    plot_top_n(result, phase, n=n, outfile=outfile, species=False, title='Surface Reactions', xlim=xlim, ylim=ylim)
