"""Water gas shift simulation, CO + H2O <=> CO2 + H2"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import copy
import logging
import os


# SIM PROPERTIES
cat_diameter = 0.013
area = np.pi * np.float_power(cat_diameter / 2.0, 2.0) * 2.0 + 0.001 * np.pi * cat_diameter  # not actually clear whether this is 2-sided or not


volume = 0.0291666667  # m^3 volume of batch reactor
cat_area = 2.65 / 100.0 / 100.0  # convert from cm^2 to m^2

P_CO_torr = 25
P_H2O_torr = 15

P_CO = P_CO_torr * ct.one_atm / 760
P_H2O = P_H2O_torr * ct.one_atm / 760
P_total = P_CO + P_H2O

x_CO = P_CO / P_total
x_H2O = P_H2O / P_total


def get_i_thing(ref_composition, phase):
    """Helper function for getting the index of a species in a Cantera phase given its composition"""
    for i in range(phase.n_species):
        if phase.species()[i].composition == ref_composition:
            return i
    assert False, f"Could not find species with composition {ref_composition} in phase {phase.name}"


def increase_enthalpy(phase, species_index, increase_enthalpy_J_per_mol):
    """Helper function for increasing the enthalpy of a species in a Cantera phase by a specified amount (in J/mol)"""
    if ct.__version__ in ['3.2.0']:
        data_copy = phase.species()[species_index].input_data.copy()
    else:    
        data_copy = copy.deepcopy(phase.species()[species_index].input_data)
    for i in range(len(data_copy['thermo']['data'])):
        data_copy['thermo']['data'][i][5] += increase_enthalpy_J_per_mol * 1000.0 / ct.gas_constant  # Cantera needs J / kmol
    new_sp = ct.Species().from_dict(data_copy)
    phase.modify_species(species_index, new_sp)




def run_simulation(
    mech_yaml=None,
    gas=None,
    surf=None,
    T=575,
    surf_thermo_perturb=None,
    surf_kinetics_perturb=None,
):
    """Run a simulation of the Flaherty WGS experiment

    Parameters
    ----------
    mech_yaml : str
        Path to the Cantera mechanism file
    gas : cantera.Solution
        Cantera gas object
    surf : cantera.Interface
        Cantera surface object
    surf_thermo_perturb : dict
        key is index of the species in the surface phase, value is the amount to perturb the enthalpy in J/mol
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
    i_CO = get_i_thing({'C': 1.0, 'O': 1.0}, gas)
    i_H2O = get_i_thing({'H': 2.0, 'O': 1.0}, gas)
    i_CO2 = get_i_thing({'C': 1.0, 'O': 2.0}, gas)
    i_H2 = get_i_thing({'H': 2.0}, gas)
    i_X = get_i_thing({'X': 1.0}, surf)


    # Initialize the reactor
    X = f'{gas.species_names[i_H2O]}: {x_H2O}, {gas.species_names[i_CO]}: {x_CO}'

    gas.TPX = T, P_total, X
    surf.TP = T, P_total

    initial_coverages = np.zeros_like(surf.coverages)
    initial_coverages[i_X] = 1.0
    surf.coverages = initial_coverages

    r = ct.IdealGasReactor(gas, energy='off')
    r.volume = volume
    rsurf = ct.ReactorSurface(surf, r, A=cat_area)
    sim = ct.ReactorNet([r])

    t_end = 22 * 60  # 10 minutes, we won't worry about catalyst poisoning? or that's not what's happening

    times = [0]
    concs = [gas.X]
    Ps = [gas.P]
    Ts = [gas.T]
    total_moles_H2 = [gas.density_mole * 1000.0 * gas.X[i_H2] * r.volume]
    molar_densities = [gas.density_mole * 1000.0]
    while sim.time < t_end:
        try:
            sim.step()
        except ct.CanteraError:
            return times, concs, Ps, Ts, molar_densities
        times.append(sim.time)
        concs.append(gas.X)
        Ps.append(gas.P)
        Ts.append(gas.T)
        total_moles_H2.append(gas.density_mole * 1000.0 * gas.X[i_H2] * r.volume)
        molar_densities.append(gas.density_mole * 1000.0)
    times = np.array(times)
    concs = np.array(concs)
    Ps = np.array(Ps)
    total_moles_H2 = np.array(total_moles_H2)
    molar_densities = np.array(molar_densities)
    Ts = np.array(Ts)

    return times, concs, Ps, Ts, molar_densities
    # TODO add error handling

def get_partial_pressures_mTorr(concs, Ps):
    partial_pressures = np.zeros_like(concs)
    for i in range(concs.shape[1]):
        partial_pressures[:, i] = np.multiply(concs[:, i], Ps) * 760000 / ct.one_atm
    return partial_pressures