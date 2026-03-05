import yaml
import pandas as pd
import wgs_simulation
import numpy as np


observed_data_y = pd.read_csv('experiment.csv').values
prior_data = pd.read_csv('priors.csv').values

with open('sim_info.yaml') as f:
    sim_info = yaml.safe_load(f)

my_prior_species_indices = sim_info['prior_species_indices']
my_prior_reaction_indices = sim_info['prior_reaction_indices']
sample_times = sim_info['sample_times']  # seconds


def simulation_wrapper(params):
    # run simulations at 575 and 625 K
    mech_yaml = 'chem_annotated.yaml'
    i_H2 = sim_info['H2_index']
    thermo_perturb = {}
    for i, param in enumerate(params[:len(my_prior_species_indices)]):
        thermo_perturb[my_prior_species_indices[i]] = param  # expecting units of J/mol

    kinetics_perturb = {}
    for i, param in enumerate(params[len(my_prior_species_indices):]):
        kinetics_perturb[my_prior_reaction_indices[i]] = np.float_power(10.0, param)  # expecting units of multiplier

    times, concs, Ps, Ts, molar_densities = wgs_simulation.run_simulation(
        mech_yaml,
        surf_thermo_perturb=thermo_perturb,
        surf_kinetics_perturb=kinetics_perturb,
        T=575
    )
    if times[-1] < sample_times[-1]:
        return np.zeros(10) + np.nan
    partial_pressure_H2_mTorr = np.multiply(concs[:, i_H2], Ps) * 760000 / 101325.0
    results_575 = []
    for i in range(len(sample_times)):
        interpolation_index = int(np.argmin(np.abs(times - sample_times[i])))
        results_575.append(partial_pressure_H2_mTorr[interpolation_index])

    
    times, concs, Ps, Ts, molar_densities = wgs_simulation.run_simulation(
        mech_yaml,
        surf_thermo_perturb=thermo_perturb,
        surf_kinetics_perturb=kinetics_perturb,
        T=625
    )
    if times[-1] < sample_times[-1]:
        return np.zeros(10) + np.nan
    partial_pressure_H2_mTorr = np.multiply(concs[:, i_H2], Ps) * 760000 / 101325.0
    results_625 = []
    for i in range(len(sample_times)):
        interpolation_index = int(np.argmin(np.abs(times - sample_times[i])))
        results_625.append(partial_pressure_H2_mTorr[interpolation_index])


    # how do I make sure this is handled in the same order as observed_data_y?
    # start by making a picture of what it's supposed to be:
    # observed_data_y goes 
    #   625_1, 575_1, 625_2, 575_2, 625_3, 575_3, 625_4, 575_4, 625_5, 575_5

    results = np.vstack([results_625, results_575]).T
    results = results.ravel()

    # print(results.shape, results)

    return results
