import yaml
import pandas as pd
import simulation
import numpy as np


observed_data_y = pd.read_csv('experiment.csv').values
prior_data = pd.read_csv('priors.csv').values

with open('sim_info.yaml') as f:
    sim_info = yaml.safe_load(f)

my_prior_species_indices = sim_info['prior_species_indices']
my_prior_reaction_indices = sim_info['prior_reaction_indices']
output_gas_species_indices = sim_info['out_gas_indices']
sample_distances = sim_info['sample_distances']

# map the sample distances to indices in the reactor
interpolation_indices = []
for i in range(len(sample_distances)):
    interpolation_indices.append(int(np.argmin(np.abs(simulation.dist_array - sample_distances[i]))))


# wrapper for PEUQSE to call monolith simulation and return results in the format expected by PEUQSE
def simulation_wrapper(params):
    # params are a list of perturbations to thermo (in units of J/mol) and kinetics (multiplier)
    # see prior.yaml for more info on priors

    mech_yaml = 'chem_annotated.yaml'

    thermo_perturb = {}
    for i, param in enumerate(params[:len(my_prior_species_indices)]):
        thermo_perturb[my_prior_species_indices[i]] = param  # expecting units of J/mol

    kinetics_perturb = {}
    for i, param in enumerate(params[len(my_prior_species_indices):]):
        kinetics_perturb[my_prior_reaction_indices[i]] = np.float_power(10.0, param)  # expecting units of multiplier

    gas_out, surf_out, gas_rates, surf_rates = simulation.run_simulation(
        mech_yaml,
        surf_thermo_perturb=thermo_perturb,
        surf_kinetics_perturb=kinetics_perturb,
    )

    # extract the results for the species of interest and the distances
    results = []
    for i in interpolation_indices:
        distance_results = []  # results at a given distance in the reactor
        for j in output_gas_species_indices:
            distance_results.append(gas_out[i, j])
        results.append(distance_results)
    
    # how do I make sure this is handled in the same order as observed_data_y?
    # start by making a picture of what it's supposed to be:
    # O2, CH4 H2, CO, CO2, H2O
    # y1, y1, y1, y1, y1, y1
    # y2, y2, y2, y2, y2, y2
    # ...

    # This really needs a good test

    results = np.vstack(results)
    results = results.ravel()

    print(results.shape, results)

    return results
