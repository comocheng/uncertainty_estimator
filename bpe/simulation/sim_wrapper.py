import yaml
import simulation
import numpy as np


with open('experiment.yaml') as f:
    experimental_data = yaml.safe_load(f)
with open('prior.yaml') as f:
    prior_data = yaml.safe_load(f)

my_prior_species_indices = prior_data['species_indices']
my_prior_reaction_indices = prior_data['reaction_indices']
output_gas_species_indices = experimental_data['out_gas_indices']

# Make hardcoded list of distance indices to copy into your simulation_wrapper (maybe this can be a global later...)
interpolation_indices = []
for i in range(len(experimental_data['Distance (m)'])):
    interpolation_indices.append(int(np.argmin(np.abs(simulation.dist_array - experimental_data['Distance (m)'][i]))))

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
    for i in output_gas_species_indices:
        species_results = []
        for j in interpolation_indices:
            species_results.append(gas_out[j, i])
        results.append(species_results)
        
    results = np.vstack(results)

    return results
