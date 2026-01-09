import simulation
import numpy as np

# TODO move this hard coding to its proper place in the setup script

# hard code the indices of the parameters? Don't want to be looking this up inside the simulation wrapper
my_prior_species_indices = [5, 3, 1, 7, 9]  # shouldn't be hard coded here
my_prior_reaction_indices = [0, 1, 3]
my_output_gas_species_indices = [3, 4, 6, 7, 8, 5]

# Make hardcoded list of distance indices to copy into your simulation_wrapper (maybe this can be a global later...)
# interpolation_indices = []
# for i in range(len(dist_array)):
#     interpolation_indices.append(int(np.argmin(np.abs(simulation.dist_array - dist_array[i]))))
# print(interpolation_indices)
interpolation_indices = [211, 238, 266, 294, 321, 349, 377, 404, 432, 460, 488, 515, 543, 571, 598, 626, 654, 681, 709, 737]


# wrapper for PEUQSE to call monolith simulation and return results in the format expected by PEUQSE
def simulation_wrapper(params):
    # params are a list of perturbations to thermo (in units of J/mol) and kinetics (multiplier)
    # see prior.yaml for more info on priors

    mech_yaml = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20251229/cantera/chem_annotated.yaml'
    mech_yaml = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20251229/cantera/chem_annotated_noCH4X.yaml'

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
    for i in my_output_gas_species_indices:
        species_results = []
        for j in interpolation_indices:
            species_results.append(gas_out[j, i])
        results.append(species_results)
        
    results = np.vstack(results)

    return results