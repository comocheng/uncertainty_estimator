import PEUQSE
import PEUQSE.UserInput
import numpy as np
import sim_wrapper
import yaml
import os

working_dir = '/home/moon/uncertainty_estimator/bpe/simulation/'

experimental_yaml_file = os.path.join(working_dir, 'experiment.yaml')
prior_yaml_file = os.path.join(working_dir, 'prior.yaml')
with open(experimental_yaml_file) as f:
    experimental_data = yaml.safe_load(f)

with open(prior_yaml_file) as f:
    prior_data = yaml.safe_load(f)


# get global info about the indices for key species
# TODO make sure this is only defined in one place
# mech_yaml = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20251229/cantera/chem_annotated.yaml'
mech_yaml = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20251229/cantera/chem_annotated_noCH4X.yaml'
with open(mech_yaml) as f:
    mech_data = yaml.safe_load(f)
species_compositions = [x['composition'] for x in mech_data['species']]

def get_i_thing(ref_composition, composition_list):
    """Helper function for getting the index of a species in a Cantera phase given its composition"""
    for i in range(len(composition_list)):
        if composition_list[i] == ref_composition:
            return i
    assert False, f"Could not find species with composition {ref_composition} in mechanism {mech_yaml}"

i_CH4 = get_i_thing({'C': 1.0, 'H': 4.0}, species_compositions)
i_O2 = get_i_thing({'O': 2.0}, species_compositions)
i_H2O = get_i_thing({'H': 2.0, 'O': 1.0}, species_compositions)
i_H2 = get_i_thing({'H': 2.0}, species_compositions)
i_CO = get_i_thing({'C': 1.0, 'O': 1.0}, species_compositions)
i_CO2 = get_i_thing({'C': 1.0, 'O': 2.0}, species_compositions)

gas_index_order = [i_CH4, i_O2, i_H2O, i_H2, i_CO, i_CO2]
print("Gas index order:", gas_index_order)
gas_names_order = ['CH4 (mol/min)', 'O2 (mol/min)', 'H2O (mol/min)', 'H2 (mol/min)', 'CO (mol/min)', 'CO2 (mol/min)']

# convert experimental data to the format expected by PEUQSE
experimental_y_values = []
experimental_y_uncertainties = []
for key in gas_names_order:
    experimental_y_values.append(experimental_data[key])
    experimental_y_uncertainties.append(experimental_data['Uncertainty ' + key])
experimental_y_values = np.vstack(experimental_y_values)
experimental_y_uncertainties = np.vstack(experimental_y_uncertainties)

# Set up PEUQSE basic info
PEUQSE.UserInput.responses['responses_abscissa'] = experimental_data['Distance (m)']
PEUQSE.UserInput.responses['responses_observed'] = experimental_y_values
PEUQSE.UserInput.responses['responses_observed_uncertainties'] = experimental_y_uncertainties

# plot options
PEUQSE.UserInput.simulated_response_plot_settings['x_label'] = "x"
PEUQSE.UserInput.simulated_response_plot_settings['y_label'] = ["y1", "y2", "y3", "y4", "y5", "y6"]


PEUQSE.UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = prior_data['species_names']

# model priors
PEUQSE.UserInput.model['InputParameterPriorValues'] = [prior_data[key] for key in prior_data['species_names']]
PEUQSE.UserInput.model['InputParametersPriorValuesUncertainties'] = prior_data['cov_J2_mol2']


PEUQSE.UserInput.model['simulateByInputParametersOnlyFunction'] = sim_wrapper.simulation_wrapper

PEUQSE.UserInput.parameter_estimation_settings['mcmc_length'] = 3

PE_object = PEUQSE.parameter_estimation(PEUQSE.UserInput)
PE_object.doMetropolisHastings()
# PE_object.doEnsembleSliceSampling()

# PE_object.doOptimizeSSR(method="BFGS") #Note: using method=Nelder-Mead doesn't change the final result much.


PE_object.createAllPlots()

# PE_object.save_to_dill('PE_object_00')
