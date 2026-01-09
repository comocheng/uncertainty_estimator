import PEUQSE
import PEUQSE.UserInput
import numpy as np
import sys
sys.path.append('/projects/westgroup/harris.se/uncertainty_estimator/bpe/simulation/')
import sim_wrapper
import yaml
import os

working_dir = os.path.dirname(os.path.abspath(__file__))

# convert experimental data to the format expected by PEUQSE
experimental_yaml_file = os.path.join(working_dir, 'experiment.yaml')
with open(experimental_yaml_file) as f:
    experimental_data = yaml.safe_load(f)
gas_index_order = experimental_data['out_gas_indices']
gas_names_order = experimental_data['out_gas_names']
experimental_y_values = []
experimental_y_uncertainties = []
for key in gas_names_order:
    experimental_y_values.append(experimental_data[key])
    experimental_y_uncertainties.append(experimental_data['Uncertainty ' + key])
experimental_y_values = np.vstack(experimental_y_values)
experimental_y_uncertainties = np.vstack(experimental_y_uncertainties)

prior_yaml_file = os.path.join(working_dir, 'prior.yaml')
with open(prior_yaml_file) as f:
    prior_data = yaml.safe_load(f)
prior_labels = prior_data['species_names'] + prior_data['reaction_equations']


# Set up PEUQSE basic info
PEUQSE.UserInput.responses['responses_abscissa'] = experimental_data['Distance (m)']
PEUQSE.UserInput.responses['responses_observed'] = experimental_y_values
PEUQSE.UserInput.responses['responses_observed_uncertainties'] = experimental_y_uncertainties

# plot options
PEUQSE.UserInput.simulated_response_plot_settings['x_label'] = "Distance (m)"
PEUQSE.UserInput.simulated_response_plot_settings['y_label'] = gas_names_order


PEUQSE.UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = prior_labels

# model priors
PEUQSE.UserInput.model['InputParameterPriorValues'] = [prior_data[key] for key in prior_data['species_names'] + prior_data['reaction_equations']]
PEUQSE.UserInput.model['InputParametersPriorValuesUncertainties'] = prior_data['cov_J2_mol2']


PEUQSE.UserInput.model['simulateByInputParametersOnlyFunction'] = sim_wrapper.simulation_wrapper

PEUQSE.UserInput.parameter_estimation_settings['mcmc_length'] = 3

PE_object = PEUQSE.parameter_estimation(PEUQSE.UserInput)
PE_object.doMetropolisHastings()
# PE_object.doEnsembleSliceSampling()

# PE_object.doOptimizeSSR(method="BFGS") #Note: using method=Nelder-Mead doesn't change the final result much.


PE_object.createAllPlots()

# PE_object.save_to_dill('PE_object_00')
