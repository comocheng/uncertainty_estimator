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

# Trying to make Cantera do sensitivity CSVs like the SimulationProfileWriter in the uncertainty tool


# +
import rmgpy.chemkin
import os

import csv

import cantera as ct
import numpy as np
import logging

import rmgpy.solver  # import SimpleReactor, TerminationTime, SurfaceReactor
import rmgpy.quantity  #import Quantity
import rmgpy.rmg.listener # import SimulationProfileWriter, SimulationProfilePlotter
import rmgpy.rmg.settings  #import ModelSettings, SimulatorSettings
import rmgpy.util
import rmgpy.tools.plot

import matplotlib
import matplotlib.pyplot as plt
# %matplotlib inline

logging.basicConfig(level=logging.DEBUG)


# -

sevy_mech = '/home/moon/rmg/RMG-Py/examples/rmg/superminimal/chemkin'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict)






# # Cantera version

def same_reaction(rmg_rxn, ct_rxn):
    rmg_r = set([str(x.to_chemkin()) for x in rmg_rxn.reactants])
    rmg_p = set([str(x.to_chemkin()) for x in rmg_rxn.products])

    ct_r = set(ct_rxn.reactants.keys())
    ct_p = set(ct_rxn.products.keys())
    return rmg_r == ct_r and rmg_p == ct_p


# +
# gas = ct.Solution('/home/moon/rmg/RMG-Py/examples/rmg/superminimal/cantera/chem_annotated.yaml')
gas = ct.Solution(
    thermo='IdealGas',
    kinetics='GasKinetics',
    species=[x.to_cantera(use_chemkin_identifier=True) for x in species_list],
    reactions=[x.to_cantera(use_chemkin_identifier=True) for x in reaction_list],
)

if len(reaction_list) != gas.n_reactions:
    raise NotImplementedError

for i in range(len(species_list)):
#     print(species_list[i].to_chemkin(), gas.species_names[i])
    assert str(species_list[i].to_chemkin()) == str(gas.species_names[i])

for i in range(len(reaction_list)):
#     print(reaction_list[i], gas.reaction_equations()[i])
    assert same_reaction(reaction_list[i], gas.reactions()[i])


# +
output_directory = 'cantera_csvs'
os.makedirs(output_directory, exist_ok=True)

os.makedirs(os.path.join(output_directory, 'solver'), exist_ok=True)
# -



# +
# gas = ct.Solution('/home/moon/rmg/RMG-Py/examples/rmg/superminimal/cantera/chem_annotated.yaml')
gas = ct.Solution(
    thermo='IdealGas',
    kinetics='GasKinetics',
    species=[x.to_cantera(use_chemkin_identifier=True) for x in species_list],
    reactions=[x.to_cantera(use_chemkin_identifier=True) for x in reaction_list],
)



T = 1500.0
P = 100000.0
initial_mole_fractions='H2(1):0.67, O2(2):0.33'
gas.TPX = T, P, initial_mole_fractions

reaction_system_index = 0
sens_worksheet = []
sensitive_species = [species_list[4]]  # RMG species H2

for spec in sensitive_species:
    csvfile_path = os.path.join(output_directory, 'solver', f'sensitivity_1_SPC_{spec.index}.csv')
    sens_worksheet.append(csvfile_path)

    
    
reactor = ct.IdealGasConstPressureReactor(gas, energy='off')  # isothermal
net = ct.ReactorNet([reactor])


# set the tolerances for the solution and for the sensitivity coefficients
net.rtol = 1.0e-6
net.atol = 1.0e-15
net.rtol_sensitivity = 1.0e-6
net.atol_sensitivity = 1.0e-6

sensitivity_threshold = net.atol_sensitivity

# find out what is the most sensitive parameter to the H2 concentration
for i in range(gas.n_reactions):
    reactor.add_sensitivity_reaction(i)

for i in range(gas.n_species):
    reactor.add_sensitivity_species_enthalpy(i)

ts_ct = [net.time]
ys_ct = [reactor.thermo.X]
senss = [np.zeros(gas.n_reactions + gas.n_species)]

while net.time < 1e-1:
    net.step()
    ts_ct.append(net.time)
    time_array = np.array(ts_ct)
    ys_ct.append(reactor.thermo.X)
    sens_mat = np.zeros(gas.n_reactions + gas.n_species)
    for i in range(gas.n_reactions):
        sens_mat[i] = net.sensitivity('H2(1)', i)
    for i in range(gas.n_species):
        sens_mat[gas.n_reactions + i] = net.sensitivity('H2(1)', gas.n_reactions + i)
    senss.append(sens_mat)
    
    
    # write sensitivities -- TODO add multiple sensitive species
    
    with open(sens_worksheet[0], 'w') as outfile:
        species_name = rmgpy.chemkin.get_species_identifier(sensitive_species[0])
        headers = ['Time (s)']
        
        worksheet = csv.writer(outfile)
        reactions_above_threshold = []
        for j in range(gas.n_reactions + gas.n_species):
            for k in range(len(ys_ct)):
                if abs(net.sensitivity('H2(1)', j)) > sensitivity_threshold:
                    reactions_above_threshold.append(j)
                    break
        
        # need conversion from Cantera to RMG and back
        headers.extend([f'dln[{species_name}]/dln[k{j + 1}]: {reaction_list[j].to_chemkin(kinetics=False)}' if j < gas.n_reactions
                        else f'dln[{species_name}]/dG[{rmgpy.chemkin.get_species_identifier(gas.species()[j - gas.n_reactions])}]' for j in reactions_above_threshold])
        worksheet.writerow(headers)

        for k in range(len(time_array)):
            row = [time_array[k]]
            row.extend([senss[k][j] for j in reactions_above_threshold])
#             row.extend([norm_sens_array[i][k][j] for j in reactions_above_threshold])
            worksheet.writerow(row)
    

# plot results
for i in range(len(ys_ct[0])):
    y_hs = [x[i] for x in ys_ct]
    plt.plot(ts_ct, y_hs, label=gas.species_names[i])
plt.xlabel('time (s)')
plt.ylabel('concentration')
plt.legend()
# -



# plot reaction and species sensitivity results
for i in range(len(senss[0])):
    if i < gas.n_reactions:
        label = gas.reaction_equations()[i]
        style = 'solid'
    else:
        label = gas.species_names[i - gas.n_reactions]
        style = 'dashed'
    y_hs = [x[i] for x in senss]
    plt.plot(ts_ct, y_hs, label=label, linestyle=style)
plt.xlabel('time (s)')
plt.ylabel('sensitivity')
plt.legend(bbox_to_anchor=(1.1, 1.05))

# plot sensitivity results
for i in range(len(senss[0])):
    y_hs = [x[i] for x in senss]
    plt.plot(ts_ct, y_hs, label=gas.reaction_equations()[i])
plt.xlabel('time (s)')
plt.ylabel('sensitivity')
plt.legend(bbox_to_anchor=(1.1, 1.05))

gas.reaction_equations

senss



# +
# can we repeat with species?

gas = ct.Solution('/home/moon/rmg/RMG-Py/examples/rmg/superminimal/cantera/chem_annotated.yaml')

T = 1500.0
P = 100000.0
initial_mole_fractions='H2(1):0.67, O2(2):0.33'
gas.TPX = T, P, initial_mole_fractions





reactor = ct.IdealGasConstPressureReactor(gas, energy='off')  # isothermal
net = ct.ReactorNet([reactor])


# set the tolerances for the solution and for the sensitivity coefficients
net.rtol = 1.0e-6
net.atol = 1.0e-15
net.rtol_sensitivity = 1.0e-6
net.atol_sensitivity = 1.0e-6


# find out what is the most sensitive parameter to the H2 concentration
for i in range(gas.n_reactions):
    reactor.add_sensitivity_reaction(i)


ts_ct = [net.time]
ys_ct = [reactor.thermo.X]
senss = [np.zeros(gas.n_reactions)]

while net.time < 1e-1:
    net.step()
    ts_ct.append(net.time)
    ys_ct.append(reactor.thermo.X)
    sens_mat = np.zeros(gas.n_reactions)
    for i in range(gas.n_reactions):
        sens_mat[i] = net.sensitivity('H2(1)', i)
    senss.append(sens_mat)

# plot results
for i in range(len(ys_ct[0])):
    y_hs = [x[i] for x in ys_ct]
    plt.plot(ts_ct, y_hs, label=gas.species_names[i])
plt.xlabel('time (s)')
plt.ylabel('concentration')
plt.legend()
# -

net.sensitivity_parameter_name(4)

net.get_state()





len(ts_ct)

len(reaction_list)



# # RMG simple.pyx

# +
T = (1500, 'K')
P = (1, 'bar')
termination_time = (1.0e-1, 's')
termination = [rmgpy.solver.TerminationTime(rmgpy.quantity.Quantity(termination_time))]
sensitivity_threshold=1e-3

# display(species_list[4], species_list[5])
initial_mole_fractions={
    species_list[4]: 0.67,  # H2
    species_list[5]: 0.33,  # O2
}
sensitive_species = [species_list[4]]
# sensitive_species = []


reaction_system = rmgpy.solver.SimpleReactor(
    T=T,
    P=P,
    initial_mole_fractions=initial_mole_fractions,
    termination=termination,
    sensitive_species=sensitive_species,
    sensitivity_threshold=sensitivity_threshold
)


reaction_system.initialize_model(
    species_list,
    reaction_list,
    [],
    [],
)

# +
# Create the csv worksheets for logging sensitivity
output_directory = 'sensitivity_calcs'
os.makedirs(output_directory, exist_ok=True)
rmgpy.util.make_output_subdirectory(output_directory, 'solver')
sens_worksheet = []
reaction_system_index = 0
for spec in reaction_system.sensitive_species:
    csvfile_path = os.path.join(output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1, spec.index))
    sens_worksheet.append(csvfile_path)

reaction_system.attach(rmgpy.rmg.listener.SimulationProfileWriter(
    output_directory, reaction_system_index, species_list))
reaction_system.attach(rmgpy.rmg.listener.SimulationProfilePlotter(
    output_directory, reaction_system_index, species_list))

simulator_settings = rmgpy.rmg.settings.SimulatorSettings()  # defaults

model_settings = rmgpy.rmg.settings.ModelSettings()  # defaults
model_settings.tol_move_to_core = 0.1
model_settings.tol_interrupt_simulation = 1.0
model_settings.tol_keep_in_edge = 0.0

print('starting simulation')
reaction_system.simulate(
    core_species=species_list,
    core_reactions=reaction_list,
    edge_species=[],
    edge_reactions=[],
    surface_species=[x for x in species_list if x.contains_surface_site()],
    surface_reactions=[x for x in reaction_list if x.is_surface_reaction()],
    model_settings=model_settings,
    simulator_settings=simulator_settings,
    sensitivity=True,
    sens_worksheet=sens_worksheet,
)
# -

reaction_system_index = 0
for sens_species in sensitive_species:
    csvfile_path = os.path.join(output_directory, 'solver',
                                'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index+1,
                                                                     sens_species.index))
    time, data_list = rmgpy.tools.plot.parse_csv_data(csvfile_path)

times = time.data

for data in data_list:
    if 'dG' in data.label:
        continue
#     print(data.label)

    plt.plot(times, data.data, label=data.label.split()[-1])
plt.xlabel('time (s)')
plt.ylabel('sensitivity')
plt.legend(bbox_to_anchor=(1.1, 1.05))







# +

# Options controlling the ODE simulations:
initial_time = 1e-12  # The time at which to initiate the simulation, in seconds
time_step = 10 ** 0.1  # The multiplicative factor to use between consecutive time points
time_step = 10 ** 0.05  # The multiplicative factor to use between consecutive time points
abs_tol = 1e-16  # The absolute tolerance to use in the ODE simluations
rel_tol = 1e-8  # The relative tolerance to use in the ODE simulations

terminated = False
next_time = initial_time
step_time = initial_time

times = [reaction_system.t]
# ys = [reaction_system.y]
ys = [reaction_system.core_species_concentrations]
ddts = [reaction_system.compute_rate_derivative()]


while reaction_system.t < termination_time[0]:
    # Integrate forward in time to the next time point
    reaction_system.step(step_time)

    if reaction_system.t >= 0.9999 * next_time:
        next_time *= time_step
        times.append(reaction_system.t)
        ys.append(reaction_system.core_species_concentrations)
        
        norm_sens[j] = 1/volume * (mole_sens[j*num_core_species+sens_species_indices[i]]-c*dVdk[j]) * forward_rate_coefficients[j]/c
        ddts.append(reaction_system.compute_rate_derivative())
#         ys.append(reaction_system.y)

    # Increment destination step time if necessary
    if reaction_system.t >= 0.9999 * step_time:
        step_time *= 10.0
# -

len(ddts)

ddts[0].shape

len(times)

species_list[4]

reaction_system.V

for i in range(len(reaction_list)):
    y_hs = [x[4, i] / reaction_system.V for x in ddts]  #it's reactions first, then spec
    plt.plot(times, y_hs, label=str(reaction_list[i]))
plt.xlabel('time (s)')
plt.ylabel('sensitivity')
plt.legend(bbox_to_anchor=(1.1, 1.05))








        # Create the csv worksheets for logging sensitivity
        util.make_output_subdirectory(self.output_directory, 'solver')
        sens_worksheet = []
        reaction_system_index = 0
        for spec in reaction_system.sensitive_species:
            csvfile_path = os.path.join(self.output_directory, 'solver',
                                        'sensitivity_{0}_SPC_{1}.csv'.format(reaction_system_index + 1, spec.index))
            sens_worksheet.append(csvfile_path)

        reaction_system.attach(SimulationProfileWriter(
            self.output_directory, reaction_system_index, self.species_list))
        reaction_system.attach(SimulationProfilePlotter(
            self.output_directory, reaction_system_index, self.species_list))

        simulator_settings = SimulatorSettings()  # defaults

        model_settings = ModelSettings()  # defaults
        model_settings.tol_move_to_core = 0.1
        model_settings.tol_interrupt_simulation = 1.0
        model_settings.tol_keep_in_edge = 0.0

        print('starting simulation')
        reaction_system.simulate(
            core_species=self.species_list,
            core_reactions=self.reaction_list,
            edge_species=[],
            edge_reactions=[],
            surface_species=[x for x in self.species_list if x.contains_surface_site()],
            surface_reactions=[x for x in self.reaction_list if x.is_surface_reaction()],
            model_settings=model_settings,
            simulator_settings=simulator_settings,
            sensitivity=True,
            sens_worksheet=sens_worksheet,
        )

reaction_system.compute_rate_derivative().shape

len(reaction_list)

len(species_list)

species_list[4]

for i in range(len(ys[0])):
    y_hs = [x[i] for x in ys]
    plt.plot(times, y_hs / np.sum(ys[0]), label=str(species_list[i]))
plt.xlabel('time (s)')
plt.ylabel('concentration')
plt.legend()

# +
# plot on top of each other

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

assert len(ys[0]) == len(ys_ct[0])

for i in range(len(ys[0])):
    y_hs = [x[i] for x in ys]
    plt.plot(times, y_hs / np.sum(ys[0]), marker='.', color=colors[i % len(colors)], label=str(species_list[i]), linestyle='dotted')
    
    y_hs = [x[i] for x in ys_ct]
    plt.plot(ts_ct, y_hs, color=colors[i % len(colors)], linestyle='dashed', alpha=0.5)
    
    
    
plt.xlabel('time (s)')
plt.ylabel('concentration')
plt.legend()
# -




