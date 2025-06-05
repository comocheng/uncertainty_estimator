import os
import sys
import time
import cantera as ct
import numpy as np
import pandas as pd
import concurrent.futures
import rmgpy.chemkin
import subprocess




working_dir = '/home/moon/uncertainty_estimator/nam29_abstract/cpox_pt/perturbed_mechs'
results_dir = os.path.join(working_dir, 'results')
os.makedirs(results_dir, exist_ok=True)


base_yaml_path = os.path.join(working_dir, 'base.yaml')
perturbed_yaml_path = os.path.join(working_dir, 'perturbed.yaml')

def run_simulation(gas, surf):
    """
    Verbose prints out values as you go along
    Sens is for sensitivity, in the form [perturbation, reaction #]
    """

    ratio = 1.0
    fo2 = 1 / (2. * ratio + 1 + 79 / 21)
    fch4 = 2 * fo2 * ratio
    far = 79 * fo2 / 21
    ratio_in = [fch4, fo2, far]

    
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

    temp = t_in
    mol_in = ratio_in

    verbose = False
    sens = False
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

        
    gas_out = np.array(gas_out)
    surf_out = np.array(surf_out)
    gas_names = np.array(gas_names)
    surf_names = np.array(surf_names)
    data_out = gas_out, surf_out, gas_names, surf_names, dist_array, T_array
    return data_out



for i in range(200):

    print(f'running {i}')
    base_yaml_path = os.path.join(working_dir, f'chem_{i:04}.yaml')
    spec_delay_file_gas = os.path.join(results_dir, f'chem_{i:04}_gas_out.npy')
    spec_delay_file_surf = os.path.join(results_dir, f'chem_{i:04}_surf_out.npy')

    if os.path.exists(spec_delay_file_gas) and os.path.exists(spec_delay_file_surf):
        print(f'{i} already exists')
        continue


    # load the base gas
    base_gas = ct.Solution(base_yaml_path, 'gas')
    base_surf = ct.Interface(base_yaml_path, 'surface1', [base_gas])

    try:
        a = run_simulation(base_gas, base_surf)
    except ct._cantera.CanteraError:
        print(f'Error in {i}')
        continue
    gas_out, surf_out, gas_names, surf_names, dist_array, T_array = a

    np.save(spec_delay_file_gas, gas_out)
    np.save(spec_delay_file_surf, surf_out)
    