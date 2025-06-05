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

# +
# trying to figure out which line is causing the holdup in my uncertainty example
# something is trying to load Julia and there is no need for it until you get to the actual reactor
import os
import sys
import pickle
import copy
import numpy as np
import rmgpy.chemkin
import rmgpy.tools.uncertainty
import rmgpy.kinetics.uncertainties

import rmgpy.tools.canteramodel
import random

import rmgpy.kinetics
import matplotlib.pyplot as plt
# %matplotlib inline


import rmgpy.data.rmg
import rmgpy


import psutil

# +
# Load the model

# Must use annotated chemkin file
chemkin_file = '/home/moon/uncertainty_estimator/kmc_presentation_20241115/ethane/chemkin/chem_annotated.inp'
dict_file = '/home/moon/uncertainty_estimator/kmc_presentation_20241115/ethane/chemkin/species_dictionary.txt'

# Run Gao estimation of input parameters (takes a long time to load database)
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='uncertainty_calculations')
uncertainty.load_model(chemkin_file, dict_file)
# -

uncertainty.database = rmgpy.data.rmg.RMGDatabase()


# +
kinetics_depositories = None
thermo_libraries = ['primaryThermoLibrary', 'BurkeH2O2']
reaction_libraries = ['BurkeH2O2inN2']
kinetics_families = 'default'


if not kinetics_depositories:
    kinetics_depositories = ['training']
if not thermo_libraries:
    thermo_libraries = ['primaryThermoLibrary']
if not reaction_libraries:
    reaction_libraries = []

uncertainty.database.load(
    rmgpy.settings['database.directory'],
    kinetics_families=kinetics_families,
    kinetics_depositories=kinetics_depositories,
    thermo_libraries=thermo_libraries,
    reaction_libraries=reaction_libraries,
)
# -

for familyLabel, family in uncertainty.database.kinetics.families.items():
    if not family.auto_generated:
        if familyLabel in ['1,4_Cyclic_birad_scission', 'Birad_R_Recombination', 'Intra_ene_reaction',
                          'R_Addition_COm', 'R_Addition_MultipleBond']:
            print(f'skipping {familyLabel}')
            continue
        print(familyLabel)
        family.add_rules_from_training(thermo_database=uncertainty.database.thermo)
        family.fill_rules_by_averaging_up(verbose=True)

psutil.cpu_count()



# +
maxproc = psutil.cpu_count()

def determine_procnum_from_ram():
    """
    Get available RAM (GB)and procnum dependent on OS.
    """
    if sys.platform.startswith("linux"):
        # linux
        memory_available = psutil.virtual_memory().free / (1000.0**3)
        memory_use = psutil.Process(os.getpid()).memory_info()[0] / (1000.0**3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
    elif sys.platform == "darwin":
        # OS X
        memory_available = psutil.virtual_memory().available / (1000.0**3)
        memory_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1000.0**3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
    else:
        # Everything else
        procnum = 1

    # Return the maximal number of processes for multiprocessing
    return procnum


# -

determine_procnum_from_ram()

__name__


