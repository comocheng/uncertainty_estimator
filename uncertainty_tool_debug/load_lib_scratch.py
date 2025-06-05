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
# import os
# import sys
# import pickle
# import copy
# import numpy as np
import rmgpy.chemkin
import rmgpy
# import rmgpy.tools.uncertainty
# import rmgpy.kinetics.uncertainties
# import rmgpy.exceptions
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
logging.debug("test")
# import random

import rmgpy.kinetics
# import matplotlib.pyplot as plt
# # %matplotlib inline

# sys.path.append('/home/moon/autoscience/reaction_calculator/database')
# import database_fun
# -



# +
database = rmgpy.data.rmg.RMGDatabase()

database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = [
        'primaryThermoLibrary',
        'FFCM1(-)',
        
    ],
    reaction_libraries = [
#         'FFCM1(-)',
        'GRI-Mech3.0',
#         'JetSurF2.0'
    ],
    kinetics_families = 'default',
#     kinetics_depositories = ['training'],
#     depository = False, # Don't bother loading the depository information, as we don't use it
)

# +

libraryName = 'FFCM1(-)'
print(f'{libraryName}.py')
# -



with open(input_file, encoding='windows-1252') as f:

open(input_file,'rb') as f:
book = f.read().decode(errors='replace')



my_lib = rmgpy.data.kinetics.library.KineticsLibrary()

database.kinetics.libraries['mine'] = my_lib

database.kinetics.libraries['mine'].load('/home/moon/rmg/RMG-database/input/kinetics/libraries/JetSurF2.0/reactions.py')





rmgpy.settings['database.directory']

dir(database)


