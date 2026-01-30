#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################


# import pytest
import os

import rmgpy.rmg.input  # note this method is currently broken for multiple reactors in an input file
import rmgpy.rmg.main

import rmgpy.tools.uncertainty
import rmgpy.kinetics.uncertainties
import logging


import matplotlib.pyplot as plt

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

CH4 = rmgpy.species.Species().from_smiles('C')
CH3 = rmgpy.species.Species().from_smiles('[CH3]')
C2H4 = rmgpy.species.Species().from_smiles('C=C')


KineticUncertainty = rmgpy.tools.uncertainty.KineticParameterUncertainty()
ThermoUncertainty = rmgpy.tools.uncertainty.ThermoParameterUncertainty()

def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    raise AssertionError(f'Could not find {thing} in list of {len(thing_list)} things')


class TestMech00Ethane():

    def setup_class(self):
        """A function that is run ONCE before all unit tests in this class."""

        # Must use annotated chemkin file
        chemkin_file = './mech_00_ethane/chemkin/chem_annotated.inp'
        mech_dir = os.path.dirname(chemkin_file)
        dict_file = os.path.join(mech_dir, 'species_dictionary.txt')
        input_file = os.path.join(os.path.dirname(mech_dir), 'input.py')

        # Run Gao estimation of input parameters (takes a long time to load database)
        logging.info('Loading Mech 00 Ethane for uncertainty analysis')
        self.uncertainty = rmgpy.tools.uncertainty.Uncertainty(
            output_directory=os.path.join(mech_dir, 'uncertainty')
        )
        self.uncertainty.load_model(chemkin_file, dict_file)

        global rmg0
        rmg0 = rmgpy.rmg.main.RMG()
        rmg0.output_directory = '.'
        rmgpy.rmg.input.read_input_file(input_file, rmg0)

        thermo_libs = rmg0.thermo_libraries
        kinetics_libs = [x[0] if isinstance(x, tuple) else x for x in rmg0.reaction_libraries]
        kinetics_families = rmg0.kinetics_families

        logging.info(f'Thermo libraries: {thermo_libs}')
        logging.info(f'Kinetics libraries: {kinetics_libs}')
        logging.info(f'Kinetics families: {kinetics_families}')

        assert set(thermo_libs) == set(['primaryThermoLibrary'])
        assert set(kinetics_libs) == set(['BurkeH2O2inN2'])
        assert kinetics_families == 'default'

        logging.info('Loading database for Mech 00 Ethane')
        self.uncertainty.load_database(  # this call does the averaging up of kinetics families
            thermo_libraries=thermo_libs,
            kinetics_families=kinetics_families,
            reaction_libraries=kinetics_libs,
        )

        # make sure CH3 hasn't been added to primaryThermoLibrary
        thermo_lib_item_list = [self.uncertainty.database.thermo.libraries['primaryThermoLibrary'].entries[key].item \
                                    for key in self.uncertainty.database.thermo.libraries['primaryThermoLibrary'].entries.keys()]
        assert not any([CH3.is_isomorphic(item) for item in thermo_lib_item_list]), "You need a new test if CH3 gets into the primaryThermoLibrary"


    def test_uncorrelated_species(self):
        # Get the different kinetic and thermo sources
        self.uncertainty.extract_sources_from_model()
        self.uncertainty.assign_parameter_uncertainties()
        self.uncertainty.compile_all_sources()

        i_CH4 = get_i_thing(CH4, self.uncertainty.species_list)
        i_CH3 = get_i_thing(CH3, self.uncertainty.species_list)
        i_C2H4 = get_i_thing(C2H4, self.uncertainty.species_list)

        # check that the sources makes sense
        CH4_source_entry = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_CH4]]
        CH3_source_entry = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_CH3]]
        C2H4_source_entry = self.uncertainty.species_sources_dict[self.uncertainty.species_list[i_C2H4]]

        assert set(['Library']) == set(CH4_source_entry.keys())
        assert set(['Library', 'GAV']) == set(CH3_source_entry.keys())
        assert set(['GAV']) == set(C2H4_source_entry.keys())

        # check that the calculation makes sense
        assert self.uncertainty.thermo_input_uncertainties[i_CH4] == ThermoUncertainty.dG_library


if __name__ == '__main__':
    print('Running TestMech00Ethane')
    a = TestMech00Ethane()
    a.setup_class()
    a.test_uncorrelated_species()