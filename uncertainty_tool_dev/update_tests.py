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

import os
import sys

sys.path.append('/home/moon/rmg/RMG-Py/test/rmgpy/data/kinetics/')

# +
import importlib

import kineticsTest
importlib.reload(kineticsTest)
# -

kineticsTest.setUpModule()

kineticsTest.database.kinetics.families

degen_test = kineticsTest.TestReactionDegeneracy()
degen_test.setup_class()
degen_test.test_r_addition_multiple_bond_benzene()
degen_test.test_degeneracy_same_reactant_different_resonance_structure()



# +

def test_degeneracy_same_reactant_different_resonance_structure(self):
    """Test if degeneracy is correct when reacting different resonance structures."""
    family_label = "Disproportionation"
    reactants = ["CC=C[CH2]", "CC=C[CH2]"]
    products = ["CC=CC", "C=CC=C"]

    correct_rxn_num = 1
    correct_degeneracy = {3}

    reaction_list = self.assert_correct_reaction_degeneracy(reactants, correct_rxn_num, correct_degeneracy, family_label, products)
    
    assert set(reaction_list[0].template) == {"C_rad/H2/Cd", "Cmethyl_Csrad/H/Cd"}
# -




