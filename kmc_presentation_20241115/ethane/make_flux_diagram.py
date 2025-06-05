import os
import rmgpy.tools.fluxdiagram 

basedir = '/home/moon/uncertainty_estimator/kmc_presentation_20241115/ethane'
spec_dict = os.path.join(basedir, 'chemkin', 'species_dictionary.txt')
chemkin = os.path.join(basedir, 'chemkin', 'chem_annotated.inp')
input_file = os.path.join(basedir, 'input.py')
output_path = os.path.join(basedir, 'flux_diagram')
rmgpy.tools.fluxdiagram.create_flux_diagram(input_file, chemkin, spec_dict, save_path=output_path)



