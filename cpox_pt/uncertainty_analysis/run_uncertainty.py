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
# script to run Gao uncertainty using Rh cpox mechanism

# +
import os

import rmgpy.tools.uncertainty

# +
chemkin_file = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin/chem_annotated-gas.inp'
dict_file = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin/species_dictionary.txt'
surface_path = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin/chem_annotated-surface.inp'

# Initialize the Uncertainty class instance and load the model
uncertainty = rmgpy.tools.uncertainty.Uncertainty(output_directory='rmg/uncertainty')
uncertainty.load_model(chemkin_file, dict_file, surface_path=surface_path)


# uncertainty.load_database(
#     thermo_libraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo'],
#     kinetics_families=[
#         'default',
#         'Surface_Abstraction',
#         'Surface_Abstraction_Single_vdW',
#         'Surface_Abstraction_vdW',
#         'Surface_Adsorption_Bidentate',
#         'Surface_Adsorption_Dissociative',
#         'Surface_Adsorption_Dissociative_Double',
#         'Surface_Adsorption_Double',
#         'Surface_Adsorption_Single',
#         'Surface_Adsorption_vdW',
#         'Surface_Bidentate_Dissociation',
#         'Surface_Dissociation',
#         'Surface_Dissociation_Beta',
#         'Surface_Dissociation_Double_vdW',
#         'Surface_Dissociation_vdW',
#         'Surface_Dissociation_to_Bidentate',
#         'Surface_EleyRideal_Addition_Multiple_Bond',
#         'Surface_Migration',
#         'Surface_vdW_to_Bidentate'
#     ],
#     reaction_libraries=['Surface/CPOX_Pt/Deutschmann2006_adjusted','BurkeH2O2inArHe'],
#     kinetics_depositories = ['training'],
# )

# uncertainty.extract_sources_from_model()
# uncertainty.assign_parameter_uncertainties()


# +
initial_mole_fractions = {
    uncertainty.species_list[0]: 0.131246,  # Ar
    uncertainty.species_list[4]: 0.03488,  # O2
    uncertainty.species_list[3]: 0.041866,  # CH4
}
initial_surface_coverages = {
    uncertainty.species_list[22]: 1.0
}

surface_volume_ratio = 1600  # Emily says she made this up. I respect that
surface_site_density = (0.7200E-09, 'mol/cm^2')

T = (1000, 'K')
P = (1.0, 'bar')
termination_time = (0.5, 'ms')

uncertainty.sensitivity_analysis(
    initial_mole_fractions,
    [uncertainty.species_list[3]],
    T,
    P,
    termination_time,
    sensitivity_threshold=1e-3,
    number=10,
    fileformat='.png',
    initial_surface_coverages=initial_surface_coverages,
    surface_volume_ratio=surface_volume_ratio,
    surface_site_density=surface_site_density
)
# -









# +
uncertainty.thermo_input_uncertainties
uncertainty.kinetic_input_uncertainties

for i in range(len(uncertainty.reaction_list)):
#     print(uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]])
    source_entry = uncertainty.reaction_sources_dict[uncertainty.reaction_list[i]]
    
    if 'Rate Rules' in source_entry.keys():
        if source_entry['Rate Rules'][1]['node']:
            print(source_entry)
            print()
        else:
            print(i, source_entry['Rate Rules'][0], source_entry['Rate Rules'][1]['template'])
# -

any([x.contains_surface_site() for x in uncertainty.species_list])

uncertainty.species_list[3]









uncertainty.species_list[27].contains_surface_site()

uncertainty.reaction_list[0].is_surface_reaction()

[x for x in uncertainty.reaction_list if not x.is_surface_reaction()]

[x for x in uncertainty.species_list if not x.contains_surface_site()]



uncertainty.reaction_sources_dict[uncertainty.reaction_list[100]]['Rate Rules'][1]['rules'][0]

uncertainty.reaction_sources_dict[uncertainty.reaction_list[86]]['Rate Rules']

uncertainty.reaction_sources_dict[uncertainty.reaction_list[104]]['Rate Rules']

# +
# sensitivity analysis

# Define the reaction conditions
initial_mole_fractions = {mapping[ethane]: 1.0}
T = (1300, 'K')
P = (1, 'bar')
termination_time = (0.5, 'ms')
sensitive_species=[mapping[ethane], mapping[C2H4]]
# -


