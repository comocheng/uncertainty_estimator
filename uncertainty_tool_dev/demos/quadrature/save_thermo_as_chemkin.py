# Notebook to grab a set of species we can use thermo libs to compare to GAV
# save the results in chemkin formats for later comparison

import rmgpy.species
import rmgpy.data.rmg
import rmgpy.chemkin


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


def species_in_list(sp, sp_list):
    return get_i_thing(sp, sp_list) != -1


run_variations = {
    'combined': [
        'NOx2018',
        'thermo_DFT_CCSDTF12_BAC',
        'DFT_QCI_thermo',
        'CBS_QB3_1dHR',
        'Klippenstein_Glarborg2016',
        'Elliott_OOQOOH',
    ],
    'NOx2018': ['NOx2018'],
    'thermo_DFT_CCSDTF12_BAC': ['thermo_DFT_CCSDTF12_BAC'],
    'DFT_QCI_thermo': ['DFT_QCI_thermo'],
    'CBS_QB3_1dHR': ['CBS_QB3_1dHR'],
    'Klippenstein_Glarborg2016': ['Klippenstein_Glarborg2016'],
    'Elliott_OOQOOH': ['Elliott_OOQOOH'],
    'GAV': [],
}

skip_list = [
    rmgpy.species.Species(smiles='[Ar]'),  # GAV fails for this
    rmgpy.species.Species(smiles='[He]'),  # GAV fails for this
    rmgpy.species.Species(smiles='[Ne]'),  # GAV fails for this
    rmgpy.species.Species(smiles='[H][H]'),  # GAV produces no comment for this
    rmgpy.species.Species(smiles='[C-]=[N+]=O'),  # GAV fails for this
    rmgpy.species.Species(smiles='[C]=N'),  # GAV fails on this
    rmgpy.species.Species(smiles='[N-]=[NH2+]'),  # GAV fails on this
]
for run, thermo_libraries in run_variations.items():
    print(f'Running with thermo libraries: {thermo_libraries}')
    if run == 'GAV':
        thermo_libraries = run_variations['combined']

    # Reset database
    database = rmgpy.data.rmg.RMGDatabase()
    database.load(
        path=rmgpy.settings['database.directory'],
        thermo_libraries=thermo_libraries,
        transport_libraries=[],
        reaction_libraries=[],
        seed_mechanisms=[],
        kinetics_families=[],
        kinetics_depositories=[],
        depository=False,
    )

    # Compile list of all species in this set of libraries
    all_species = []
    for i in range(len(thermo_libraries)):
        lib = thermo_libraries[i]
        for sp_key, entry in database.thermo.libraries[lib].entries.items():
            ref_sp = rmgpy.species.Species(molecule=[entry.item])
            if species_in_list(ref_sp, skip_list):
                continue
            if not species_in_list(ref_sp, all_species):
                ref_sp.label = ref_sp.molecule[0].get_formula() + f'({len(all_species)})'
                all_species.append(ref_sp)
    # print(len(all_species), 'species found')

    if run == 'GAV':
        database = rmgpy.data.rmg.RMGDatabase()
        database.load(
            path=rmgpy.settings['database.directory'],
            thermo_libraries=[],
            transport_libraries=[],
            reaction_libraries=[],
            seed_mechanisms=[],
            kinetics_families=[],
            kinetics_depositories=[],
            depository=False,
        )

    for i in range(len(all_species)):
        all_species[i].thermo = database.thermo.get_thermo_data(all_species[i])
        if not isinstance(all_species[i].thermo, rmgpy.thermo.NASA):
            all_species[i].thermo = all_species[i].thermo.to_nasa(Tmin=298, Tmax=3000, Tint=1000)

        if run == 'GAV':
            assert 'group additivity' in all_species[i].thermo.comment

    # save the species in a chemkin file
    rmgpy.chemkin.save_chemkin_file(f'thermo_chemkin_files/{run}.inp', all_species, [])
    rmgpy.chemkin.save_species_dictionary(f'thermo_chemkin_files/{run}_species_dictionary.txt', all_species)
