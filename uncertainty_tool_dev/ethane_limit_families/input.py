# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary', 'BurkeH2O2'],
    reactionLibraries = ['BurkeH2O2inN2'],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = [
        'Disproportionation',
        'H_Abstraction',
        'intra_H_migration',
        'R_Recombination',
        'Intra_Disproportionation',
    ],
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

# Reaction systems
simpleReactor(
    temperature=(1350,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "ethane": 1.0,
    },
    terminationConversion={
        'ethane': 0.9,
    },
    terminationTime=(1e6,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
    filterReactions=True,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
