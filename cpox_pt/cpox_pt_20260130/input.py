# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo'],
    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006_adjusted', True), 'BurkeH2O2inArHe'],
    seedMechanisms = ['O2_dissociative_adsorption'],
    kineticsDepositories = ['training'],
    kineticsFamilies =[
        'default',
        'Surface_Adsorption_Single',
        'Surface_Adsorption_vdW',
        'Surface_Adsorption_Dissociative',
        'Surface_Dissociation',
        'Surface_Abstraction',
        'Surface_Dissociation_Double_vdW',
        'Surface_Dissociation_vdW',
        'Surface_Abstraction_vdW',
        'Surface_Dissociation_Beta',
        'Surface_Adsorption_Bidentate',
        'Surface_Bidentate_Dissociation',
        'Surface_Monodentate_to_Bidentate',
        'Surface_Dissociation_to_Bidentate', 
        'Surface_vdW_to_Bidentate',
        'Surface_Adsorption_Dissociative_Double',
        'Surface_Abstraction_Beta',
        # 'Surface_Abstraction_Beta_double_vdW',
        'Surface_Dissociation_Double',
        'Surface_Dissociation_Beta_vdW',
        'Surface_Abstraction_Beta_vdW',
        'Surface_Abstraction_Single_vdW',
    ],
    kineticsEstimator = 'rate rules',
)

catalystProperties(
    metal = 'Pt111'
)

# List of species
species(
    label='X',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
    label='CH4',
    reactive=True,
    structure=SMILES("[CH4]"),
)
species(
    label='O2',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u1 p2 c0 {2,S}
    2 O u1 p2 c0 {1,S}
    """),
)

species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]"),
)

species(
    label='CO2',
    reactive=True,
    structure=SMILES("O=C=O"),
)

species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
)

species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)

species(
    label='C2H6',
    reactive=True,
    structure=SMILES("CC"),
)

species(
    label='CH2O',
    reactive=True,
    structure=SMILES("C=O"),
)

species(
    label='CH3',
    reactive=True,
    structure=SMILES("[CH3]"),
)

species(
    label='C3H8',
    reactive=True,
    structure=SMILES("CCC"),
)

species(
    label='H',
    reactive=True,
    structure=SMILES("[H]"),
)

species(
    label='C2H5',
    reactive=True,
    structure=SMILES("C[CH2]"),
)

species(
    label='CH3OH',
    reactive=True,
    structure=SMILES("CO"),
)

species(
    label='HCO',
    reactive=True,
    structure=SMILES("[CH]=O"),
)

species(
    label='CH3CHO',
    reactive=True,
    structure=SMILES("CC=O"),
)

species(
    label='OH',
    reactive=True,
    structure=SMILES("[OH]"),
)

species(
    label='C2H4',
    reactive=True,
    structure=SMILES("C=C"),
)

species(
    label='CH3CH',
    reactive=True,
    structure=SMILES("[CH]C"),
)

species(
    label='CH3OO',
    reactive=True,
    structure=SMILES("CO[O]"),
)

species(
    label='HX',
    reactive=True,
    structure=adjacencyList(
    """
    1 H u0 p0 c0 {2,S}
    2 X u0 p0 c0 {1,S}
    """),
)

species(
    label='CO2X',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {3,D}
    2 O u0 p2 c0 {3,D}
    3 C u0 p0 c0 {1,D} {2,D}
    4 X u0 p0 c0
    """),
)

species(
    label='COX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,D}
    2 C u0 p0 c0 {1,D} {3,D}
    3 X u0 p0 c0 {2,D}
    """),
)

species(
    label='CH4X',
    reactive=True,
    structure=adjacencyList(
    """
    1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}
    6 X u0 p0 c0
    """),
)

species(
    label='OX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,D}
    2 X u0 p0 c0 {1,D}
    """),
)

species(
    label='CH2X',
    reactive=True,
    structure=adjacencyList(
    """
    1 C u0 p0 c0 {2,S} {3,S} {4,D}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 X u0 p0 c0 {1,D}
    """),
)

species(
    label='CH3X',
    reactive=True,
    structure=adjacencyList(
    """
    1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 X u0 p0 c0 {1,S}
    """),
)

species(
    label='CHX',
    reactive=True,
    structure=adjacencyList(
    """
    1 C u0 p0 c0 {2,S} {3,T}
    2 H u0 p0 c0 {1,S}
    3 X u0 p0 c0 {1,T}
    """),
)

species(
    label='CX',
    reactive=True,
    structure=adjacencyList(
    """
    1 C u0 p0 c0 {2,Q}
    2 X u0 p0 c0 {1,Q}
    """),
)

species(
    label='H2X',
    reactive=True,
    structure=adjacencyList(
    """
    1 H u0 p0 c0 {2,S}
    2 H u0 p0 c0 {1,S}
    3 X u0 p0 c0
    """),
)

species(
    label='OHX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,S} {3,S}
    2 H u0 p0 c0 {1,S}
    3 X u0 p0 c0 {1,S}
    """),
)

species(
    label='H2OX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,S} {3,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 X u0 p0 c0
    """),
)

species(
    label='CHOX',
    reactive=True,
    structure=adjacencyList(
    """
    1 O u0 p2 c0 {2,D}
    2 C u0 p0 c0 {1,D} {3,S} {4,S}
    3 H u0 p0 c0 {2,S}
    4 X u0 p0 c0 {2,S}
    """),
)

#----------
# Reaction systems
surfaceReactor(
    temperature=(600,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CH4": 0.041866,
        "O2": 0.03488,
        "Ar": 0.131246,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CH4":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=(600,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CH4": 0.108574,
        "O2": 0.02088,
        "Ar": 0.78547,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CH4":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=(2000,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CH4": 0.041866,
        "O2": 0.03488,
        "Ar": 0.131246,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CH4":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=(2000,'K'),
    initialPressure=(1.0, 'bar'),
    initialGasMoleFractions={
        "CH4": 0.108574,
        "O2": 0.02088,
        "Ar": 0.78547,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CH4":0.95,},
    terminationTime=(10., 's'),
#    terminationConversion={'O2': 0.99,},
    terminationRateRatio=0.01
)

simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=1e8,
    maximumEdgeSpecies=500000,
    maxNumSpecies=50,  # stop after 50 species
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)

generatedSpeciesConstraints(
    allowed=['input species', 'reaction libraries', 'seed mechanisms'],
    maximumSurfaceSites=2,
    maximumCarbonAtoms=2,
    maximumOxygenAtoms=2,
    maximumSurfaceBondOrder=4,
)




