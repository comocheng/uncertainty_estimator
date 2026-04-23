#!/usr/bin/env python
# encoding: utf-8


name = "shortened_surfaceThermoPt111"
shortDesc = u"Surface adsorbates on Pt(111)"
longDesc = u"""
Surface species adsorbed on Pt(111). The thermochemistry of all adsorbates with up to 2 heavy atoms was calculated by Katrin Blondal at Brown University around 2018,
based on DFT calculations by Jelena Jelic at KIT. See https://doi.org/10.1021/acs.iecr.9b01464 for the details on the computational methods as well as the results.
This database was extended with DFT calculations for larger adsorbates by Bjarne Kreitz (Brown University).  
The computational methods for the extension are explained in detail in https://doi.org/10.1021/acscatal.2c03378. If you use this database in your work, please cite the publications mentioned above. 
Note: X indicates a bond to the surface. It is always on the left hand site of an atom that is bonded to the surface e.g. XCCH2 it means that C is bonded to the surface.
If the X is on the right hand side and at the end of a label, it means that this species is physisorbed. 
"""

entry(
    index = 1,
    label = "vacant",
    molecule =
"""
1 X  u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(1000.0,'K'), Tmax=(3000.0, 'K')),
        ],
        Tmin = (298.0, 'K'),
        Tmax = (3000.0, 'K'),
    ),
    metal = "Pt",
    facet = "111",
)

entry(
    index = 2,
    label = "XH",
    molecule =
"""
1 X  u0 p0 c0 {2,S}
2 H  u0 p0 c0 {1,S}
""",
    thermo = NASA(
       polynomials = [
        NASAPolynomial(coeffs=[-2.07570125E+00, 1.73580835E-02, -2.60920784E-05, 1.89282268E-08, -5.38835643E-12, -3.16618959E+03, 8.15361518E+00], Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
        NASAPolynomial(coeffs=[2.72248139E+00, -1.06817206E-03, 1.98653790E-06, -1.12048461E-09, 2.09811636E-13, -4.21823896E+03, -1.53207470E+01], Tmin=(1000.0, 'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin = (298.0,'K'),
        Tmax = (2000.0,'K'),
    ),
    longDesc = u"""Calculated by Bjarne Kreitz at Brown University using statistical mechanics (file: ThermoPt111.py).
            Based on DFT calculations by Bjarne Kreitz from Brown University. DFT calculations were performed with Quantum Espresso
            using PAW pseudopotentials and the BEEF-vdW functional for an optimized 3x3 supercell (1/9ML coverage)
            following the procedure outlined by Blondal et al (DOI:10.1021/acs.iecr.9b01464). The following settings were applied:
            kpoints=(5x5x1), 4 layers (2 bottom layers fixed), ecutwfc=60 Ry, smearing='mazari-vanderbilt', mixing_mode='local-TF',
            fmax=2.5e-2. DFT binding energy: -2.481 eV.
""",
    metal = "Pt",
    facet = "111",
)

entry(
    index = 5,
    label = "XOH",
    molecule =
"""
1 X  u0 p0 c0 {2,S}
2 O  u0 p2 c0 {1,S} {3,S}
3 H  u0 p0 c0 {2,S}
""",
    thermo=NASA(
        polynomials=[
            NASAPolynomial(coeffs=[1.42359828E+00, 1.57830409E-02, -2.91659307E-05, 2.50432531E-08, -8.04088046E-12,
                                   -1.89993029E+04, -3.15230478E+00], Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[5.03574444E+00, -1.34422013E-03, 2.25915779E-06, -1.08547653E-09, 1.77875516E-13,
                                   -1.96344169E+04, -2.00345222E+01], Tmin=(1000.0, 'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin=(298.0, 'K'),
        Tmax=(2000.0, 'K'),
    ),
    longDesc=u"""Calculated by Bjarne Kreitz at Brown University using statistical mechanics (file: ThermoPt111.py).
                Based on DFT calculations by Bjarne Kreitz from Brown University. DFT calculations were performed with Quantum Espresso
                using PAW pseudopotentials and the BEEF-vdW functional for an optimized 3x3 supercell (1/9ML coverage)
                following the procedure outlined by Blondal et al (DOI:10.1021/acs.iecr.9b01464). The following settings were applied:
                kpoints=(5x5x1), 4 layers (2 bottom layers fixed), ecutwfc=60 Ry, smearing='mazari-vanderbilt', mixing_mode='local-TF',
                fmax=2.5e-2. DFT binding energy: -1.628 eV.
    """,
    metal="Pt",
    facet="111",
)

entry(
    index = 9,
    label = "XO",
    molecule =
"""
1 X  u0 p0 c0 {2,D}
2 O  u0 p2 c0 {1,D}
""",
    thermo=NASA(
        polynomials=[
            NASAPolynomial(coeffs=[-2.94475701E-01, 1.44162624E-02, -2.61322704E-05, 2.19005957E-08, -6.98019420E-12,
                                   -1.64619234E+04, -1.99445623E-01], Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[2.90244691E+00, -3.38584457E-04, 6.43372619E-07, -3.66326660E-10, 6.90093884E-14,
                                   -1.70497471E+04, -1.52559728E+01], Tmin=(1000.0, 'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin=(298.0, 'K'),
        Tmax=(2000.0, 'K'),
    ),
    longDesc=u"""Calculated by Bjarne Kreitz at Brown University using statistical mechanics (file: ThermoPt111.py).
                Based on DFT calculations by Bjarne Kreitz from Brown University. DFT calculations were performed with Quantum Espresso
                using PAW pseudopotentials and the BEEF-vdW functional for an optimized 3x3 supercell (1/9ML coverage)
                following the procedure outlined by Blondal et al (DOI:10.1021/acs.iecr.9b01464). The following settings were applied:
                kpoints=(5x5x1), 4 layers (2 bottom layers fixed), ecutwfc=60 Ry, smearing='mazari-vanderbilt', mixing_mode='local-TF',
                fmax=2.5e-2. DFT binding energy: -4.101 eV.
    """,
    metal="Pt",
    facet="111",
)



entry(
    index = 38,
    label = "XCH",
    molecule =
"""
1 X  u0 p0 c0 {2,T}
2 C  u0 p0 c0 {1,T} {3,S}
3 H  u0 p0 c0 {2,S}
""",
    thermo=NASA(
        polynomials=[
            NASAPolynomial(coeffs=[-2.66805027E+00, 2.90693485E-02, -4.82653552E-05, 3.87589256E-08, -1.19749384E-11,
                                   -2.91815541E+03, 9.72941427E+00], Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[4.90429859E+00, -2.63865042E-03, 4.71729273E-06, -2.51267002E-09, 4.49659283E-13,
                                   -4.46440812E+03, -2.67107945E+01], Tmin=(1000.0, 'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin=(298.0, 'K'),
        Tmax=(2000.0, 'K'),
    ),
    longDesc=u"""Calculated by Bjarne Kreitz at Brown University using statistical mechanics (file: ThermoPt111.py).
                Based on DFT calculations by Bjarne Kreitz from Brown University. DFT calculations were performed with Quantum Espresso
                using PAW pseudopotentials and the BEEF-vdW functional for an optimized 3x3 supercell (1/9ML coverage)
                following the procedure outlined by Blondal et al (DOI:10.1021/acs.iecr.9b01464). The following settings were applied:
                kpoints=(5x5x1), 4 layers (2 bottom layers fixed), ecutwfc=60 Ry, smearing='mazari-vanderbilt', mixing_mode='local-TF',
                fmax=2.5e-2. DFT binding energy: -4.820 eV.
    """,
    metal="Pt",
    facet="111",
)

entry(
    index = 43,
    label = "XCH3",
    molecule =
"""
1 X  u0 p0 c0 {2,S}
2 C  u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H  u0 p0 c0 {2,S}
4 H  u0 p0 c0 {2,S}
5 H  u0 p0 c0 {2,S}
""",
    thermo=NASA(
        polynomials=[
            NASAPolynomial(coeffs=[-4.44549060E-02, 1.94367665E-02, -1.91028733E-05, 1.11269371E-08, -2.73735895E-12,
                                   -6.38803762E+03, -1.73375969E-01], Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[8.65704524E+00, -7.90307778E-03, 1.40100438E-05, -7.40016074E-09, 1.31516592E-12,
                                   -8.63598517E+03, -4.43352558E+01], Tmin=(1000.0, 'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin=(298.0, 'K'),
        Tmax=(2000.0, 'K'),
    ),
    longDesc=u"""Calculated by Bjarne Kreitz at Brown University using statistical mechanics (file: ThermoPt111.py).
                Based on DFT calculations by Bjarne Kreitz from Brown University. DFT calculations were performed with Quantum Espresso
                using PAW pseudopotentials and the BEEF-vdW functional for an optimized 3x3 supercell (1/9ML coverage)
                following the procedure outlined by Blondal et al (DOI:10.1021/acs.iecr.9b01464). The following settings were applied:
                kpoints=(5x5x1), 4 layers (2 bottom layers fixed), ecutwfc=60 Ry, smearing='mazari-vanderbilt', mixing_mode='local-TF',
                fmax=2.5e-2. DFT binding energy: -1.721 eV.
    """,
    metal="Pt",
    facet="111",
)

entry(
    index = 49,
    label = "XCO",
    molecule =
"""
1 X  u0  p0 c0 {2,D}
2 C  u0  p0 c0 {1,D} {3,D}
3 O  u0  p2 c0 {2,D}
""",
    thermo=NASA(
        polynomials=[
            NASAPolynomial(coeffs=[1.42895000E+00, 1.40374509E-02, -2.21178920E-05, 1.78659581E-08, -5.71478802E-12,
                                   -3.45688484E+04, -7.78265517E+00], Tmin=(298.0, 'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[5.48656312E+00, -1.68118543E-03, 3.09030310E-06, -1.71186643E-09, 3.15864598E-13,
                                   -3.54815495E+04, -2.76788365E+01], Tmin=(1000.0, 'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin=(298.0, 'K'),
        Tmax=(2000.0, 'K'),
    ),
    longDesc=u"""Calculated by Bjarne Kreitz at Brown University using statistical mechanics (file: ThermoPt111.py).
                Based on DFT calculations by Bjarne Kreitz from Brown University. DFT calculations were performed with Quantum Espresso
                using PAW pseudopotentials and the BEEF-vdW functional for an optimized 3x3 supercell (1/9ML coverage)
                following the procedure outlined by Blondal et al (DOI:10.1021/acs.iecr.9b01464). The following settings were applied:
                kpoints=(5x5x1), 4 layers (2 bottom layers fixed), ecutwfc=60 Ry, smearing='mazari-vanderbilt', mixing_mode='local-TF',
                fmax=2.5e-2. DFT binding energy: -1.415 eV.
    """,
    metal="Pt",
    facet="111",
)
