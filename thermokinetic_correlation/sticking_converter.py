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

# example script to handle sticking coefficients
import pickle
import numpy as np

# define constants
R = 8.314472
kB = 1.3806504e-23

# +
# Read in data files
species_dict_file = 'saved_data_files_v2.4/species_dict.pickle'
with open(species_dict_file, 'rb') as f:
    my_thermo_data = pickle.load(f)

reaction_dict_file = 'saved_data_files_v2.4/reaction_dict.pickle'
with open(reaction_dict_file, 'rb') as f:
    my_kinetic_data = pickle.load(f)

# -

my_example_reaction = my_kinetic_data['H2 + X_3 + X_3 <=> HX_3 + HX_3']

my_example_reaction

my_thermo_data['HX_3']


def get_thermo_from_NASA(NASA0, NASA1, T):
    # compute thermo properties from nasa polynomials  units are Joules and mols
    # NASA0 is the lower temperature range and NASA1 is the higher
    # expecting NASA polynomials in the following dictionary format:
#     {'coeffs': array([ 3.53732118e+00, -1.21570202e-03,  5.31615358e-06, -4.89440364e-09,
#          1.45843807e-12, -1.03858843e+03,  4.68368633e+00]),
#      'Tmin': (100,'K'),
#      'Tmax': (1074.56,'K')}
    
    assert T >= NASA0['Tmin']
    assert T <= NASA1['Tmax']
    
    a_low = NASA0['coeffs']
    a_high = NASA1['coeffs']
    
    if T < NASA0['Tmax']:
        cp = a_low[0] + a_low[1] * T + a_low[2] * T**2.0 + a_low[3] * T**3.0 + a_low[4] * T**4.0
        h = a_low[0] * T + a_low[1] / 2.0 * T**2.0 + a_low[2] / 3.0 * T**3.0 + a_low[3] / 4.0 * T**4.0 + a_low[4] / 5.0 * T**5.0 + a_low[5]
        s = a_low[0] * np.log(T) + a_low[1] * T + a_low[2] / 2.0 * T**2.0 + a_low[3] / 3.0 * T**3.0 + a_low[4] / 4.0 * T**4.0 + a_low[6]
    else:
        cp = a_high[0] + a_high[1] * T + a_high[2] * T**2.0 + a_high[3] * T**3.0 + a_high[4] * T**4.0
        h = a_high[0] * T + a_high[1] / 2.0 * T**2.0 + a_high[2] / 3.0 * T**3.0 + a_high[3] / 4.0 * T**4.0 + a_high[4] / 5.0 * T**5.0 + a_high[5]
        s = a_high[0] * np.log(T) + a_high[1] * T + a_high[2] / 2.0 * T**2.0 + a_high[3] / 3.0 * T**3.0 + a_high[4] / 4.0 * T**4.0 + a_high[6]

    
    cp *= R
    h *= R
    s *= R

    return cp, h, s


def get_enthalpy_of_reaction(reaction_entry, Tref=1000):
    Hrxn = 0
    for product in reaction_entry['products']:
        NASA0 = my_thermo_data[product]['NASA0']
        NASA1 = my_thermo_data[product]['NASA1']
        cp, h, s = get_thermo_from_NASA(NASA0, NASA1, Tref)  # in J/mol
        Hrxn += h
    for reactant in reaction_entry['reactants']:
        NASA0 = my_thermo_data[reactant]['NASA0']
        NASA1 = my_thermo_data[reactant]['NASA1']
        cp, h, s = get_thermo_from_NASA(NASA0, NASA1, Tref)  # in J/mol
        Hrxn -= h
    return Hrxn
get_enthalpy_of_reaction(my_kinetic_data['H2 + X_3 + X_3 <=> HX_3 + HX_3'])


def get_activation_energy(reaction_entry, Hrxn):
    alpha = reaction_entry['kinetics']['alpha']
    E0 = reaction_entry['kinetics']['E0'] * 1000
    assert reaction_entry['kinetics']['E0_units'] == 'kJ/mol'
    Ea = alpha * Hrxn + E0
    
    if E0 > 0:
        if Hrxn < 0.0 and Ea < 0.0:  # Activation barrier shouldn't be negative if intrinsic barrier E0 is > 0
            Ea = 0.0
        elif Hrxn > 0.0 and Ea < Hrxn:  # Raise activation barrier so reverse won't be negative
            Ea = Hrxn
    
    return Ea  # J/mol


Hrxn = get_enthalpy_of_reaction(my_example_reaction, Tref=1000)
get_activation_energy(my_example_reaction, Hrxn)



my_kinetic_data['H2 + X_3 + X_3 <=> HX_3 + HX_3']['kinetics']


# # Dimensionless sticking coefficient: $\gamma = a T^b e^{-c/RT}$

def get_sticking_coefficient(reaction_entry, T, Tref=None):
    """
    Return the sticking coefficient (dimensionless) at
    temperature `T` in K and enthalpy of reaction `dHrxn` in J/mol.
    Not supposed to exceed 1.0.
    
    Tref is the temperature at which we calculate H_rxn
    """
    
    if 'stickingcoefficientbep' in reaction_entry['parameterization'].lower():
        if Tref is None:
            Tref = T
        Hrxn = get_enthalpy_of_reaction(my_example_reaction, Tref=Tref)
        Ea = get_activation_energy(reaction_entry, Hrxn)
    else:
        Ea = reaction_entry['kinetics']['Ea']
    A = reaction_entry['kinetics']['A']
    n = reaction_entry['kinetics']['n']
    
    sticking_coefficient = A * T ** n * np.exp(-Ea / (R * T))
    assert sticking_coefficient >= 0
    return min(sticking_coefficient, 1.0)


get_sticking_coefficient(my_example_reaction, 1000)



# # Reaction Rate: $k_f = \frac{\gamma}{\Gamma_{tot}^m}\sqrt{\frac{RT}{2\pi W}}$

def get_sticking_rate_coefficient(reaction_entry, T, surface_site_density):
    """
    Return the overall surface rate coefficient for the forward reaction at
    temperature `T` in K with surface site density `surface_site_density` in mol/m2.
    Value is returned in combination of [m,mol,s]
    """

    rate_coefficient = get_sticking_coefficient(reaction_entry, T)
    
    # detect which species is the gas-phase adsorbate
    adsorbate = None
    for r in reaction_entry['reactants']:
        reactant = my_thermo_data[r]
        if reactant['is_surface_species']:
            rate_coefficient /= surface_site_density
            sites = reactant['n_surface_sites']
            if sites > 1:
                rate_coefficient /= sites
        else:
            if adsorbate is None:
                adsorbate = reactant
            else:
                raise ValueError("More than one adsorbate detected")
    
    rate_coefficient *= np.sqrt(kB * T / (2 * np.pi * adsorbate['molecular_weight_kg']))
    
    # Multidentate adsorption requires multiplication of the sticking coefficient
    # with the number of binding sites**stoichiometric coefficients (it is 1 for monodentates)
    # It was already integrated in the loop above for the reactants
    for p in reaction_entry['products']:
        product = my_thermo_data[p]
        sites = product['n_surface_sites']
        if sites > 1:
            rate_coefficient *= sites
    return rate_coefficient


surface_site_density = 2.7200E-05  # mol/m^2
get_sticking_rate_coefficient(my_example_reaction, 1000, surface_site_density)


# # Figure out the units on the rate coefficient

# +
def write_units_str(units_dict):
    units_str = ''
    
    num_items = []
    den_items = []
    for unit, power in units_dict.items():
        if power > 1:
            num_items.append(f'{unit}^{power}')
        elif power == 1:
            num_items.append(f'{unit}')
        elif power == -1:
            den_items.append(f'{unit}')
        elif power < -1:
            den_items.append(f'{unit}^{abs(power)}')

    num_string = '*'.join(num_items)
    den_string = '*'.join(den_items)
    if len(den_items) > 1:
        den_string = f'({den_string})'
    
    if len(num_items) == 0 and len(den_items) == 0:
        return ''
    
    elif len(num_items) > 0 and len(den_items) == 0:
        return num_string
    elif len(num_items) == 0 and len(den_items) > 0:
        return f'1/{den_string}'
    else:
        return f'{num_string}/{den_string}'
    
def get_species_units(species):
    if species['is_surface_species']:
        return {
            'mol': 1,
            'm': -2,
            's': 0,
        }
    return {
        'mol': 1,
        'm': -3,
        's': 0,
    }

def get_forward_rate_units(reaction_entry, ref_species=None):
    if ref_species is None:
        reactant_label = reaction_entry['reactants'][0]
        ref_species = my_thermo_data[reactant_label]
    else:
        assert ref_species['label'] in reaction_entry['reactants'], 'reference species not in reaction'
    
    # make the equation
    units_dref_dt = get_species_units(ref_species)
    units_dref_dt['s'] = units_dref_dt['s'] - 1
    
#     print('ref units: ', write_units_str(units_dref_dt))
    
    forward_eq_units = {
        'mol': 0,
        'm': 0,
        's': 0,
    }
    
    for i in range(len(reaction_entry['reactants'])):
        reactant = reaction_entry['reactants'][i]
        sp_units = get_species_units(my_thermo_data[reactant])
        forward_eq_units['mol'] += sp_units['mol']
        forward_eq_units['m'] += sp_units['m']
        forward_eq_units['s'] += sp_units['s']
#     print('eq units: ', write_units_str(forward_eq_units))
    
    k_units = {
        'mol': units_dref_dt['mol'] - forward_eq_units['mol'],
        'm': units_dref_dt['m'] - forward_eq_units['m'],
        's': units_dref_dt['s'] - forward_eq_units['s'],
    }
    return write_units_str(k_units)


# -

get_forward_rate_units(my_example_reaction)








