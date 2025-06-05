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
import copy
import rmgpy.reaction
import rmgpy.kinetics
import rmgpy.chemkin
import numpy as np

import matplotlib.pyplot as plt
# %matplotlib inline
# -



kinetics = rmgpy.kinetics.surface.StickingCoefficientBEP(
        A = 0.1,
        n = 0,
        alpha = 0.69,
        E0 = (107.9, 'kJ/mol'),
        Tmin = (200, 'K'),
        Tmax = (3000, 'K'),
    )

kinetics.get_rate_coefficient(1000)



# +
# load example 

gas_chemkin = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241020/chem_annotated-gas.inp'
gas_surface = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241020/chem_annotated-surface.inp'
sp_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_pt_20241020/species_dictionary.txt'
test_species_list, test_reaction_list = rmgpy.chemkin.load_chemkin_file(gas_chemkin, sp_dict, surface_path=gas_surface)
# -





def sticking2arrhenius(my_sticking_reaction):
    SDEN = 2.7200E-05  # mol/m^2
    my_arrhenius_reaction = copy.deepcopy(my_sticking_reaction)
    my_arrhenius_reaction.kinetics = rmgpy.kinetics.surface.SurfaceArrhenius()
    
    Tmin = my_sticking_reaction.kinetics.Tmin
    if Tmin is None:
        Tmin = 200
    Tmax = my_sticking_reaction.kinetics.Tmax
    if Tmax is None:
        Tmax = 3000

    Tlist = np.linspace(Tmin, Tmax, 51)
    klist = np.zeros(len(Tlist))
    for i, T in enumerate(Tlist):
        klist[i] = my_sticking_reaction.get_rate_coefficient(T, surface_site_density=SDEN)
    
    kunits = get_forward_rate_units(my_sticking_reaction)
    
    my_arrhenius_reaction.kinetics.fit_to_data(Tlist, klist, kunits)
    return my_arrhenius_reaction


new_reaction = sticking2arrhenius(test_reaction_list[108])

new_reaction.get_rate_coefficient(800)

test_reaction_list[108].get_rate_coefficient(800, surface_site_density=SDEN)

test_reaction_list[108]

test_reaction_list[0]



my_dict = {
            'mol': 1,
            'm': 2,
            's': 0,
        }

for key, item in my_dict.items():
    print(item)



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
    if species.contains_surface_site():
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

def get_forward_rate_units(my_reaction, ref_species=None):
    if ref_species is None:
        ref_species = my_reaction.reactants[0]
    else:
        assert ref_species in my_reaction.reactants
    
    # make the equation
    units_dref_dt = get_species_units(ref_species)
    units_dref_dt['s'] = units_dref_dt['s'] - 1
    
#     print('ref units: ', write_units_str(units_dref_dt))
    
    forward_eq_units = {
        'mol': 0,
        'm': 0,
        's': 0,
    }
    
    for i in range(len(my_reaction.reactants)):
        sp_units = get_species_units(my_reaction.reactants[i])
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

get_forward_rate_units(test_reaction_list[108])



write_units_str(get_species_units(test_reaction_list[108].reactants[0]))

test_reaction_list[108].reactants[2]

Tmin = 200 if test_reaction_list[108].kinetics.Tmin is None else test_reaction_list[108].kinetics.Tmin

Tmin

test_reaction_list[108].kinetics

test_reaction_list[108].get_rate_coefficient(1000, surface_site_density=SDEN)

test_reaction_list[108].get_uni

test_reaction_list[108].kinetics



# +
# SDEN
# -

SDEN = 2.7200E-05  # mol/m^2

test_reaction_list[108].get_rate_coefficient(1000, surface_site_density=SDEN)

test_reaction_list[108].kinetics.Tmin

test_reaction_list[108].kinetics.Tmax

print([-1 for x in test_reaction_list[108].reactants].extend([1 for x in test_reaction_list[108].products]))

[-1 for x in test_reaction_list[108].reactants].extend([1 for x in test_reaction_list[108].products])

[1, 2, 3].extend([1, 2])

my_list = [1, 2, 3]
my_list.extend([1, 2])
my_list

new_arrhenius_kinetics = rmgpy.kinetics.surface.SurfaceArrhenius()

new_arrhenius_kinetics.fit_to_data()



test_reaction_list[108].kinetics





test_reaction_list[108].get_rate_coefficient(1000, surface_site_density=2.2e-9)

for i in range(len(test_reaction_list)):
    print(i, type(test_reaction_list[i].kinetics))





def get_thermo_from_NASA(NASA0, NASA1, T):
    # compute thermo properties from nasa polynomials
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

    R = 8.314472
    cp *= R
    h *= R
    s *= R

    return cp, h, s
