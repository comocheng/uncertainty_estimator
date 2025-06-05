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

import rmgpy.chemkin
import numpy as np
import os
import copy
import matplotlib.pyplot as plt
# %matplotlib inline


# +
emily_mech = './cpox_rh_emily'
e_chemkin = os.path.join(emily_mech, 'chem_annotated-gas.inp')
e_surface = os.path.join(emily_mech, 'chem_annotated-surface.inp')
e_dict = os.path.join(emily_mech, 'species_dictionary.txt')
e_sp_gas, e_rxn_gas = rmgpy.chemkin.load_chemkin_file(e_chemkin, e_dict, read_comments=False, use_chemkin_names=True)
e_sp_surf, e_rxn_surf = rmgpy.chemkin.load_chemkin_file(e_surface, e_dict, read_comments=False, use_chemkin_names=True)

# e_sp, e_rxn = rmgpy.chemkin.load_chemkin_file(e_chemkin, e_dict, surface_path=e_surface, read_comments=False)

sevy_mech = './cpox_rh_20241028'
s_chemkin = os.path.join(sevy_mech, 'chem_annotated-gas.inp')
s_surface = os.path.join(sevy_mech, 'chem_annotated-surface.inp')
s_dict = os.path.join(sevy_mech, 'species_dictionary.txt')
s_sp_gas, s_rxn_gas = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, use_chemkin_names=True)
s_sp_surf, s_rxn_surf = rmgpy.chemkin.load_chemkin_file(s_surface, s_dict, use_chemkin_names=True)
# s_sp, s_rxn = rmgpy.chemkin.load_chemkin_file(s_chemkin, s_dict, surface_path=s_surface)
# -

display(e_rxn_surf[0])

# +
# get_i_thing(e_rxn_surf[0].reactants[2], e_sp_rxn)

# +
# check which species are missing from the gas list
for i in range(len(e_rxn_surf)):
#     print(i)
    for reactant in e_rxn_surf[i].reactants:
        j = get_i_thing(reactant, e_sp_surf)
        k = get_i_thing(reactant, e_sp_gas)
        assert j >=0 or k>=0
    for product in e_rxn_surf[i].products:
        j = get_i_thing(product, e_sp_surf)
        k = get_i_thing(product, e_sp_gas)
        assert j >=0 or k>=0 
    

# -

for s in e_sp_gas:
    print(s)



def get_i_thing(my_thing, thing_list):
    for j in range(len(thing_list)):
        if thing_list[j].is_isomorphic(my_thing):
            return j
    return -1


# +
new_sp_list_gas = []
new_rxn_list_gas = []
new_sp_list_surf = []
new_rxn_list_surf = []

# pick the gases
for i in range(len(e_sp_gas)):
    j = get_i_thing(e_sp_gas[i], s_sp_gas)
    if j >= 0:
        new_sp_list_gas.append(s_sp_gas[j])
    else:
        display(e_sp_gas[i])
        
for i in range(len(e_rxn_gas)):
    j = get_i_thing(e_rxn_gas[i], s_rxn_gas)
    if j >= 0:
        new_rxn_list_gas.append(s_rxn_gas[j])
    else:
        display(e_rxn_gas[i])

# pick the surfaces
for i in range(len(e_sp_surf)):
    j = get_i_thing(e_sp_surf[i], s_sp_surf)
    if j >= 0:
        new_sp_list_surf.append(s_sp_surf[j])
    else:
        display(e_sp_surf[i])
        
for i in range(len(e_rxn_surf)):
    j = get_i_thing(e_rxn_surf[i], s_rxn_surf)
    if j >= 0:
        new_rxn_list_surf.append(s_rxn_surf[j])
    else:
        display(e_rxn_surf[i])

# -

new_rxn_list_gas

# +
# save the result
total_sp_list = new_sp_list_gas + new_sp_list_surf
rmgpy.chemkin.save_chemkin_file('my_mech-gas.inp', new_sp_list_gas, new_rxn_list_gas)
rmgpy.chemkin.save_species_dictionary('my_species.txt', total_sp_list)

rmgpy.chemkin.save_chemkin_surface_file('my_mech-surface.inp', new_sp_list_surf, new_rxn_list_surf)


# +
# Splice emily's thermo into my RMG model

new_sp_list_surf = []  # use emily's


# pick the surface thermo from Emily's if it exists
for i in range(len(s_sp_surf)):
    new_sp_list_surf.append(s_sp_surf[i])
    j = get_i_thing(s_sp_surf[i], e_sp_surf)
    if j >= 0:
        new_sp_list_surf[i].thermo = e_sp_surf[j].thermo
        print(f'replacing thermo for {new_sp_list_surf[i].smiles}')

rmgpy.chemkin.save_chemkin_file('emily_thermo-gas.inp', s_sp_gas, s_rxn_gas)
rmgpy.chemkin.save_chemkin_surface_file('emily_thermo-surface.inp', new_sp_list_surf, s_rxn_surf)

# -

# Check the workflow. Just choose emily's everything
rmgpy.chemkin.save_chemkin_file('emily_copy-gas.inp', e_sp_gas, e_rxn_gas)
rmgpy.chemkin.save_chemkin_surface_file('emily_copy-surface.inp', e_sp_surf, e_rxn_surf)


# +
# Splice emily's thermo into my RMG model, but also only use the reactions that are in Emily's model

new_sp_list_surf = []  # use emily's

new_rxn_list_gas = []
new_rxn_list_surf = []

# pick the surface thermo from Emily's if it exists
for i in range(len(s_sp_surf)):
    new_sp_list_surf.append(s_sp_surf[i])
    j = get_i_thing(s_sp_surf[i], e_sp_surf)
    if j >= 0:
        new_sp_list_surf[i].thermo = e_sp_surf[j].thermo
        print(f'replacing thermo for {new_sp_list_surf[i].smiles}')

        
for i in range(len(e_rxn_gas)):
    j = get_i_thing(e_rxn_gas[i], s_rxn_gas)
    if j >= 0:
        new_rxn_list_gas.append(s_rxn_gas[j])
    else:
        display(e_rxn_gas[i])
        
for i in range(len(e_rxn_surf)):
    j = get_i_thing(e_rxn_surf[i], s_rxn_surf)
    if j >= 0:
        new_rxn_list_surf.append(s_rxn_surf[j])
    else:
        display(e_rxn_surf[i])

        
        
rmgpy.chemkin.save_chemkin_file('emily_thermo_intersect-gas.inp', s_sp_gas, new_rxn_list_gas)
rmgpy.chemkin.save_chemkin_surface_file('emily_thermo_intersect-surface.inp', new_sp_list_surf, new_rxn_list_surf)

# -



# +
# Use emily's kinetics and my thermo
new_sp_list_gas = copy.deepcopy(s_sp_gas)
new_sp_list_surf = copy.deepcopy(s_sp_surf)

new_rxn_list_gas = []
new_rxn_list_surf = []

# add the species from Emily's if it's not in the model already
for i in range(len(e_sp_gas)):
    if e_sp_gas[i].smiles not in [x.smiles for x in new_sp_list_gas]:
        new_sp_list_gas.append(e_sp_gas[i])
for i in range(len(e_sp_surf)):
    if e_sp_surf[i].smiles not in [x.smiles for x in new_sp_list_surf]:
        new_sp_list_surf.append(e_sp_surf[i])
        
new_sp_list_total = new_sp_list_gas + new_sp_list_surf

# SPLICE - use my species to fill out the kinetics
for i in range(len(e_rxn_gas)):
    my_rxn = copy.deepcopy(e_rxn_gas[i])
    products = []
    for p in my_rxn.products:
        j = get_i_thing(p, new_sp_list_total)
        assert j >= 0
        products.append(new_sp_list_total[j])
    my_rxn.products = products
    
    reactants = []
    for r in my_rxn.reactants:
        j = get_i_thing(r, new_sp_list_total)
        assert j >= 0
        reactants.append(new_sp_list_total[j])
    my_rxn.reactants = reactants
    
    new_rxn_list_gas.append(my_rxn)

for i in range(len(e_rxn_surf)):
    my_rxn = copy.deepcopy(e_rxn_surf[i])
    products = []
    for p in my_rxn.products:
        j = get_i_thing(p, new_sp_list_total)
        assert j >= 0
        products.append(new_sp_list_total[j])
    my_rxn.products = products
    
    reactants = []
    for r in my_rxn.reactants:
        j = get_i_thing(r, new_sp_list_total)
        assert j >= 0
        reactants.append(new_sp_list_total[j])
    my_rxn.reactants = reactants
    
    new_rxn_list_surf.append(my_rxn)


        
rmgpy.chemkin.save_chemkin_file('emily_kinetics-gas.inp', new_sp_list_gas, new_rxn_list_gas)
rmgpy.chemkin.save_chemkin_surface_file('emily_kinetics-surface.inp', new_sp_list_surf, new_rxn_list_surf)


# +
[x.smiles for x in new_sp_list_gas]

for i in range(len(new_sp_list_gas)):
    print(i, new_sp_list_gas[-1].is_isomorphic(new_sp_list_gas[i]))
    print(i, new_sp_list_gas[-1].smiles == new_sp_list_gas[i].smiles)
# -







# +
new_sp_list_gas = []
new_rxn_list_gas = []
new_sp_list_surf = []
new_rxn_list_surf = []

# pick the gases
for i in range(len(e_sp_gas)):
    j = get_i_thing(e_sp_gas[i], s_sp_gas)
    if j >= 0:
        new_sp_list_gas.append(s_sp_gas[j])
    else:
        display(e_sp_gas[i])
        new_sp_list_gas.append(e_sp_gas[i])
        
for i in range(len(e_rxn_gas)):
    j = get_i_thing(e_rxn_gas[i], s_rxn_gas)
    if j >= 0:
        new_rxn_list_gas.append(s_rxn_gas[j])
    else:
        display(e_rxn_gas[i])
        new_rxn_list_gas.append(e_rxn_gas[i])

# pick the surfaces
for i in range(len(e_sp_surf)):
    j = get_i_thing(e_sp_surf[i], s_sp_surf)
    if j >= 0:
        new_sp_list_surf.append(s_sp_surf[j])
    else:
        display(e_sp_surf[i])
        new_sp_list_surf.append(e_sp_surf[i])
        
for i in range(len(e_rxn_surf)):
    j = get_i_thing(e_rxn_surf[i], s_rxn_surf)
    if j >= 0:
        new_rxn_list_surf.append(s_rxn_surf[j])
    else:
        display(e_rxn_surf[i])
        new_rxn_list_surf.append(e_rxn_surf[i])

        
# save the result
total_sp_list = new_sp_list_gas + new_sp_list_surf
rmgpy.chemkin.save_chemkin_file('my_mech2-gas.inp', new_sp_list_gas, new_rxn_list_gas)
# rmgpy.chemkin.save_species_dictionary('my_species.txt', total_sp_list)

rmgpy.chemkin.save_chemkin_surface_file('my_mech2-surface.inp', new_sp_list_surf, new_rxn_list_surf)

# -



new_sp_list_surf

new_sp_list_gas

len(new_sp_list)

len(new_rxn_list)

len(e_sp)

len(e_rxn)

for r in s_rxn:
    display(r)


