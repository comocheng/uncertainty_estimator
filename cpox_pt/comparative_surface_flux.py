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

# # Comparative Flux Diagram for Surface Mechanisms

# +
import os
import sys
import copy
import numpy as np
import pydot
import rmgpy.tools.fluxdiagram
import rmgpy.chemkin

import matplotlib
import matplotlib.pyplot as plt
# %matplotlib inline

# -

# # Settings

# +
# Options controlling the individual flux diagram renderings:
program = 'dot'  # The program to use to lay out the nodes and edges
max_node_count = 35  # The maximum number of nodes to show in the diagram
max_edge_count = 35  # The maximum number of edges to show in the diagram
concentration_tol = 1e-6  # The lowest fractional concentration to show (values below this will appear as zero)
species_rate_tol = 1e-6  # The lowest fractional species rate to show (values below this will appear as zero)
max_node_pen_width = 7.0  # The thickness of the border around a node at maximum concentration
max_edge_pen_width = 9.0  # The thickness of the edge at maximum species rate
radius = 1  # The graph radius to plot around a central species
central_reaction_count = None  # The maximum number of reactions to draw from each central species (None draws all)
# If radius > 1, then this is the number of reactions from every species

# Options controlling the ODE simulations:
initial_time = 1e-12  # The time at which to initiate the simulation, in seconds
time_step = 10 ** 0.1  # The multiplicative factor to use between consecutive time points
abs_tol = 1e-16  # The absolute tolerance to use in the ODE simluations
rel_tol = 1e-8  # The relative tolerance to use in the ODE simulations

# Options controlling the generated movie:
video_fps = 6  # The number of frames per second in the generated movie
initial_padding = 5  # The number of seconds to display the initial fluxes at the start of the video
final_padding = 5  # The number of seconds to display the final fluxes at the end of the video
# -

# # Load Mechanisms

# +
# species_path = None
java = False            # always False
settings = None
chemkin_output = ''     # this will be generated automatically
central_species_list = None
superimpose = False     # this will always be false, delete it
save_states = False
read_states = False     # fine to keep this always false and delete relevant code below
diffusion_limited = True
check_duplicates = True

rmg_input_file = '/home/moon/rmg/PR/surface_flux_diagram/input.py'  # for conditions T = 1000

# diagram_base_name = 'surface_comparison'
# os.makedirs(diagram_base_name, exist_ok=True)
# mech_1_inp = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin/chem_annotated-gas.inp'
# mech_1_surf = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin/chem_annotated-surface.inp'
# mech_1_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/mincat/chemkin/species_dictionary.txt'
# mech_1_label = 'cpox_pt_20241112_mincat'

# mech_2_inp = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/chem_annotated-gas.inp'
# mech_2_surf = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/chem-surface.inp'
# mech_2_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/species_dictionary.txt'
# mech_2_label = 'cpox_rh_emily'

# diagram_base_name = 'surface_comparison_katrinPt111'
# os.makedirs(diagram_base_name, exist_ok=True)
# mech_1_inp = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/debug_OX_scaling/chemkin/chem_annotated-gas.inp'
# mech_1_surf = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/debug_OX_scaling/chemkin/chem_annotated-surface.inp'
# mech_1_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/debug_OX_scaling/chemkin/species_dictionary.txt'
# mech_1_label = 'cpox_pt_20241112_mincat'

# mech_2_inp = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/chem_annotated-gas.inp'
# mech_2_surf = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/chem-surface.inp'
# mech_2_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/species_dictionary.txt'
# mech_2_label = 'cpox_rh_emily'

diagram_base_name = 'surface_comparison_stick'
os.makedirs(diagram_base_name, exist_ok=True)
mech_1_inp = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/handpicked/chem_annotated-gas.inp'
mech_1_surf = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/handpicked/chem_annotated-surface.inp'
mech_1_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_20241112/handpicked/species_dictionary.txt'
mech_1_label = 'cpox_pt_stick_0p001'

mech_2_inp = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/chem_annotated-gas.inp'
mech_2_surf = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/chem-surface.inp'
mech_2_dict = '/home/moon/uncertainty_estimator/cpox_pt/cpox_rh_emily/species_dictionary.txt'
mech_2_label = 'cpox_rh_emily'


species_path = '/home/moon/rmg/PR/surface_flux_diagram/species'

generate_images = False
print('Loading RMG job 1...')
rmg_job1 = rmgpy.tools.fluxdiagram.load_rmg_job(
    rmg_input_file,
    mech_1_inp,
    mech_1_dict,
    generate_images=generate_images,
    check_duplicates=check_duplicates,
    surface_path=mech_1_surf
)

print('Loading RMG job 2...')
rmg_job2 = rmgpy.tools.fluxdiagram.load_rmg_job(
    rmg_input_file,
    mech_2_inp,
    mech_2_dict,
    generate_images=generate_images,
    check_duplicates=check_duplicates,
    surface_path=mech_2_surf
)

# -

# # Simulation

# +
print('Conducting simulation of reaction 1')
times1, concentrations1, reaction_rates1 = rmgpy.tools.fluxdiagram.simulate(
    rmg_job1.reaction_model,
    rmg_job1.reaction_systems[0],
    settings
)

print('Conducting simulation of reaction 2')
times2, concentrations2, reaction_rates2 = rmgpy.tools.fluxdiagram.simulate(
    rmg_job2.reaction_model,
    rmg_job2.reaction_systems[0],
    settings
)

# -

# # Compute/assemble species concentrations and fluxes into convenient form

# +
# Get the RMG species and reactions objects
species_list1 = rmg_job1.reaction_model.core.species[:]
reaction_list1 = rmg_job1.reaction_model.core.reactions[:]
num_species1 = len(species_list1)

species_list2 = rmg_job2.reaction_model.core.species[:]
reaction_list2 = rmg_job2.reaction_model.core.reactions[:]
num_species2 = len(species_list2)

# +
# Compute the rates between each pair of species (build up big matrices)
species_rates1 = np.zeros((len(times1), num_species1, num_species1), float)
for index1, reaction1 in enumerate(reaction_list1):
    rate1 = reaction_rates1[:, index1]
    if not reaction1.pairs: reaction1.generate_pairs()
    for reactant1, product1 in reaction1.pairs:
        reactant_index1 = species_list1.index(reactant1)
        product_index1 = species_list1.index(product1)
        species_rates1[:, reactant_index1, product_index1] += rate1
        species_rates1[:, product_index1, reactant_index1] -= rate1
        
species_rates2 = np.zeros((len(times2), num_species2, num_species2), float)
for index2, reaction2 in enumerate(reaction_list2):
    rate2 = reaction_rates2[:, index2]
    if not reaction2.pairs: reaction2.generate_pairs()
    for reactant2, product2 in reaction2.pairs:
        reactant_index2 = species_list2.index(reactant2)
        product_index2 = species_list2.index(product2)
        species_rates2[:, reactant_index2, product_index2] += rate2
        species_rates2[:, product_index2, reactant_index2] -= rate2

# +
# Determine the maximum concentration for each species and the maximum overall concentration
max_concentrations1 = np.max(np.abs(concentrations1), axis=0)
max_concentration1 = np.max(max_concentrations1)

# Determine the maximum reaction rates
max_reaction_rates1 = np.max(np.abs(reaction_rates1), axis=0)

# Determine the maximum rate for each species-species pair and the maximum overall species-species rate
max_species_rates1 = np.max(np.abs(species_rates1), axis=0)
max_species_rate1 = np.max(max_species_rates1)
species_index1 = max_species_rates1.reshape((num_species1 * num_species1)).argsort()


max_concentrations2 = np.max(np.abs(concentrations2), axis=0)
max_concentration2 = np.max(max_concentrations2)

# Determine the maximum reaction rates
max_reaction_rates2 = np.max(np.abs(reaction_rates2), axis=0)

# Determine the maximum rate for each species-species pair and the maximum overall species-species rate
max_species_rates2 = np.max(np.abs(species_rates2), axis=0)
max_species_rate2 = np.max(max_species_rates2)
species_index2 = max_species_rates2.reshape((num_species2 * num_species2)).argsort()


max_species_rate_total = max(max_species_rate1, max_species_rate2)
max_concentration_total = max(max_concentration1, max_concentration2)
# -

# # Determine nodes and edges to include for each model

# +
nodes1 = []
edges1 = []
for i in range(num_species1 * num_species1):
    product_index1, reactant_index1 = divmod(species_index1[-i - 1], num_species1)
    if reactant_index1 > product_index1:
        # Both reactant -> product and product -> reactant are in this list, so only keep one of them
        continue
    if max_species_rates1[reactant_index1, product_index1] == 0:
        break
    if reactant_index1 not in nodes1 and len(nodes1) < max_node_count: nodes1.append(reactant_index1)
    if product_index1 not in nodes1 and len(nodes1) < max_node_count: nodes1.append(product_index1)
    if [reactant_index1, product_index1] not in edges1 and [product_index1, reactant_index1] not in edges1:
        edges1.append([reactant_index1, product_index1])
    if len(nodes1) > max_node_count:
        break
    if len(edges1) >= max_edge_count:
        break
        
        
nodes2 = []
edges2 = []
for i in range(num_species2 * num_species2):
    product_index2, reactant_index2 = divmod(species_index2[-i - 1], num_species2)
    if reactant_index2 > product_index2:
        # Both reactant -> product and product -> reactant are in this list, so only keep one of them
        continue
    if max_species_rates2[reactant_index2, product_index2] == 0:
        break
    if reactant_index2 not in nodes2 and len(nodes2) < max_node_count: nodes2.append(reactant_index2)
    if product_index2 not in nodes2 and len(nodes2) < max_node_count: nodes2.append(product_index2)
    if [reactant_index2, product_index2] not in edges2 and [product_index2, reactant_index2] not in edges2:
        edges2.append([reactant_index2, product_index2])
    if len(nodes2) > max_node_count:
        break
    if len(edges2) >= max_edge_count:
        break
         
# -

# # Create mapping of species indices between the two models

# +
species2_to_1 = {}
species2_to_1 = {key: value for key, value in zip([x for x in range(len(species_list1))], [-1] * len(species_list1))}
for i in range(len(species_list2)):
    for j in range(len(species_list1)):
        if species_list2[i].is_isomorphic(species_list1[j]):
            species2_to_1[i] = j
            break
    else:
        species2_to_1[i] = -1

species1_to_2 = {}
species1_to_2 = {key: value for key, value in zip([x for x in range(len(species_list1))], [-1] * len(species_list1))}
for i in range(len(species_list1)):
    for j in range(len(species_list2)):
        if species_list1[i].is_isomorphic(species_list2[j]):
            species1_to_2[i] = j
            break
    else:
        species1_to_2[i] = -1
        
# -

# Function to grab and generate the image for a given species
def get_image_path(species):
    species_index = str(species) + '.png'
    image_path = ''
    if not species_path or not os.path.exists(species_path):  # species_path is defined while loading the mechanism
        raise OSError
    for root, dirs, files in os.walk(species_path):
        for f in files:
            if f == species_index:
                image_path = os.path.join(root, f)
                break
    if not image_path:
        image_path = os.path.join(species_path, species_index)
        species.molecule[0].draw(image_path)
    return image_path


[x for x in np.logspace(-15, 0, 20)]

# # Build the Graph

# +
# Create the graph
mech1_color = matplotlib.colors.to_hex((1.0, 0.92, 0.0))
mech2_color = matplotlib.colors.to_hex((0.18627451, 0.48823529, 0.94117647))

# add alpha specified by adding hex digits before RGB values
alpha = 0.9
alpha_hex = hex(int(alpha * 255))[2:]
mech1_color = f'#{mech1_color[1:]}{alpha_hex}'
mech2_color = f'#{mech2_color[1:]}{alpha_hex}'

# Grab the fluxes from the time closest (without going over) to each snapshot time
t1_snapshot = 0
t2_snapshot = 0


snapshot_times = np.logspace(-13, -9, 40)
for snapshot_time in snapshot_times:
    t1_snapshot = snapshot_time
    t2_snapshot = snapshot_time

    t1 = np.abs(times1 - t1_snapshot).argmin()
    if times1[t1] > t1_snapshot and t1 != 0:
        t1 -= 1
    t2 = np.abs(times2 - t2_snapshot).argmin()
    if times2[t2] > t2_snapshot and t2 != 0:
        t2 -= 1

    assert -1 not in nodes1
    assert -1 not in nodes2


    slope = -max_node_pen_width / np.log10(concentration_tol)
    graph = pydot.Dot('flux_diagram', graph_type='digraph', overlap="false")
    graph.set_fontname('sans')
    graph.set_fontsize('10')


    # ----------------------------ADD NODES ------------------------------#
    # For Mechanism 1
    for index1 in nodes1:
        nodewidths = np.zeros(2)  # keep track of species concentrations/nodewidths for both mechanisms
        species1 = species_list1[index1]
        node1 = pydot.Node(name=str(species1))
        concentration1 = concentrations1[t1, index1] / max_concentration_total
        if concentration1 < concentration_tol:
            penwidth = 0.0
        else:
            penwidth = round(slope * np.log10(concentration1) + max_node_pen_width, 3)
            nodewidths[0] = penwidth
        node1.set_penwidth(penwidth)
        node1.set_fillcolor('white')
        node1.set_color(mech1_color)
        image_path1 = get_image_path(species1)
        if os.path.exists(image_path1):
            node1.set_image(image_path1)
            node1.set_label("")

        index2 = species1_to_2[index1]
        if index2 >= 0:
            concentration2 = concentrations2[t2, index2] / max_concentration_total
            if concentration2 < concentration_tol:
                penwidth = 0.0
            else:
                penwidth = round(slope * np.log10(concentration2) + max_node_pen_width, 3)
                nodewidths[1] = penwidth
                if node1.get_penwidth() > 0:
                    node1.set_color('black')
                else:
                    node1.set_color(mech2_color)
                    node1.set_penwidth(penwidth)

        if node1.get_color() == 'black':
            node1.set_penwidth(np.average(nodewidths))
        graph.add_node(node1)


    # For Mechanism 2
    for index2 in nodes2:
        if species2_to_1[index2] in nodes1:
            continue  # already took care of it above

        nodewidths = np.zeros(2)
        species2 = species_list2[index2]
        node2 = pydot.Node(name=str(species2))
        concentration2 = concentrations2[t2, index2] / max_concentration_total
        if concentration2 < concentration_tol:
            penwidth = 0.0
        else:
            penwidth = round(slope * np.log10(concentration2) + max_node_pen_width, 3)
            nodewidths[1] = penwidth
        node2.set_fillcolor('white')
        node2.set_color(mech2_color)
        node2.set_penwidth(penwidth)
        # Try to use an image instead of the label
        image_path2 = get_image_path(species2)
        if os.path.exists(image_path2):
            node2.set_image(image_path2)
            node2.set_label("")

        if node2.get_color() == 'black':
            node2.set_penwidth(np.average(nodewidths))

        graph.add_node(node2)


    # ------------------------------- EDGES ------------------------------#
    # Add an edge for each species-species rate
    slope = -max_edge_pen_width / np.log10(species_rate_tol)

    # Go through edges in Mechanism 1
    for reactant_index1, product_index1 in edges1:
        if reactant_index1 in nodes1 and product_index1 in nodes1:
            reactant1 = species_list1[reactant_index1]
            product1 = species_list1[product_index1]

            edge1 = pydot.Edge(str(reactant1), str(product1), color=mech1_color)
            species_rate1 = species_rates1[t1, reactant_index1, product_index1] / max_species_rate_total
            if species_rate1 < 0:
                edge1.set_dir("back")
                species_rate1 = -species_rate1
            else:
                edge1.set_dir("forward")
            # Set the edge pen width
            if species_rate1 < species_rate_tol:
                penwidth = 0.0
                edge1.set_dir("none")
            else:
                penwidth = round(slope * np.log10(species_rate1) + max_edge_pen_width, 3)
            edge1.set_penwidth(penwidth)

            graph.add_edge(edge1)

            # add mech 2
            if species1_to_2[reactant_index1] >= 0 and species1_to_2[product_index1] >= 0:
                reactant_index2 = species1_to_2[reactant_index1]
                product_index2 = species1_to_2[product_index1]

                edge2 = pydot.Edge(str(reactant1), str(product1), color=mech2_color)
                species_rate2 = species_rates2[t2, reactant_index2, product_index2] / max_species_rate_total
                if species_rate2 < 0:
                    edge2.set_dir("back")
                    species_rate2 = -species_rate2
                else:
                    edge2.set_dir("forward")
                # Set the edge pen width
                if species_rate2 < species_rate_tol:
                    penwidth = 0.0
                    edge2.set_dir("none")
                else:
                    penwidth = round(slope * np.log10(species_rate2) + max_edge_pen_width, 3)
                edge2.set_penwidth(penwidth)
                graph.add_edge(edge2)

    # Go through edges in Mechanism 2
    for reactant_index2, product_index2 in edges2:
        # skip if this was already done in edges 1
        if [species2_to_1[reactant_index2], species2_to_1[product_index2]] in edges1 or \
            [species2_to_1[product_index2], species2_to_1[reactant_index2]] in edges1:
            continue

        if reactant_index2 in nodes2 and product_index2 in nodes2:
            # mech 2 says include this edge for all mechs
            if species2_to_1[reactant_index2] in nodes1:
                reactant2 = species_list1[species2_to_1[reactant_index2]]
            else:
                reactant2 = species_list2[reactant_index2]
            if species2_to_1[product_index2] in nodes1:
                product2 = species_list1[species2_to_1[product_index2]]
            else:
                product2 = species_list2[product_index2]
            edge2 = pydot.Edge(str(reactant2), str(product2), color=mech2_color)

            species_rate2 = species_rates2[t2, reactant_index2, product_index2] / max_species_rate_total
            if species_rate2 < 0:
                edge2.set_dir("back")
                species_rate2 = -species_rate2
            else:
                edge2.set_dir("forward")
            # Set the edge pen width
            if species_rate2 < species_rate_tol:
                penwidth = 0.0
                edge2.set_dir("none")
            else:
                penwidth = round(slope * np.log10(species_rate2) + max_edge_pen_width, 3)

            edge2.set_penwidth(penwidth)
            graph.add_edge(edge2)

            # add mech 1
            if species2_to_1[reactant_index2] >= 0 and species2_to_1[product_index2] >= 0:
                reactant_index1 = species2_to_1[reactant_index2]
                product_index1 = species2_to_1[product_index2]

                edge1 = pydot.Edge(str(reactant2), str(product2), color=mech1_color)
                species_rate1 = species_rates1[t1, reactant_index1, product_index1] / max_species_rate_total
                if species_rate1 < 0:
                    edge1.set_dir("back")
                    species_rate1 = -species_rate1
                else:
                    edge1.set_dir("forward")
                # Set the edge pen width
                if species_rate1 < species_rate_tol:
                    penwidth = 0.0
                    edge1.set_dir("none")
                else:
                    penwidth = round(slope * np.log10(species_rate1) + max_edge_pen_width, 3)
                edge1.set_penwidth(penwidth)
                graph.add_edge(edge1)


    # # General purpose graph settings
    # graph.set_nodesep(0.11)
    # graph.set_ranksep(0.35)
    # graph.set_rankdir('LR')

    # Add Legend
    graph.add_node(pydot.Node(mech_1_label + f'\nt={times1[t1]:.4e}', label=mech_1_label + f'\nt={times1[t1]:.4e}', color=mech1_color, shape='box', penwidth=max_node_pen_width))
    graph.add_node(pydot.Node(mech_2_label + f'\nt={times2[t2]:.4e}', label=mech_2_label + f'\nt={times2[t2]:.4e}', color=mech2_color, shape='box', penwidth=max_node_pen_width))


    # write in multiple formats
    graph.write_dot(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1}.dot')) # Yes this is supposed to be an index, not an actual time
    graph.write_png(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1}.png'))
    # graph.write_pdf(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1}.pdf'))
    # graph.write_pdf(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1}.svg'))


# -

def plot_thermos(thermos, labels=None):
    if type(thermos) != list:
        thermos = [thermos]
    if labels is None:
        labels = ['' for t in thermos]
    markers = ['solid', 'dashed']
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(12, 3)
    fig.tight_layout()
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('H (kJ / mol)')
    ax[0].set_title('Enthalpy vs. Temperature')
    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('S (kJ / mol K)')
    ax[1].set_title('Entropy vs. Temperature')
    ax[2].set_xlabel('Temperature (K)')
    ax[2].set_ylabel('Cp (kJ / mol K)')
    ax[2].set_title('Heat Capacity vs. Temperature')
    T = np.linspace(300, 3000, 1001)
    for z, thermo in enumerate(thermos):
        H = np.zeros(len(T))
        S = np.zeros(len(T))
        Cp = np.zeros(len(T))
        for i in range(0, len(T)):
            H[i] = thermo.get_enthalpy(T[i]) / 1000.0
            S[i] = thermo.get_entropy(T[i]) / 1000.0
            Cp[i] = thermo.get_heat_capacity(T[i]) / 1000.0
        ax[0].plot(T, H, linestyle=markers[z])
        ax[1].plot(T, S, linestyle=markers[z])
        ax[2].plot(T, Cp, linestyle=markers[z])
    ax[0].legend(labels)
    ax[1].legend(labels)
    ax[2].legend(labels)
    plt.subplots_adjust(wspace=0.25)
    plt.show()


def get_i_thing(thing, thing_list):
    for j in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[j]):
            return j
    return -1


display(species_list1[27])

print(species_list1[27].to_adjacency_list())

# # Debug O2 Adsorption

rmgpy.species.Species().from_adjacency_list(
    """
    OX
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
    """
)

# +
# compare the O2 adsorption rates

O2 = rmgpy.species.Species(smiles='[O][O]')
OX = rmgpy.species.Species().from_adjacency_list(
    """
    OX
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
    """
)
my_reactant = O2
my_product = OX

idx1 = get_i_thing(O2, species_list1)
idx2 = get_i_thing(O2, species_list2)


matches1 = []
for i in range(len(reaction_list1)):
    for reactant in reaction_list1[i].reactants:
        if my_reactant.is_isomorphic(reactant):
            for product in reaction_list1[i].products:
                if my_product.is_isomorphic(product):
                    matches1.append(i)
        if my_product.is_isomorphic(reactant):
            for product in reaction_list1[i].products:
                if my_reactant.is_isomorphic(product):
                    matches1.append(i)
matches2 = []
for i in range(len(reaction_list2)):
    for reactant in reaction_list2[i].reactants:
        if my_reactant.is_isomorphic(reactant):
            for product in reaction_list2[i].products:
                if my_product.is_isomorphic(product):
                    matches2.append(i)
        if my_product.is_isomorphic(reactant):
            for product in reaction_list2[i].products:
                if my_reactant.is_isomorphic(product):
                    matches2.append(i)                    


matches1 = list(set(matches1))
matches2 = list(set(matches2))
for i in matches1:
    display(i, reaction_list1[i])
for i in matches2:
    display(i, reaction_list2[i])

# -

get_i_thing(OX, species_list2)

# +
# compare thermo
reaction_list1[49].reactants[0].thermo


# -

reaction_list2[42].products[0].thermo

species_list1[22].to_adjacency_list()

# +
OX = rmgpy.species.Species().from_adjacency_list(
    """
    OX
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
    """
)


X = rmgpy.species.Species().from_adjacency_list(
    """
    X
1 X u0 p0 c0
    """
)
O2 = rmgpy.species.Species(smiles='[O][O]')


my_sp = OX
display(my_sp)
idx1 = get_i_thing(my_sp, species_list1)
idx2 = get_i_thing(my_sp, species_list2)


plot_thermos([species_list1[idx1], species_list2[idx2]], ['mine', 'Emily'])
# -





# +
# plot kinetics of mech1 49 reversed

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

idx1 = 49
my_rxn1 = copy.deepcopy(reaction_list1[idx1])
my_rxn1.reactants = reaction_list1[idx1].products
my_rxn1.products = reaction_list1[idx1].reactants
my_rxn1.kinetics = reaction_list1[idx1].generate_reverse_rate_coefficient()

idx2 = 42
my_rxn2 = copy.deepcopy(reaction_list2[idx2])

display(my_rxn2)

T = np.linspace(300, 3000, 1001)
P = 101325
k1 = np.zeros(len(T))
k2 = np.zeros(len(T))
for j in range(0, len(T)):
    k1[j] = my_rxn1.get_rate_coefficient(T[j], P)
    k2[j] = my_rxn2.get_rate_coefficient(T[j], P)
plt.plot(1000.0 / T, k2, label=f'Emily Rh', color=colors[0])
plt.plot(1000.0 / T, k1, label=f'My Mechanism (reversed)', color=colors[1])
plt.legend()
plt.yscale('log')


plt.show()
T = 1000
display(my_rxn1)
energies1 = [0, my_rxn1.kinetics.Ea.value_si, my_rxn1.get_enthalpy_of_reaction(T)]
energies2 = [0, my_rxn2.kinetics.Ea.value_si, my_rxn2.get_enthalpy_of_reaction(T)]

plt.plot(energies1, label='My mech (reversed)', color=colors[1])
plt.plot(energies2, label='Emily Rh', color=colors[0])
plt.legend()
# -



my_rxn1.kinetics.Ea.value_si



# +
# plot kinetics of mech1 49 reversed

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

idx1 = 49
my_rxn1 = copy.deepcopy(reaction_list1[idx1])
# my_rxn1 = copy.deepcopy(reaction_list1[idx1])
# my_rxn1.reactants = reaction_list1[idx1].products
# my_rxn1.kinetics = reaction_list1[idx1].generate_reverse_rate_coefficient()

display(my_rxn1)

idx2 = 42
# my_rxn2 = copy.deepcopy(reaction_list2[idx2])
my_rxn2 = copy.deepcopy(reaction_list2[idx2])
my_rxn2.reactants = reaction_list2[idx2].products
my_rxn2.products = reaction_list2[idx2].reactants
my_rxn2.kinetics = reaction_list2[idx2].generate_reverse_rate_coefficient()

T = np.linspace(300, 3000, 1001)
P = 101325
k1 = np.zeros(len(T))
k2 = np.zeros(len(T))
for j in range(0, len(T)):
    k1[j] = my_rxn1.get_rate_coefficient(T[j], P)
    k2[j] = my_rxn2.get_rate_coefficient(T[j], P)
plt.plot(1000.0 / T, k2, label=f'Emily Rh', color=colors[0])
plt.plot(1000.0 / T, k1, label=f'My Mechanism', color=colors[1])
plt.legend()
plt.yscale('log')


plt.show()
T = 1000
energies1 = [0, my_rxn1.kinetics.Ea.value_si, my_rxn1.get_enthalpy_of_reaction(T)]
energies2 = [0, my_rxn2.kinetics.Ea.value_si, my_rxn2.get_enthalpy_of_reaction(T)]
display(my_rxn1)
plt.plot(energies1, label='My mech')
plt.plot(energies2, label='Emily Rh (revesed)')
plt.legend()
# -



my_rxn2.products

idx2

reaction_list2[42]

# +
r2_rev = copy.deepcopy(reaction_list2[42])

r2_rev.products = reaction_list2[42].reactants
r2_rev.reactants = reaction_list2[42].products
r2_rev.kinetics = reaction_list2[42].generate_reverse_rate_coefficient()

display(r2_rev)
print(r2_rev.kinetics)
# -

# # Debug OH formation

# +
# compare the O2 adsorption rates

O2 = rmgpy.species.Species(smiles='[O][O]')
OX = rmgpy.species.Species().from_adjacency_list(
    """
    OX
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
    """
)


X = rmgpy.species.Species().from_adjacency_list(
    """
    X
1 X u0 p0 c0
    """
)

OHX = rmgpy.species.Species().from_adjacency_list(
"""
OHX
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
"""
)

my_reactant = X
my_product = OHX

idx1 = get_i_thing(O2, species_list1)
idx2 = get_i_thing(O2, species_list2)


matches1 = []
for i in range(len(reaction_list1)):
    for reactant in reaction_list1[i].reactants:
        if my_reactant.is_isomorphic(reactant):
            for product in reaction_list1[i].products:
                if my_product.is_isomorphic(product):
                    matches1.append(i)
        if my_product.is_isomorphic(reactant):
            for product in reaction_list1[i].products:
                if my_reactant.is_isomorphic(product):
                    matches1.append(i)
matches2 = []
for i in range(len(reaction_list2)):
    for reactant in reaction_list2[i].reactants:
        if my_reactant.is_isomorphic(reactant):
            for product in reaction_list2[i].products:
                if my_product.is_isomorphic(product):
                    matches2.append(i)
        if my_product.is_isomorphic(reactant):
            for product in reaction_list2[i].products:
                if my_reactant.is_isomorphic(product):
                    matches2.append(i)                    


matches1 = list(set(matches1))
matches2 = list(set(matches2))
for i in matches1:
    display(i, reaction_list1[i])
    print(reaction_rates1[0, i])
print('--------------------------------------------------')
for i in matches2:
    display(i, reaction_list2[i])
    print(reaction_rates2[0, i])

# -

print(species_list1[33].to_adjacency_list())

reaction_list2[50].kinetics

reaction_list1[57].kinetics

reaction_list1[57]

reaction_rates1.shape

len(reaction_rates1[0, :])

species_rates1[:, reactant_index1, product_index1]

species_rates1.shape

len(times1)

len(species_list1)


