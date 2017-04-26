import openmc
import openmc.mgxs as mgxs
import numpy as np
import copy


def duplicate(universe):
	cell_copy = copy.deepcopy(universe)
	cell_copy.id = None
	old_cells = copy.deepcopy(cell_copy.cells)
	cell_copy._cells.clear()
	for cell_id in old_cells:
		new_cell = copy.deepcopy(old_cells[cell_id])
		new_cell.id = None
		cell_copy.add_cell(new_cell)
	return cell_copy


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 20
inactive = 10
particles = 10000

###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide('U235', 1.)

moderator = openmc.Material(name='moderator')
moderator.set_density('g/cc', 1.0)
moderator.add_element('H', 2.)
moderator.add_element('O', 1.)
moderator.add_s_alpha_beta('c_H_in_H2O')

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([moderator, fuel])
materials_file.export_to_xml()

###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate Surfaces
min_x = openmc.XPlane(x0=-2, name='min-x')
max_x = openmc.XPlane(x0=+2, name='max-x')
min_y = openmc.YPlane(y0=-2, name='min-y')
max_y = openmc.YPlane(y0=+2, name='max-y')
min_z = openmc.ZPlane(z0=-2, name='min-z')
max_z = openmc.ZPlane(z0=+2, name='max-z')
fuel1 = openmc.ZCylinder(x0=0, y0=0, R=0.4)
fuel2 = openmc.ZCylinder(x0=0, y0=0, R=0.3)
fuel3 = openmc.ZCylinder(x0=0, y0=0, R=0.2)

min_x.boundary_type = 'reflective'
max_x.boundary_type = 'reflective'
min_y.boundary_type = 'reflective'
max_y.boundary_type = 'reflective'
min_z.boundary_type = 'reflective'
max_z.boundary_type = 'reflective'

# Instantiate Cells
cell1 = openmc.Cell(name='Cell 1')
cell2 = openmc.Cell(name='cell 2')
cell3 = openmc.Cell(name='cell 3')
cell4 = openmc.Cell(name='cell 4')
cell5 = openmc.Cell(name='cell 5')
cell6 = openmc.Cell(name='cell 6')
cell7 = openmc.Cell(name='cell 7')

# Use surface half-spaces to define regions
cell1.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
cell2.region = -fuel1
cell3.region = +fuel1
cell4.region = -fuel2
cell5.region = +fuel2
cell6.region = -fuel3
cell7.region = +fuel3

# Register Materials with Cells
cell2.fill = fuel
cell3.fill = moderator
cell4.fill = fuel
cell5.fill = moderator
cell6.fill = fuel
cell7.fill = moderator

# Instantiate Universe
univ1 = openmc.Universe(name='small pin')
univ2 = openmc.Universe(name='medium pin')
univ3 = openmc.Universe(name='big pin')
root = openmc.Universe(name='root universe')

# Register Cells with Universe
univ1.add_cells([cell2, cell3])
univ2.add_cells([cell4, cell5])
univ3.add_cells([cell6, cell7])
root.add_cell(cell1)

# Instantiate a Lattice
lattice = openmc.RectLattice()
lattice.lower_left = np.array([-2., -2.])
lattice.pitch = np.array([1., 1.])
lattice.universes = [[univ1, univ2, univ1, univ2],
                     [univ2, univ3, univ2, univ3],
                     [univ1, univ2, univ1, univ2],
                     [univ2, univ3, univ2, univ3]]

for i in lattice.indices:
	lattice.universes[i] = duplicate(lattice.universes[i])

# Fill Cell with the Lattice
cell1.fill = lattice

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()

###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-2, -2, -2, 2, 2, 2]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space = uniform_dist)

settings_file.run_mode = "eigenvalue"
settings_file.export_to_xml()

###############################################################################
#                   Exporting to OpenMC plots.xml file
###############################################################################

plot = openmc.Plot()
plot.origin = [0, 0, 0]
plot.width = [4, 4]
plot.pixels = [400, 400]
plot.color_by = 'cell'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()

###############################################################################
#                   Exporting to OpenMC tallies.xml file
###############################################################################


############################################
# Begin treat-like stuff
############################################

# Instantiate a tally Mesh
mesh = openmc.Mesh()
mesh.lower_left = lattice.lower_left
mesh.upper_right = -lattice.lower_left
mesh.type = 'regular'
mesh.dimension = lattice.shape

# Two energy groups
two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])

# Mesh MGXS library
mesh_lib = mgxs.Library(geometry)
mesh_lib.energy_groups = two_groups
mesh_lib.mgxs_types = ['fission', 'nu-fission', 'total', 'chi', 'consistent nu-scatter matrix']
mesh_lib.by_nuclide = False
mesh_lib.correction = None
mesh_lib.domain_type = "mesh"
mesh_lib.domains = [mesh]
mesh_lib.build_library()
mesh_lib.dump_to_file("mesh_lib")

# TODO: select all the nuclides from `mats`


#######################
'''
# Instantiate a tally mesh
old_mesh = openmc.Mesh(mesh_id = 1)
old_mesh.type = 'regular'
old_mesh.dimension = [4, 4]
old_mesh.lower_left = [-2, -2]
old_mesh.width = [1, 1]

two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])

mesh_lib = mgxs.Library(geometry)
mesh_lib.energy_groups = two_groups
mesh_lib.mgxs_types = ['fission', 'nu-fission', 'transport', 'chi', 'scatter matrix']
mesh_lib.by_nuclide = False
mesh_lib.domain_type = "cell"
mesh_lib.domains = geometry.get_all_material_cells().values()
mesh_lib.build_library()
'''

# Instantiate a Tallies collection and export to XML
tallies_file = openmc.Tallies()
mesh_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()


if __name__ == "__main__":
	openmc.run()

