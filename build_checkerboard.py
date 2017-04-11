import openmc
import openmc.mgxs as mgxs
import numpy as np
import copy

def duplicate(universe):
    thing2 = copy.deepcopy(universe)
    thing2.id = None
    old_cells = copy.deepcopy(thing2.cells)
    thing2._cells.clear()
    for id in old_cells:
        new_cell = copy.deepcopy(old_cells[id])
        new_cell.id = None
        thing2.add_cell(new_cell)
    return thing2

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
fuel = openmc.Material(material_id=1, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide('U235', 1.)

moderator = openmc.Material(material_id=2, name='moderator')
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
left = openmc.XPlane(surface_id=1, x0=-2, name='left')
right = openmc.XPlane(surface_id=2, x0=2, name='right')
bottom = openmc.YPlane(surface_id=3, y0=-2, name='bottom')
top = openmc.YPlane(surface_id=4, y0=2, name='top')
fuel1 = openmc.ZCylinder(surface_id=5, x0=0, y0=0, R=0.4)
fuel2 = openmc.ZCylinder(surface_id=6, x0=0, y0=0, R=0.3)
fuel3 = openmc.ZCylinder(surface_id=7, x0=0, y0=0, R=0.2)

left.boundary_type = 'vacuum'
right.boundary_type = 'vacuum'
top.boundary_type = 'vacuum'
bottom.boundary_type = 'vacuum'

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='Cell 1')
cell2 = openmc.Cell(cell_id=101, name='cell 2')
cell3 = openmc.Cell(cell_id=102, name='cell 3')
cell4 = openmc.Cell(cell_id=201, name='cell 4')
cell5 = openmc.Cell(cell_id=202, name='cell 5')
cell6 = openmc.Cell(cell_id=301, name='cell 6')
cell7 = openmc.Cell(cell_id=302, name='cell 7')

# Use surface half-spaces to define regions
cell1.region = +left & -right & +bottom & -top
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
univ1 = openmc.Universe(universe_id=1)
univ2 = openmc.Universe(universe_id=2)
univ3 = openmc.Universe(universe_id=3)
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
univ1.add_cells([cell2, cell3])
univ2.add_cells([cell4, cell5])
univ3.add_cells([cell6, cell7])
root.add_cell(cell1)

# Instantiate a Lattice
lattice = openmc.RectLattice(lattice_id=5)
lattice.lower_left = [-2., -2.]
lattice.pitch = [1., 1.]
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
bounds = [-1, -1, -1, 1, 1, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.trigger_active = True
settings_file.trigger_max_batches = 100
settings_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC plots.xml file
###############################################################################

plot = openmc.Plot(plot_id=1)
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

# Instantiate a tally mesh
mesh = openmc.Mesh(mesh_id=1)
mesh.type = 'regular'
mesh.dimension = [4, 4]
mesh.lower_left = [-2, -2]
mesh.width = [1, 1]

# Instantiate tally Filter
mesh_filter = openmc.MeshFilter(mesh)


two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])

mesh_lib = mgxs.Library(geometry)
mesh_lib.energy_groups = two_groups
mesh_lib.mgxs_types = ['fission', 'nu-fission', 'transport', 'chi', 'scatter matrix']
mesh_lib.by_nuclide = False
mesh_lib.domain_type = "cell"
mesh_lib.domains = geometry.get_all_material_cells().values()
mesh_lib.build_library()


# Instantiate a Tallies collection and export to XML
tallies_file = openmc.Tallies()
mesh_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()


openmc.run()

# Export to MOC
import openmoc
import openmoc.materialize
import openmc.openmoc_compatible
sp = openmc.StatePoint('statepoint.020.h5')
mesh_lib.load_from_statepoint(sp)

moc_geom = openmc.openmoc_compatible.get_openmoc_geometry(sp.summary.geometry)
materials = openmoc.materialize.load_openmc_mgxs_lib(mesh_lib, moc_geom)
# Generate tracks for OpenMOC
# note: increase num_azim and decrease azim_spacing for actual results (as for TREAT)
track_generator = openmoc.TrackGenerator(moc_geom, num_azim=32, azim_spacing=0.1)
track_generator.generateTracks()

# Run OpenMOC
solver = openmoc.CPUSolver(track_generator)
solver.computeEigenvalue()

