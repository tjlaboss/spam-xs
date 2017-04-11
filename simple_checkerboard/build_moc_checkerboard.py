#############
# Export to MOC
import openmoc
import openmoc.materialize
import openmc.openmoc_compatible
import openmoc.plotter as plt

import openmc.mgxs as mgxs
import numpy as np
import pandas


sp = openmc.StatePoint('statepoint.020.h5')
'''
# 2-group approximation
two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])
mesh_lib = mgxs.Library(sp.summary.geometry)
mesh_lib.energy_groups = two_groups

# For purposes of this demonstration, let's just look at the capture
# and the transport cross sections
mesh_lib.mgxs_types = ['fission', 'nu-fission', 'transport', 'chi', 'scatter matrix']
mesh_lib.by_nuclide = False
# mesh_lib.nuclides = ["U235", "U238"]
# TODO: select all the nuclides from `mats`
mesh_lib.domain_type = "mesh"
mesh_lib.load_from_statepoint(sp)
'''

mesh_lib = mgxs.Library.load_from_file(filename = "mesh_lib")
mesh_lib.load_from_statepoint(sp)
mesh = mesh_lib.domains[0]
x_width, y_width = (mesh.upper_right - mesh.lower_left)/mesh.dimension

# Build a checkerboard geometry in OpenMOC
lattice = openmoc.Lattice(name='4x4 lattice')
lattice.setWidth(x_width, y_width)

nx, ny = mesh.dimension[:2]
universes = [[None for i in range(nx)] for j in range (ny)]

for i in range(nx):
	for j in range(ny):
		c = openmoc.Cell()
		m = openmoc.Material()
		c.setFill(m)
		u = openmoc.Universe()
		u.addCell(c)
		universes[j][i] = u
print(universes)
lattice.setUniverses([universes])

print(lattice.getMinX(), lattice.getMaxX())
print(lattice.getMinY(), lattice.getMaxY())

root_universe = openmoc.Universe(name = "root universe")
root_cell = openmoc.Cell(name = "root cell")
root_cell.setFill(lattice)
# Make some boundaries
min_x, min_y = mesh.lower_left[:2]
max_x, max_y = mesh.upper_right[:2]
xleft =  openmoc.XPlane(x = min_x)
xright = openmoc.XPlane(x = max_x)
yleft =  openmoc.YPlane(y = min_y)
yright = openmoc.YPlane(y = max_y)
surfs = [xleft, xright, yleft, yright]
for s in surfs:
	s.setBoundaryType(openmoc.VACUUM)
root_cell.addSurface(+1, xleft)
root_cell.addSurface(-1, xright)
root_cell.addSurface(+1, yleft)
root_cell.addSurface(-1, yright)

root_universe.addCell(root_cell)

geom = openmoc.Geometry()
geom.setRootUniverse(root_universe)

plt.plot_cells(geom)
plt.plot_materials(geom)


#######################################
# Dataframes
#######################################
# xs_types = ['nu-fission', 'transport', 'chi', 'scatter matrix']
transport  = mesh_lib.get_mgxs(domain = mesh, mgxs_type = "transport")
chi        = mesh_lib.get_mgxs(domain = mesh, mgxs_type = "chi")
scatter    = mesh_lib.get_mgxs(domain = mesh, mgxs_type = "scatter matrix")
nu_fission = mesh_lib.get_mgxs(domain = mesh, mgxs_type = "nu-fission")
mgxs_list = [transport, chi, scatter, nu_fission]

# Get the dataframes


dataframes = [mg.get_pandas_dataframe(nuclides = "sum") for mg in mgxs_list]
#nu_fission_df = transport.get_pandas_dataframe()
#transport_df = transport.get_pandas_dataframe()
#chi_df = transport.get_pandas_dataframe(nuclides = "sum")
#scatter_df = transport.get_pandas_dataframe(nuclides = "sum")






moc_geom = openmc.openmoc_compatible.get_openmoc_geometry(sp.summary.geometry)
materials = openmoc.materialize.load_openmc_mgxs_lib(mesh_lib, moc_geom)
# Generate tracks for OpenMOC
# note: increase num_azim and decrease azim_spacing for actual results (as for TREAT)
track_generator = openmoc.TrackGenerator(moc_geom, num_azim = 32, azim_spacing = 0.1)
track_generator.generateTracks()

# Run OpenMOC
solver = openmoc.CPUSolver(track_generator)
solver.computeEigenvalue()
