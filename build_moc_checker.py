# MOC Checkerboard
#
# Given the mesh tally results from OpenMC, take the cross sections
# and plug them into an MOC model

import openmoc
import openmoc.materialize
import openmc.openmoc_compatible
import openmoc.plotter as plt
import openmc.mgxs as mgxs
import numpy as np
from build_mesh import mesh, Treat_Mesh

PLOT = False

print(type(mesh))
# Load the Monte Carlo results
#fname = "test_model/statepoint.30.h5"
fname = "test_model/statepoint_quick.h5"
sp = openmc.StatePoint(fname)
mesh_lib = mgxs.Library.load_from_file(filename = "treat_mesh_lib")
mesh_lib.domains = [mesh]
# FIXME: this next line breaks because it thinks it has a Mesh, not a Treat_Mesh
#mesh_lib.load_from_statepoint(sp)
x_width, y_width, z_width = (mesh.upper_right - mesh.lower_left)/mesh.dimension


#######################################
# Dataframes
#######################################

# FIXME: it doesn't like the "stride" attribute for "transport"...
# FIXME: ...and it is only looking for Tally 10001, when it needs Tally 1.

#total = mesh_lib.get_mgxs(domain=mesh, mgxs_type="total")
#chi = mesh_lib.get_mgxs(domain=mesh, mgxs_type="chi")
#scatter = mesh_lib.get_mgxs(domain=mesh, mgxs_type="consistent nu-scatter matrix")
nu_fission = mesh_lib.get_mgxs(domain=mesh, mgxs_type="nu-fission")


# Get the dataframes
mgxs_dfs = {}
#mgxs_dfs['total'] = total.get_pandas_dataframe(nuclides='sum')
#mgxs_dfs['nu-scatter'] = scatter.get_pandas_dataframe(nuclides='sum')
#mgxs_dfs['chi'] = chi.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['nu-fission'] = nu_fission.get_pandas_dataframe(nuclides='sum')


#######################################
# Geometry
#######################################

# Build a checkerboard geometry in OpenMOC
lattice = openmoc.Lattice(name = 'TREAT lattice')
lattice.setWidth(x_width, y_width)

nx, ny = mesh.dimension[:2]
universes = [[None for i in range(nx)] for j in range(ny)]

for i in range(nx):
	for j in range(ny):
		
		c = openmoc.Cell()
		m = openmoc.Material()
		m.setNumEnergyGroups(2)
		
		for key in mgxs_dfs:
			df = mgxs_dfs[key]
			x_df = df[df[('mesh 1', 'x')] == i + 1]
			y_at_x = x_df[x_df[("mesh 1", "y")] == j + 1]
			
			'''
			if key == "total":
				m.setSigmaT(y_at_x['mean'].values)
			elif key == "chi":
				m.setChi(y_at_x['mean'].values)
			elif key == "nu-fission":
				m.setNuSigmaF(y_at_x['mean'].values)
			elif key == "nu-scatter":
				m.setSigmaS(y_at_x['mean'].values)
			'''
			if key == "nu-fission":
				m.setNuSigmaF(y_at_x['mean'].values)
		
		c.setFill(m)
		u = openmoc.Universe()
		u.addCell(c)
		universes[j][i] = u

lattice.setUniverses([universes])

root_universe = openmoc.Universe(name="root universe")
root_cell = openmoc.Cell(name="root cell")
root_cell.setFill(lattice)

# Make some boundaries
min_x = openmoc.XPlane(x=mesh.lower_left[0])
max_x = openmoc.XPlane(x=mesh.upper_right[0])
min_y = openmoc.YPlane(y=mesh.lower_left[1])
max_y = openmoc.YPlane(y=mesh.upper_right[1])
surfs = [min_x, max_x, min_y, max_y]
for s in surfs:
	s.setBoundaryType(openmoc.REFLECTIVE)
root_cell.addSurface(+1, min_x)
root_cell.addSurface(-1, max_x)
root_cell.addSurface(+1, min_y)
root_cell.addSurface(-1, max_y)

root_universe.addCell(root_cell)
geom = openmoc.Geometry()
geom.setRootUniverse(root_universe)

if PLOT:
	plt.plot_cells(geom)
	plt.plot_materials(geom)
