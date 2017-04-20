#############
# Export to MOC
import openmoc
import openmoc.materialize
import openmc.openmoc_compatible
import openmoc.plotter as plt

import openmc.mgxs as mgxs
import numpy as np

sp = openmc.StatePoint('statepoint.20.h5')

mesh_lib = mgxs.Library.load_from_file(filename = "mesh_lib")
mesh_lib.load_from_statepoint(sp)
mesh = mesh_lib.domains[0]
x_width, y_width = (mesh.upper_right - mesh.lower_left)/mesh.dimension

#######################################
# Dataframes
#######################################

total = mesh_lib.get_mgxs(domain=mesh, mgxs_type="total")
chi = mesh_lib.get_mgxs(domain=mesh, mgxs_type="chi")
scatter = mesh_lib.get_mgxs(domain=mesh, mgxs_type="consistent nu-scatter matrix")
nu_fission = mesh_lib.get_mgxs(domain=mesh, mgxs_type="nu-fission")

# Get the dataframes
mgxs_dfs = {}
mgxs_dfs['total'] = total.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['nu-fission'] = nu_fission.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['nu-scatter'] = scatter.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['chi'] = chi.get_pandas_dataframe(nuclides='sum')

#######################################
# Geometry
#######################################

# Build a checkerboard geometry in OpenMOC
lattice = openmoc.Lattice(name='4x4 lattice')
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
			x_df = df[df[('mesh 10000', 'x')] == i + 1]
			y_at_x = x_df[x_df[("mesh 10000", "y")] == j + 1]
			
			if key == "total":
				m.setSigmaT(y_at_x['mean'].values)
			elif key == "chi":
				m.setChi(y_at_x['mean'].values)
			elif key == "nu-fission":
				m.setNuSigmaF(y_at_x['mean'].values)
			elif key == "nu-scatter":
				m.setSigmaS(y_at_x['mean'].values)
		
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

plt.plot_cells(geom)
plt.plot_materials(geom)

# Generate tracks for OpenMOC
# note: increase num_azim and decrease azim_spacing for actual results (as for TREAT)
track_generator = openmoc.TrackGenerator(geom, num_azim=32, azim_spacing=0.1)
track_generator.generateTracks()

# Run OpenMOC
solver = openmoc.CPUSolver(track_generator)
solver.computeEigenvalue()

# Compute eigenvalue bias with OpenMC
keff_mc = sp.k_combined[0]
keff_moc = solver.getKeff()
bias = (keff_moc - keff_mc) * 1e5

print('OpenMC keff: {:1.6f} +/- {:1.6f}'.format(keff_mc, sp.k_combined[1]))
print('OpenMC keff: {:1.6f}'.format(keff_moc))
print('OpenMC bias: {:.0f} [pcm]'.format(bias))

