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
import energy_groups
from build_mesh import mesh, STATEPOINT

PLOT = True
RUN = True

# Load the Monte Carlo results
sp = openmc.StatePoint(STATEPOINT)
mesh_lib = mgxs.Library.load_from_file(filename = "treat_mesh_lib")

# Set the domain of the MGXS objects and the MeshFilter to use the Treat_Mesh instance!
# This must be done before loading the statepoint
for xstype in mesh_lib.mgxs_types:
	for domain in mesh_lib.domains:
		mg = mesh_lib.get_mgxs(domain, xstype)
		mg.domain = mesh
		for tally in mg.tallies.values():
			for filt in tally.filters:
				if isinstance(filt, openmc.MeshFilter):
					filt.mesh = mesh
					
mesh_lib.load_from_statepoint(sp)

'''
# Optional: condense energy groups
# 2 group example
two_groups = energy_groups.casmo['2-group']
two_groups.group_edges *= 1E6
mesh_lib = mesh_lib.get_condensed_library(two_groups)
'''

# Loading from the statepoint overrides the domain
# Set the MGXS domains to the Treat_Mesh again.
for xstype in mesh_lib.mgxs_types:
	for domain in mesh_lib.domains:
		mg = mesh_lib.get_mgxs(domain, xstype)
		mg.domain = mesh

mesh_lib.domains = [mesh]

x_width, y_width, z_width = (mesh.upper_right - mesh.lower_left)/mesh.dimension


#######################################
# Dataframes
#######################################


total = mesh_lib.get_mgxs(domain=mesh, mgxs_type="total")
chi = mesh_lib.get_mgxs(domain=mesh, mgxs_type="chi")
scatter = mesh_lib.get_mgxs(domain=mesh, mgxs_type="consistent nu-scatter matrix")
nu_fission = mesh_lib.get_mgxs(domain=mesh, mgxs_type="nu-fission")
fission = mesh_lib.get_mgxs(domain=mesh, mgxs_type="fission")

# Get the dataframes
mgxs_dfs = {}
mgxs_dfs['total'] = total.get_pandas_dataframe(nuclides='sum')
# Warning: memory error on nu-scatter is possible!
mgxs_dfs['nu-scatter'] = scatter.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['chi'] = chi.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['nu-fission'] = nu_fission.get_pandas_dataframe(nuclides='sum')
mgxs_dfs['fission'] = fission.get_pandas_dataframe(nuclides='sum')


# TODO: New! Get the capture rate mesh tally data
# I believe this to be the difference between "absorption" and "fission"
capture_tally = sp.get_tally(name = "mesh tally")
vals = capture_tally.get_values(scores = ["fission"])
fission_rates = vals[:, 0, 0]

#capture_rates = absorption_rates - fission_rates
#capture_rates[capture_rates == 0] = np.nan
#capture_rates.shape = mesh.dimension
#capture_rates /= np.nanmean(capture_rates)

fission_rates[fission_rates == 0] = np.nan
fission_rates.shape = mesh.dimension
fission_rates /= np.nanmean(fission_rates)

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
		m.setNumEnergyGroups(mesh_lib.num_groups)
		
		for key in mgxs_dfs:
			df = mgxs_dfs[key]
			x_df = df[df[('mesh 1', 'x')] == i + 1]
			y_at_x = x_df[x_df[("mesh 1", "y")] == j + 1]
			
			if key == "total":
				m.setSigmaT(y_at_x['mean'].values)
			elif key == "chi":
				m.setChi(y_at_x['mean'].values)
			elif key == "nu-fission":
				m.setNuSigmaF(y_at_x['mean'].values)
			elif key == "fission":
				m.setSigmaF(y_at_x['mean'].values)
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

# Spatial discretization
cells = geom.getAllMaterialCells()
for c in cells:
	cells[c].setNumSectors(8)


if PLOT:
	plt.plot_cells(geom)
	plt.plot_materials(geom)

if RUN:
	# CMFD
	cmfd = openmoc.Cmfd()
	cmfd.setSORRelaxationFactor(1.5)
	cmfd.setLatticeStructure(mesh.dimension[0], mesh.dimension[1])
	#cmfd.setGroupStructure([[1, 2, 3], [4, 5, 6, 7]])
	cmfd.setKNearest(3)
	
	geom.setCmfd(cmfd)
	
	# Generate tracks for OpenMOC
	# note: increase num_azim and decrease azim_spacing for actual results (as for TREAT)
	
	# good run:
	#track_generator = openmoc.TrackGenerator(geom, num_azim = 128, azim_spacing = 0.01)
	# quick run:
	track_generator = openmoc.TrackGenerator(geom, num_azim = 32, azim_spacing = 0.5)
	track_generator.generateTracks()
	print("Tracks generated!")
	
	# plot the flat source region
	plt.plot_flat_source_regions(geom)
	
	# Run OpenMOC
	solver = openmoc.CPUSolver(track_generator)
	solver.computeEigenvalue()
	
	# Compute eigenvalue bias with OpenMC
	keff_mc = sp.k_combined[0]
	keff_moc = solver.getKeff()
	bias = (keff_moc - keff_mc) * 1e5
	
	print('OpenMC keff:  {:1.6f} +/- {:1.6f}'.format(keff_mc, sp.k_combined[1]))
	print('OpenMOC keff: {:1.6f}'.format(keff_moc))
	print('OpenMOC bias: {:.0f} [pcm]'.format(bias))
	
	#######################################
	# MOC Mesh
	#######################################
	#
	# Create OpenMOC Mesh on which to tally fission rates
	moc_mesh = openmoc.process.Mesh()
	moc_mesh.dimension = mesh.dimension
	moc_mesh.lower_left = mesh.lower_left
	moc_mesh.upper_right = mesh.upper_right
	moc_mesh.width = moc_mesh.upper_right - moc_mesh.lower_left
	moc_mesh.width /= moc_mesh.dimension
	# Tally OpenMOC fission rates on the Mesh
	moc_fission_rates = np.array(moc_mesh.tally_fission_rates(solver))
	moc_fission_rates.shape = mesh.dimension
	moc_fission_rates = np.fliplr(moc_fission_rates)
	
	np.savetxt("moc_data/moc_fission_rates", moc_fission_rates)
	np.savetxt("moc_data/montecarlo_fission_rates", fission_rates)
	
	if PLOT:
		import plot_moc_results
		plot_moc_results.plot_reaction_rates()




