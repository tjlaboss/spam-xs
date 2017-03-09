# Mesh test
#
# Learning how to use a mesh tally

import openmc
from openmc import mgxs
import numpy as np
import pylab

# Settings
EXPORT = False
PLOT = True


# Extract the geometry from an existing summary
summ = openmc.Summary("summary.h5")
geom = summ.openmc_geometry
#mats = geom.get_all_materials()
fuel = summ.get_material_by_id(90000)

# 2-group approximation
two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])
mg_lib = mgxs.Library(geom)
mg_lib.energy_groups = two_groups

# For purposes of this demonstration, let's just look at the capture
# and the transport cross sections
mg_lib.mgxs_types = ['fission', 'nu-fission', 'transport']
mg_lib.by_nuclide = False
mg_lib.domain_type = "material"
#mg_lib.domain_type = "mesh"

# Define a mesh
# Instantiate a tally Mesh
mesh = openmc.Mesh(mesh_id=1)
# Use the core lattice as a template
core_lat = summ.get_lattice_by_id(100)
xdist = -core_lat.pitch[0]*core_lat.lower_left[0]

mesh.lower_left = core_lat.lower_left
mesh.upper_right = [0, 0, 0]
mesh.type = 'regular'
mesh.dimension = [17, 17, 17]

#mesh2 = copy.deepcopy(mesh)
#mesh2.id = 2

#mg_lib.domains = [mesh,]
mg_lib.build_library()

# Instantiate tally Filter
mesh_filter = openmc.MeshFilter(mesh)

# Instantiate the Tally
tally = openmc.Tally(name='mesh tally')
tally.filters = [mesh_filter]
tally.scores = ["fission", "nu-fission"]

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()
mg_lib.add_to_tallies_file(tallies_file, merge=True)




def plot_mgxs(nuc, xs_lib, xs_df, g, groups, x0 = 0, x1 = xdist, n = 17):
	"""Plotting a single energy group as a function of space
	
	Inputs:
		nuc:        str; name of nuclide
		xs_lib:     instance of mgxs.MGXS
		xs_df:      instance of pandas dataframe
		g:          int; group number
		groups:     instance of mgxs.EnergyGroups()
	
	Outputs:
		???
	"""
	
	xlist = pylab.linspace(x0, x1, n)
	xs_scale = "macro"
	#nuc_xs = xs_lib.get_xs(order_groups = "decreasing", xs_type = xs_scale, groups = g)
	
	# FIXME: This returns an empty array
	nuc_df = xs_df[xs_df['nuclide'] == nuc]['mean']
	print(nuc_df)
	nuc_matrix = nuc_df.as_matrix()
	print(nuc_matrix)
	
	# plotting stuff for later
	ylist = xlist #debug
	pylab.plot(xlist, ylist, drawstyle = "steps")
	title_string = "{0} {1}scopic Cross Section".format(nuc, xs_scale.title())
	pylab.xlabel("Radial distance (cm)")
	pylab.ylabel("$\Sigma$ (cm$^{-1}$)")
	pylab.title(title_string, {"fontsize":14})
	#pylab.show()
	




if __name__ == "__main__":
	# Add tally to collection
	tallies_file.append(tally)
	# Filter
	mesh_filter = openmc.Filter(mesh)
	
	if EXPORT:
		tallies_file.export_to_xml("test_model/tallies.xml")
	
	# Examine the data after the run
	sp = openmc.StatePoint('test_model/statepoint.50.h5')
	mg_lib.load_from_statepoint(sp)
	print(mg_lib.all_mgxs.keys())
	
	"""
	I think the main problem is that fission_mgxs should be
	on a mesh, not on a material.
	
	I'm not sure how to plot the fission cross section for a group
	as a function of radial position without this.
	"""
	fission_mgxs = mg_lib.get_mgxs(fuel, "fission")
	# fission_mgxs = mg_lib.get_mgxs(mesh, "fission")
	# fission_mgxs = mg_lib.all_mgxs[mesh.id]["fission"]
	fission_df = fission_mgxs.get_pandas_dataframe()
	
	fission_mgxs.print_xs()
	
	if PLOT:
		# Plot stuff
		plot_mgxs("U235", fission_mgxs, fission_df, 0, two_groups)
	



