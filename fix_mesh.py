# Mesh test
#
# Learning how to use a mesh tally


import openmc
from openmc import mgxs
import numpy as np
import pylab
from treat_mesh import Treat_Mesh

# Settings
EXPORT = True
PLOT = True


# Extract the geometry from an existing summary
summ = openmc.Summary("summary.h5")
geom = summ.geometry
mats = geom.get_all_materials()
fuel = mats[90000]

# 2-group approximation
two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])
mesh_lib = mgxs.Library(geom)
mesh_lib.energy_groups = two_groups

# For purposes of this demonstration, let's just look at the capture
# and the transport cross sections
mesh_lib.mgxs_types = ['fission', 'nu-fission', 'transport']
mesh_lib.by_nuclide = True
#mesh_lib.nuclides = ["U235", "U238"]
# TODO: select all the nuclides from `mats`
mesh_lib.domain_type = "mesh"


# Define a mesh
# Instantiate a tally Mesh
mesh = Treat_Mesh(mesh_id=1)
# Use the core lattice as a template
#core_lat = summ.get_lattice_by_id(100)
core_lat = geom.get_all_lattices()[100]
xdist = -core_lat.pitch[0]*core_lat.lower_left[0]

mesh.lower_left = core_lat.lower_left
mesh.upper_right = [0, 0, 0]
mesh.type = 'regular'
mesh.dimension = [17, 17, 17]

mesh_lib.domains = [mesh,]
mesh_lib.build_library()

# Instantiate tally Filter
mesh_filter = openmc.MeshFilter(mesh)


# Create the material lib
material_lib = mgxs.Library(geom)
material_lib.energy_groups = two_groups
material_lib.mgxs_types = ['fission', 'nu-fission', 'transport']
material_lib.domain_type = "material"
material_lib.domains = mats.values()
material_lib.by_nuclide = True
material_lib.build_library()

# Instantiate the Tally
tally = openmc.Tally(name='mesh tally')
tally.filters = [mesh_filter]
tally.scores = ["fission", "nu-fission"]

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()


mesh_lib.add_to_tallies_file(tallies_file, merge=True)
material_lib.add_to_tallies_file(tallies_file, merge = True)
print(material_lib.all_mgxs)



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
	
	xlist = pylab.linspace(x1, x0, n)
	xs_scale = "macro"
	#nuc_xs = xs_lib.get_xs(order_groups = "decreasing", xs_type = xs_scale, groups = g).
	
	group_df = xs_df[xs_df["group in"] == g]
	y_df = group_df[group_df[('mesh 1', 'y')] == 17]
	y_at_z = y_df[y_df[("mesh 1", "z")] == 17]
	yvals = y_at_z['mean'].values
	print(yvals)
	
	# FIXME: This returns an empty array
	#nuc_df = xs_df[xs_df['nuclide'] == nuc]['mean']
	#print(nuc_df)
	#nuc_matrix = nuc_df.as_matrix()
	#print(nuc_matrix)
	
	# plotting stuff for later
	ylist = yvals
	pylab.plot(xlist, ylist, drawstyle = "steps")
	title_string = "{0} {1}scopic Cross Section".format("Total", xs_scale.title())
	pylab.xlabel("Radial distance (cm)")
	pylab.ylabel("$\Sigma$ (cm$^{-1}$)")
	pylab.title(title_string, {"fontsize":14})
	pylab.show()
	




if __name__ == "__main__":
	# Add tally to collection
	tallies_file.append(tally)
	# Filter
	mesh_filter = openmc.Filter(mesh)
	
	if EXPORT:
		tallies_file.export_to_xml("test_model/tallies.xml")
	
	# Examine the data after the run
	sp = openmc.StatePoint('test_model/statepoint.50.h5')
	mesh_lib.load_from_statepoint(sp)
	material_lib.load_from_statepoint(sp)
	
	"""
	I think the main problem is that fission_mgxs should be
	on a mesh, not on a material.
	
	I'm not sure how to plot the fission cross section for a group
	as a function of radial position without this.
	"""
	#fission_mgxs = material_lib.get_mgxs(fuel, "fission")
	fission_mgxs = mesh_lib.get_mgxs(mesh, "fission")
	# fission_mgxs = mg_lib.get_mgxs(mesh, "fission")
	# fission_mgxs = mg_lib.all_mgxs[mesh.id]["fission"]
	fission_df = fission_mgxs.get_pandas_dataframe()
	#print(fission_df.head(17), fission_df.tail(17))
	#fission_mgxs.print_xs()
	
	if PLOT:
		# Plot stuff
		plot_mgxs("U235", fission_mgxs, fission_df, 2, two_groups)
	



