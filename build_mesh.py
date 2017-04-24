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
# mesh_lib.nuclides = ["U235", "U238"]
# TODO: select all the nuclides from `mats`
mesh_lib.domain_type = "mesh"

# Define a mesh
# Instantiate a tally Mesh
mesh = Treat_Mesh(mesh_id = 1, geometry = geom)
# FIXME: For some reason, in tallies.xml, the mesh gets exported 3 times!!
core_lat = geom.get_all_lattices()[100]
xdist = -core_lat.lower_left[0]
zbot = mesh.surfaces[20009].z0  # bottom of active fuel region
ztop = mesh.surfaces[20010].z0  # top of active fuel

mesh.lower_left = core_lat.lower_left
mesh.lower_left[-1] = zbot
mesh.upper_right = -core_lat.lower_left
mesh.upper_right[-1] = ztop
mesh.type = 'regular'
mesh.dimension = core_lat.shape

mesh_lib.domains = [mesh]
mesh_lib.build_library()
mesh_lib.dump_to_file("treat_mesh_lib")

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


def make_tallies():
	# Instantiate the Tally
	tally = openmc.Tally(name = 'mesh tally')
	tally.filters = [mesh_filter]
	tally.scores = ["fission", "nu-fission"]
	
	# Create a "tallies.xml" file for the MGXS Library
	tallies_file = openmc.Tallies()
	tallies_file.append(tally)
	
	mesh_lib.add_to_tallies_file(tallies_file, merge = True)
	material_lib.add_to_tallies_file(tallies_file, merge = True)
	return tallies_file


def plot_mgxs(nuc, xs_lib, xs_df, g, groups, x0 = -xdist, x1 = xdist, n = 19):
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
	# nuc_xs = xs_lib.get_xs(order_groups = "decreasing", xs_type = xs_scale, groups = g).
	
	group_df = xs_df[xs_df["group in"] == g]
	y_df = group_df[group_df[('mesh 1', 'y')] == 9]
	# y_at_z = y_df[y_df[("mesh 1", "z")] == 9]
	yvals = y_df['mean']
	uncert = y_df['std. dev.']
	
	# FIXME: This returns an empty array
	# nuc_df = xs_df[xs_df['nuclide'] == nuc]['mean']
	# print(nuc_df
	# nuc_matrix = nuc_df.as_matrix()
	# print(nuc_matrix)
	
	# plotting stuff for later
	ylist = yvals
	pylab.plot(xlist, ylist, drawstyle = "steps", label = "$\Sigma$")
	pylab.plot(xlist, ylist + uncert, "red", drawstyle = "steps", alpha = 0.5, label = "+/- 1sigma")
	pylab.plot(xlist, ylist - uncert, "red", drawstyle = "steps", alpha = 0.5)
	pylab.legend(loc = "upper left")
	pylab.grid()
	pylab.xticks(pylab.linspace(x0, x1, n))
	title_string = "{0} {1}scopic Cross Section".format("Total", xs_scale.title())
	pylab.xlabel("Radial distance (cm)")
	pylab.ylabel("$\Sigma$ (cm$^{-1}$)")
	pylab.title(title_string, {"fontsize": 14})
	pylab.show()


if __name__ == "__main__":
	# Add tally to collection
	tallies_file = make_tallies()
	if EXPORT:
		tallies_file.export_to_xml("test_model/tallies.xml")
	
	# Filter
	mesh_filter = openmc.Filter(mesh)
	
	# Examine the data after the run
	# sp = openmc.StatePoint('test_model/statepoint.50.h5')
	sp = openmc.StatePoint('test_model/statepoint.30.h5')
	
	mesh_lib.load_from_statepoint(sp)
	mesh_lib.domains = [mesh]
	for domain in mesh_lib.domains:
		for mgxs_type in mesh_lib.mgxs_types:
			xs = mesh_lib.get_mgxs(domain, mgxs_type)
			xs.domain = mesh
	
	material_lib.load_from_statepoint(sp)
	
	fission_mgxs = mesh_lib.get_mgxs(mesh, "transport")
	fission_df = fission_mgxs.get_pandas_dataframe(nuclides = ["C0"])
	# print(fission_df.head(17), fission_df.tail(17))
	# fission_mgxs.print_xs()
	
	if PLOT:
		# Plot stuff
		plot_mgxs("U235", fission_mgxs, fission_df, 1, two_groups)
