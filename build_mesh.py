# Mesh test
#
# Learning how to use a mesh tally


import openmc
from openmc import mgxs
import pylab
import energy_groups
from treat_mesh import Treat_Mesh
from copy import deepcopy

# Settings
EXPORT = True
PLOT = True
STATEPOINT = 'treat2d/statepoint_8groups.h5'

# Extract the geometry from an existing summary
summ = openmc.Summary("summary.h5")
geom = summ.geometry
mesh_lib = mgxs.Library(geom)
mats = geom.get_all_materials()
fuel = mats[90000]

# 8 Energy Groups
groups = mgxs.EnergyGroups()
# Convert from MeV to eV
groups.group_edges = energy_groups.casmo["8-group"].group_edges*1E6
mesh_lib.energy_groups = groups



# The four most important cross sections to tally right now
mesh_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'capture', 'chi', 'consistent nu-scatter matrix']
mesh_lib.by_nuclide = True
mesh_lib.domain_type = "mesh"
mesh_lib.correction = None

# Material
material_lib = mgxs.Library(geom)
material_lib.energy_groups = groups
material_lib.mgxs_types = mesh_lib.mgxs_types
material_lib.domain_type = "material"
material_lib.domains = mats.values()
material_lib.by_nuclide = False
material_lib.build_library()

# Define a mesh
# Instantiate a tally Mesh
mesh = Treat_Mesh(mesh_id = 1, geometry = geom)
mesh.mesh_size = (1, 1, 1)
#mesh.mesh_size = (2, 2, 1)

# FIXME: For some reason, in tallies.xml, the mesh gets exported 4 times!!
core_lat = geom.get_all_lattices()[100]
xdist = -core_lat.lower_left[0]
zbot = mesh._surfaces[20009].z0  # bottom of active fuel region
ztop = mesh._surfaces[20010].z0  # top of active fuel

mesh.lower_left = deepcopy(core_lat.lower_left)
mesh.lower_left[-1] = zbot
mesh.upper_right = -deepcopy(core_lat.lower_left)
mesh.upper_right[-1] = ztop
mesh.type = 'regular'
mesh.dimension = deepcopy(core_lat.shape)

mesh_lib.domains = [mesh]
mesh_lib.build_library()
# Turn off by_nuclide for nu-scatter
cnsm_mgxs = mesh_lib.get_mgxs(mesh, 'consistent nu-scatter matrix')
cnsm_mgxs.by_nuclide = False

mesh_lib.dump_to_file("treat_mesh_lib")
material_lib.dump_to_file("treat_material_lib")

def make_tallies():
	# Instantiate tally Filter
	mesh_filter = openmc.MeshFilter(mesh)
	
	# Instantiate the Tally
	fission_tally = openmc.Tally(name = 'mesh tally')
	fission_tally.filters = [mesh_filter]
	fission_tally.scores = ["fission"]
	
	capture_tally = openmc.Tally(name = "U238 capture tally")
	capture_tally.filters = [mesh_filter]
	capture_tally.scores = ["absorption", "fission"]
	capture_tally.nuclides = ["U238"]
	
	# Create a "tallies.xml" file for the MGXS Library
	tallies_file = openmc.Tallies()
	tallies_file.extend([fission_tally, capture_tally])
	
	mesh_lib.add_to_tallies_file(tallies_file, merge = True)
	material_lib.add_to_tallies_file(tallies_file, merge = True)
	return tallies_file


def plot_mgxs(nuc, xstype, xs_df, g, groups, x0 = -xdist, x1 = xdist, n = mesh.dimension[1]):
	"""Plotting a single energy group as a function of space
	
	Inputs:
		nuc:        str; name of nuclide
		xstype:     str; name of reaction type
		xs_df:      instance of pandas dataframe containing the cross sections
		g:          int; group number
		
		x0:         float, optional; x-value to start plotting at
					[Default: `lower_left` x coordinate of Treat lattice]
		x1:         float, optional; x-value to stop plotting at
					[Default: `upper_right` x coordinate of Treat lattice]
		n:          int, optional; number of values to plot
					[Default: x `dimension` of Treat lattice]
	
	Outputs:
		None
	"""
	row = int(pylab.ceil(n/2))
	xlist = pylab.linspace(x0, x1, n)
	xs_scale = "macro"
	# nuc_xs = xs_lib.get_xs(order_groups = "decreasing", xs_type = xs_scale, groups = g).
	pitch = (x1 - x0)/(n - 1)
	
	group_df = xs_df[xs_df["group in"] == g]
	y_df = group_df[group_df[('mesh 1', 'y')] == row]
	# y_at_z = y_df[y_df[("mesh 1", "z")] == row]
	yvals = y_df['mean']
	uncert = y_df['std. dev.']
	
	# plotting stuff
	ylist = pylab.array(yvals)
	ulist = pylab.array(uncert)
	xtvals = pylab.linspace(x0, x1, n)
	
	if n % 2:
		# Odd number: assemblies are offset by a halfwidth
		style = "steps-mid"
		xtvals -= pitch/2
	else:
		# Even: assemblies are aligned with default grid
		style = "steps"
	
	pylab.grid()
	pylab.xticks(xtvals)
	pylab.xlim(min(xlist) - pitch, max(xlist) + pitch)
	
	pylab.plot(xlist, ylist + ulist, "red", drawstyle = style, alpha = 0.5, label = "+/- 1sigma")
	pylab.plot(xlist, ylist - ulist, "red", drawstyle = style, alpha = 0.5)
	pylab.plot(xlist, ylist, drawstyle = style, label = "$\Sigma$")
	
	pylab.legend(loc = "lower center")
	title_string = "{} {}scopic Cross Section for {}".format(xstype.title(), xs_scale.title(), nuc)
	pylab.xlabel("Radial distance (cm)")
	pylab.ylabel("$\Sigma$ (cm$^{-1}$)")
	pylab.title(title_string, {"fontsize": 14})
	pylab.suptitle("Group {} of {}".format(g, groups.num_groups))
	pylab.show()


if __name__ == "__main__":
	# Add tally to collection
	tallies_xml = make_tallies()
	if EXPORT:
		tallies_xml.export_to_xml("treat2d/tallies.xml")
	
	# Examine the data after the run
	sp = openmc.StatePoint(STATEPOINT)
	
	mesh_lib.load_from_statepoint(sp)
	mesh_lib.domains = [mesh]
	# Reassign the loaded data to be on the Treat_Mesh
	for domain in mesh_lib.domains:
		for mgxs_type in mesh_lib.mgxs_types:
			xs = mesh_lib.get_mgxs(domain, mgxs_type)
			xs.domain = mesh
	
	#nuc = "C0"
	#xstype = "capture"
	nuc = "U235"
	xstype = "nu-fission"
	fission_mgxs = mesh_lib.get_mgxs(mesh, xstype)
	fission_df = fission_mgxs.get_pandas_dataframe(nuclides = [nuc])
	
	if PLOT:
		# Plot stuff
		plot_mgxs(nuc, xstype, fission_df, 7, groups)
