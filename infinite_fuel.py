# Infinite fuel
#
# Used for k-inf calculations of a fuel element

import openmc
import openmc.mgxs as mgxs
import energy_groups

EXPORT = True
DESTINATION = "treat2d/kinf/"
ADDITIONAL_STATEPOINTS = 3
STATEPOINT_INTERVAL = 20
TALLY_MGXS = True
MESH_DIVISIONS = 40

# Extract the fuel element geometry from an existing summary
summ = openmc.Summary("treat2d/summary.h5")
geom = summ.geometry
element = geom.get_all_cells()[20110]
# Get and set the boundaries for the cell
core_lat = geom.get_all_lattices()[100]
px, py = core_lat.pitch[0:2]
left = openmc.XPlane(x0=-px/2, boundary_type="reflective")
right = openmc.XPlane(x0=px/2, boundary_type="reflective")
front = openmc.YPlane(y0=px/2, boundary_type="reflective")
back = openmc.YPlane(y0=-py/2, boundary_type="reflective")
bot = openmc.ZPlane(z0=-1.0, boundary_type="reflective")
top = openmc.ZPlane(z0=+1.0, boundary_type="reflective")
element.region = +left & -right & +back & -front & +bot & -top
# Let that be the root universe
fuel_verse = openmc.Universe(name="Individual element universe")
fuel_verse.add_cell(element)
geometry_xml = openmc.Geometry()
geometry_xml.root_universe = fuel_verse

mats = element.fill.get_all_materials()
materials_xml = openmc.Materials(mats.values())

settings_xml = openmc.Settings()
settings_xml.batches = 100
settings_xml.inactive = 35
settings_xml.particles = int(1E7)
if ADDITIONAL_STATEPOINTS:
	n = ADDITIONAL_STATEPOINTS + 1
	batch_nos = [None]*n
	for i in range(n):
		batch_nos[-1 - i] = settings_xml.batches - i*STATEPOINT_INTERVAL
	settings_xml.statepoint = {"batches": batch_nos}

plots_xml = openmc.Plots()
for axis in ("xy", "yz"):
	new_plot = openmc.Plot()
	new_plot.basis = axis
	new_plot.color_by = "material"
	new_plot.pixels = (400, 400)
	if axis == "xy":
		new_plot.width = [px, py]
	elif axis == "yz":
		new_plot.width = [py, 2]
	plots_xml.append(new_plot)

if TALLY_MGXS:
	mesh_lib = mgxs.Library(geom)
	groups = mgxs.EnergyGroups()
	# RattleSNake normally uses 11 energy groups
	groups.group_edges = energy_groups.treat["11-group"].group_edges*1E6
	mesh_lib.energy_groups = groups
	# The four most important cross sections to tally right now
	mesh_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'capture', 'chi', 'consistent nu-scatter matrix']
	mesh_lib.by_nuclide = False
	mesh_lib.domain_type = "mesh"
	mesh_lib.correction = None
	# Build the tally mesh
	mesh = openmc.Mesh(1, "Elemental Mesh")
	mesh.lower_left  = (-px/2, -px/2)
	mesh.upper_right = (+px/2, +px/2)
	mesh.dimension = (MESH_DIVISIONS, MESH_DIVISIONS)
	# Finalize for tallies
	mesh_lib.domains = [mesh]
	mesh_lib.build_library()
	tallies_xml = openmc.Tallies()
	mesh_lib.add_to_tallies_file(tallies_xml, merge=True)


if EXPORT:
	file_dict = {"geometry.xml" : geometry_xml,
	             "materials.xml": materials_xml,
	             "settings.xml" : settings_xml,
	             "plots.xml"    : plots_xml,
	             }
	if TALLY_MGXS:
		file_dict["tallies.xml"] = tallies_xml
		mesh_lib.dump_to_file("treat_mesh_lib", DESTINATION)
	for filename, xml in file_dict.items():
		xml.export_to_xml(DESTINATION + filename)
	print("\nExported to", DESTINATION)
