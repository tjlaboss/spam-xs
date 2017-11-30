# Infinite fuel
#
# Used for k-inf calculations of a fuel element

import openmc

EXPORT = True
DESTINATION = "treat2d/kinf/"
ADDITIONAL_STATEPOINTS = 3
STATEPOINT_INTERVAL = 20

# Extract the fuel element geometry from an existing summary
summ = openmc.Summary("treat2d/summary.h5")
geom = summ.geometry
element = geom.get_all_cells()[20110]
# Get and set the boundaries for the cell
core_lat = geom.get_all_lattices()[100]
px, py = core_lat.pitch[0:2]
left = openmc.XPlane(x0 = -px/2, boundary_type = "reflective")
right = openmc.XPlane(x0 = px/2, boundary_type = "reflective")
front = openmc.YPlane(y0 = px/2, boundary_type = "reflective")
back = openmc.YPlane(y0 = -py/2, boundary_type = "reflective")
bot = openmc.ZPlane(z0 = -1.0, boundary_type = "reflective")
top = openmc.ZPlane(z0 = +1.0, boundary_type = "reflective")
element.region = +left & -right & +back & -front & +bot & -top
# Let that be the root universe
fuel_verse = openmc.Universe(name = "Individual element universe")
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
		batch_nos[-1-i] = settings_xml.batches - i*STATEPOINT_INTERVAL
	print(batch_nos)
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

if EXPORT:
	file_dict = {"geometry.xml": geometry_xml,
	             "materials.xml": materials_xml,
	             "settings.xml": settings_xml,
	             "plots.xml": plots_xml,
	             }
	for filename, xml in file_dict.items():
		xml.export_to_xml(DESTINATION + filename)
	print("\nExported to", DESTINATION)
