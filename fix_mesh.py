# Mesh test
#
# Learning how to use a mesh tally

import openmc
from openmc import mgxs
import numpy as np


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

if __name__ == "__main__":
	# Add tally to collection
	tallies_file.append(tally)
	# Filter
	mesh_filter = openmc.Filter(mesh)
	
	tallies_file.export_to_xml("test_model/tallies.xml")
	
	# Examine the data after the run
	sp = openmc.StatePoint('test_model/statepoint.50.h5')
	mg_lib.load_from_statepoint(sp)
	# print(mg_lib.get_mgxs(mesh, "fission"))
	print(mg_lib.all_mgxs.keys())

	#fission_mgxs = mg_lib.get_mgxs(mesh, "fission")
	fission_mgxs = mg_lib.get_mgxs(fuel, "fission")
	print(type(fission_mgxs))
	# fission_mgxs = mg_lib.all_mgxs[mesh.id]["fission"]
	fission_df = fission_mgxs.get_pandas_dataframe()
	#print(type(fission_df))
	print(fission_df)
	fission_mgxs.print_xs()
	



