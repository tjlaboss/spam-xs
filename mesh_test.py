# Mesh test
#
# Learning how to use a mesh tally

import openmc
from openmc import mgxs
import numpy as np
import matplotlib.pyplot as plt


# Extract the geometry from an existing summary
summ = openmc.Summary("summary.h5")
geom = summ.openmc_geometry
mats = summ.materials

# 2-group approximation
two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])
mg_lib = mgxs.Library(geom)
mg_lib.energy_groups = two_groups

# For purposes of this demonstration, let's just look at the capture
# and the transport cross sections
mg_lib.mgxs_types = ["transport", "capture"]
mg_lib.by_nuclide = True
mg_lib.domain_type = "material"
#mg_lib.domain_type = "mesh"
mg_lib.build_library()

# Define a mesh
# Use the core lattice as a template
core_lat = summ.get_lattice_by_id(100)
test_mesh = openmc.Mesh(name = "Coarse mesh (test)")
test_mesh.lower_left = core_lat.lower_left
print(test_mesh.lower_left)
test_mesh.upper_right = [0, 0, 0]
test_mesh.type = "regular"
#test_mesh.upper_right = tuple(map(lambda l: -l, core_lat.lower_left))
test_mesh.dimension = [10, 10, 10]
#test_mesh.dimension = core_lat.universes.shape
#test_mesh.width = [10, 10, 50]
print(test_mesh)
print(type(test_mesh))

# Filter
mesh_filter = openmc.Filter(test_mesh)
# Tallies
mesh_tally = openmc.Tally(name = "mesh tally")
mesh_tally.scores = mg_lib.mgxs_types
mesh_tally.filters = [mesh_filter,]



tallies_xml = openmc.Tallies()
mg_lib.add_to_tallies_file(tallies_file = tallies_xml, merge = True)
tallies_xml.append(mesh_tally)
tallies_xml.export_to_xml()



