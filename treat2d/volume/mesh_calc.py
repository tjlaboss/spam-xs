import sys; sys.path.append("../..")
from copy import deepcopy
from collections import OrderedDict
import numpy
import openmc
import openmc.mgxs as mgxs
import energy_groups
from treat_mesh import Treat_Mesh

# Universe IDs
FUEL = 9
CRD = 5
ZR_DUMMY = 3
AL_DUMMY = 4

# Settings
EXPORT = False
RUN = False
MESH_DIVISIONS = 3

summ = openmc.Summary("../summary.h5")
geom = summ.geometry

# Settings
num = MESH_DIVISIONS*19
ROOT = 'treat2d/{0}x{0}/'.format(num)
STATEPOINT = ROOT + 'statepoint_11groups.h5'

#
lat = geom.get_all_lattices()[100]
universes = geom.get_all_universes()
cells = geom.get_all_cells()
fuel_verse1 = universes[9]
fuel_verse2 = universes[98]

# Instantiate a tally Mesh
mesh = Treat_Mesh(1, geometry=geom)
mesh.mesh_size = (MESH_DIVISIONS, MESH_DIVISIONS, 1)
xdist = -lat.lower_left[0]
zbot = mesh._surfaces[20009].z0  # bottom of active fuel region
ztop = mesh._surfaces[20010].z0  # top of active fuel
mesh.lower_left = deepcopy(lat.lower_left)
mesh.lower_left[-1] = zbot
mesh.upper_right = -deepcopy(lat.lower_left)
mesh.upper_right[-1] = ztop
mesh.type = 'regular'
mesh.dimension = deepcopy(lat.shape)

'''
mesh_lib = mgxs.Library(geom)
groups = mgxs.EnergyGroups()
groups.group_edges = energy_groups.treat["11-group"].group_edges*1E6
mesh_lib.energy_groups = groups
mesh_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'capture', 'chi', 'consistent nu-scatter matrix']
mesh_lib.by_nuclide = True
mesh_lib.domain_type = "mesh"
mesh_lib.correction = None
cnsm_mgxs = mesh_lib.get_mgxs(mesh, 'consistent nu-scatter matrix')
cnsm_mgxs.by_nuclide = False
'''

# OK: Now try volume calculating on the mesh
#universes_to_tally = list(set(lat.universes.flatten()))
ll0 = mesh.lower_left/mesh.dimension
ur0 = mesh.upper_right/mesh.dimension
verse_coordinates = OrderedDict({FUEL: (0, 0),
			                     CRD: (4, 2),
			                     ZR_DUMMY: (0, 6),
			                     AL_DUMMY: (0, 7)})
all_calcs = OrderedDict()
px, py, pz = lat.pitch
for uid in verse_coordinates:
	print("\nUniverse", uid)
	nx, ny = verse_coordinates[uid]
	x0 = round((nx - 0.5)*px, 4)
	x1 = round((nx + 0.5)*px, 4)
	y0 = round((ny - 0.5)*py, 4)
	y1 = round((ny + 0.5)*py, 4)
	ll = (x0, y0, -0.5)
	ur = (x1, y1, +0.5)
	all_calcs[uid] = (openmc.VolumeCalculation([universes[uid]], int(1E4), ll, ur))
	
if EXPORT:
	settings = openmc.Settings()
	settings.volume_calculations = list(all_calcs)
	settings.export_to_xml()
if RUN:
	openmc.calculate_volumes(threads=12)

for u in verse_coordinates:
	all_calcs[u].load_results("volume_u{}.h5".format(u))

print(all_calcs[FUEL].atoms[FUEL].keys())
