# Plot 2D Cross Sections
#
# Do a heat map for a given mgxs

MGXS_TYPE = "fission"
DIRECTORY = "kinf/tmp/"
STATEPOINT = DIRECTORY + "statepoint.100.h5"

import openmc
import openmc.mgxs as mgxs
from pylab import *

sp = openmc.StatePoint(STATEPOINT)
mesh = sp.meshes[1]
mesh_lib = mgxs.Library.load_from_file(filename="treat_mesh_lib", directory="kinf/")
mesh_lib.load_from_statepoint(sp)

#mesh_lib.load_from_statepoint(sp)
#mesh_lib.domains = [mesh]
mgxs_data = mesh_lib.get_mgxs(mesh, MGXS_TYPE)
df = mgxs_data.get_pandas_dataframe()
tally26 = sp.tallies[26]
scores = tally26._scores
print(scores)
print(mgxs_data)

for g in range(1, 11 + 1):
#for g in range(11,12):
	group_df = df[df["group in"] == g]
	#uncert = group_df['std. dev.']
	xsvals = array(group_df['mean'])
	xsvals.shape = mesh.dimension
	# Set the zero xs to NaN
	indices = xsvals <= 1E-6
	xsvals[indices] = NaN
	figure()
	title("{} macro xs, group {}".format(MGXS_TYPE, g))
	imshow(xsvals.squeeze(), interpolation='none', cmap='jet')
	colorbar()

show()
