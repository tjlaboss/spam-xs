# Test script for loading the geometry

import openmc
from openmc import mgxs
import numpy as np

summ = openmc.Summary("summary.h5")
geom = summ.openmc_geometry

openmc_cells = geom.get_all_material_cells()
print("Number of cells found:", len(openmc_cells))

# Instantiate a "coarse" 2-group EnergyGroups object
two_groups = mgxs.EnergyGroups()
two_groups.group_edges = np.array([0., 0.625, 20.0e6])

# Create dictionary to store multi-group cross sections for all cells
xs_lib = {}

# Instantiate 8-group cross sections for each cell
#for cell in openmc_cells:
cell = openmc_cells[0]
if True:	# done to preserve indentation
    xs_lib[cell.id] = {}
    xs_lib[cell.id]['transport']  = mgxs.TransportXS(groups=two_groups)
    xs_lib[cell.id]['fission'] = mgxs.FissionXS(groups=two_groups)
    xs_lib[cell.id]['nu-fission'] = mgxs.NuFissionXS(groups=two_groups)
    xs_lib[cell.id]['nu-scatter'] = mgxs.NuScatterMatrixXS(groups=two_groups)
    xs_lib[cell.id]['chi'] = mgxs.Chi(groups=two_groups)

print(xs_lib[cell.id])

