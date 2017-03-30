# Area calculator for 2D TREAT lattice cells


import openmc
from math import sqrt, pi

rt = sqrt(2)/2  # "root two" (useful shorthand)

# Global variables
# Extract the geometry from an existing summary
geometry = openmc.Summary("summary.h5").geometry


def __setup(geom):
	"""Get some of the basics out of the way
	
	Inputs:
		:param geom: instance of openmc.Geometry for the TREAT model
	
	Outputs:
		:return surfs:    OrderedDict of all the surfaces in `geom`
		:return key_list: tuple of the strings for each of the directions
		:return n:        int; length of `key_list` (should be 8)
		:return xpitch:   float; assembly pitch in the x-direction (cm)
		:return ypitch:   float; assembly pitch in the y-direction (cm)
	"""
	# The main core lattice which the cells appear in
	core_lat = geom.get_all_lattices()[100]
	xpitch, ypitch = core_lat.pitch[0:2]
	
	# Create dictionaries for each region
	key_list = ("e", "s", "w", "n", "se", "sw", "nw", "ne")
	n = len(key_list)
	
	# Get an OrderedDict of all the surfaces in the geometry
	surfs = geometry.get_all_surfaces()
	return surfs, key_list, n, xpitch, ypitch


def fuel_cell_by_material(geom, display = False):
	"""Calculate the area of each material a fuel lattice cell
	
	Inputs:
		:param geom: instance of openmc.Geometry for the TREAT model
		:param display: Boolean; whether to print the areas for each region
						[Default: False]
	
	Outputs:
		:return: fuel_area, gap_area, clad_area, outer_area
				 Areas in cm^2 for the regions indicated
	"""
	surfs, key_list, n, xpitch, ypitch = __setup(geom)
	
	fuel_base = 90001
	fuel_list = [None, ]*n
	for i in range(n):
		fuel_list[i] = surfs[fuel_base + i]
	fuel = dict(zip(key_list, fuel_list))
	
	gap_base = 90011
	gap_list = [None, ]*n
	for i in range(n):
		gap_list[i] = surfs[gap_base + i]
	gap = dict(zip(key_list, gap_list))
	
	clad_base = 90021
	clad_list = [None, ]*n
	for i in range(n):
		clad_list[i] = surfs[clad_base + i]
	clad = dict(zip(key_list, clad_list))
	
	# Get the surface area (which is all this code does now)
	# Corners. Subscript 'm' = midpoint, 'c' = corner
	"""
	._______
	|c\    /
	|  \./
	|  / ^ m
	|/
	|

	Let "h" be the distance from the corner to the midpoint: the
	perpendicular bisector of the right angle in the corner.
		h = sqrt(x^2 + y^2) - midpoint

	Area of this triangle:
		A = h^2
	"""
	
	# Start with the outer gap (cooling air)
	# Total area of the cell
	cell_area = xpitch*ypitch
	
	# Fuel area: Calculate the middle area and subtract the corners
	# Sides (for middle)
	right = fuel["e"].x0
	left = fuel["w"].x0
	top = fuel["n"].y0
	bottom = fuel["s"].y0
	fuel_middle_area = (right - left)*(top - bottom)
	# Corner
	m = fuel["nw"].coefficients['D']*rt
	c = sqrt(top**2 + left**2)
	h = c - m
	fuel_corner_area = h**2
	
	# Air gap calculation: same as fuel, but will ultimately subtract fuel
	# Middle
	left = gap["w"].x0
	top = gap["n"].y0
	right = gap["e"].x0
	bottom = gap["s"].y0
	gap_middle_area = (right - left)*(top - bottom)
	# Corner
	m = gap["nw"].coefficients['D']*rt
	c = sqrt(top**2 + left**2)
	h = c - m
	gap_corner_area = h**2
	
	# Clad: Same deal as the air gap
	# Middle
	left = clad["w"].x0
	top = clad["n"].y0
	right = clad["e"].x0
	bottom = clad["s"].y0
	clad_middle_area = (right - left)*(top - bottom)
	# Corner
	m = clad["nw"].coefficients['D']*rt
	c = sqrt(top**2 + left**2)
	h = c - m
	clad_corner_area = h**2
	
	# Now that we've got the areas of individual pieces of the puzzle,
	# get those for each of the entire regions
	fuel_area = fuel_middle_area - 4*fuel_corner_area
	gap_area = (gap_middle_area - 4*gap_corner_area) - fuel_area
	clad_area = (clad_middle_area - 4*clad_corner_area) - (fuel_area + gap_area)
	outer_area = cell_area - (fuel_area + gap_area + clad_area)
	
	if display:
		print("\tFuel area:          {0:.4} cm^2".format(fuel_area))
		print("\tGap area:           {0:.4} cm^2".format(gap_area))
		print("\tClad area:          {0:.4} cm^2".format(clad_area))
		print("\tOuter cooling area: {0:.4} cm^2".format(outer_area))
		
		# Check
		print()
		print("\tCELL AREA:", cell_area)
		print("\tDifference", cell_area - fuel_area - gap_area - clad_area - outer_area)
	
	return fuel_area, gap_area, clad_area, outer_area


def control_cell_by_material(geom, display = False):
	"""Calculate the area of each material a control lattice cell.
	The control cells are the same as the fuel cells, but with 5 concentric
	rings on the inside.
	
	Where it gets complicated is the axial zoning of the control cells,
	but fortunately, the math is easy there.
	
	Inputs:
		:param geom: instance of openmc.Geometry for the TREAT model
		:param display: Boolean; whether to print the areas for each region
						[Default: False]
	
	Outputs:
		:return: fuel_area, gap_area, clad_area, outer_area
				 Areas in cm^2 for the regions indicated
	"""
	
	fuel_area, gap_area, clad_area, outer_area = fuel_cell_by_material(geom, False)
	surfs = geom.get_all_surfaces()
	
	# Get each of the radii from the relevant rings
	control_base = 50001
	n = 5
	# Control rod
	crd_radius = surfs[50005].r
	crd_area = pi*crd_radius**2
	
	areas = [None, ]*n
	areas[0] = crd_area
	for i in range(n):
		# Note that these are in reverse: from largest to smallest,
		# due to the way they are ordered in geometry.xml
		j = n - 1 - i
		r = surfs[control_base + j].r
		areas[i] = pi*r**2 - sum(areas[0:i])
	crd_area, crd_clad_area, crd_gap_area, channel_clad_area, channel_gap_area = areas
	# And finally, account for the loss of fuel area due to the control rod channel
	fuel_area -= sum(areas)
	
	if display:
		print("\tControl rod area:   {0:.4} cm^2".format(crd_area))
		print("\tCrd clad area:      {0:.4} cm^2".format(crd_clad_area))
		print("\tInner gas gap area: {0:.4} cm^2".format(crd_gap_area))
		print("\tChannel clad area:  {0:.4} cm^2".format(channel_clad_area))
		print("\tChannel gap area:   {0:.4} cm^2".format(channel_gap_area))
		print("\tFuel area:          {0:.4} cm^2".format(fuel_area))
		print("\tGap area:           {0:.4} cm^2".format(gap_area))
		print("\tClad area:          {0:.4} cm^2".format(clad_area))
		print("\tOuter cooling area: {0:.4} cm^2".format(outer_area))
		
		# Check
		xpitch, ypitch = __setup(geom)[-2:]
		cell_area = xpitch*ypitch
		print()
		print("\tCELL AREA:", cell_area)
		print("\tDifference", cell_area - fuel_area - gap_area - clad_area - outer_area - \
		      crd_area - crd_clad_area - crd_gap_area - channel_clad_area - channel_gap_area)
	
	return crd_area, crd_clad_area, crd_gap_area, channel_clad_area, channel_gap_area, \
	       fuel_area, gap_area, clad_area, outer_area


def reflector_cell_by_material(geom, display = False):
	"""Calculate the area of each material a graphite reflector lattice cell

	Inputs:
		:param geom: instance of openmc.Geometry for the TREAT model
		:param display: Boolean; whether to print the areas for each region
						[Default: False]

	Outputs:
		:return: refl_area, gap_area, clad_area, outer_area
				 Areas in cm^2 for the materials indicated
	"""
	surfs, key_list, n, xpitch, ypitch = __setup(geom)
	
	refl_base = 90061
	refl_list = [None, ]*n
	for i in range(n):
		refl_list[i] = surfs[refl_base + i]
	refl = dict(zip(key_list, refl_list))
	
	gap_base = 90041
	gap_list = [None, ]*n
	for i in range(n):
		gap_list[i] = surfs[gap_base + i]
	gap = dict(zip(key_list, gap_list))
	
	clad_base = 90051
	clad_list = [None, ]*n
	for i in range(n):
		clad_list[i] = surfs[clad_base + i]
	clad = dict(zip(key_list, clad_list))
	
	# Graphite reflector area: Calculate the middle area and subtract the corners
	# Sides (for middle)
	right = refl["e"].x0
	left = refl["w"].x0
	top = refl["n"].y0
	bottom = refl["s"].y0
	refl_middle_area = (right - left)*(top - bottom)
	# Corner
	m = refl["nw"].coefficients['D']*rt
	c = sqrt(top**2 + left**2)
	h = c - m
	refl_corner_area = h**2
	
	# Air gap between the graphite and the cladding:
	# same shape as the reflector, but will ultimately subtract reflector
	# Middle
	left = gap["w"].x0
	top = gap["n"].y0
	right = gap["e"].x0
	bottom = gap["s"].y0
	gap_middle_area = (right - left)*(top - bottom)
	# Corner
	m = gap["nw"].coefficients['D']*rt
	c = sqrt(top**2 + left**2)
	h = c - m
	gap_corner_area = h**2
	
	# Clad: Same deal as the air gap
	# Middle
	left = clad["w"].x0
	top = clad["n"].y0
	right = clad["e"].x0
	bottom = clad["s"].y0
	clad_middle_area = (right - left)*(top - bottom)
	# Corner
	m = clad["nw"].coefficients['D']*rt
	c = sqrt(top**2 + left**2)
	h = c - m
	clad_corner_area = h**2
	
	# And then the outer gap, which is determined by the cell area
	cell_area = xpitch*ypitch
	
	# Now that we've got the areas of individual pieces of the puzzle,
	# get those for each of the entire regions
	refl_area = refl_middle_area - 4*refl_corner_area
	gap_area = (gap_middle_area - 4*gap_corner_area) - refl_area
	clad_area = (clad_middle_area - 4*clad_corner_area) - (refl_area + gap_area)
	outer_area = cell_area - (refl_area + gap_area + clad_area)
	
	if display:
		print("\tReflector area:     {0:.4} cm^2".format(refl_area))
		print("\tGap area:           {0:.4} cm^2".format(gap_area))
		print("\tClad area:          {0:.4} cm^2".format(clad_area))
		print("\tOuter cooling area: {0:.4} cm^2".format(outer_area))
		
		# Check
		print()
		print("\tCELL AREA:", cell_area)
		print("\tDifference", cell_area - refl_area - gap_area - clad_area - outer_area)
	
	return refl_area, gap_area, clad_area, outer_area


if __name__ == "__main__":
	# Test
	print("Fuel Cell:")
	fuel_cell_by_material(geometry, True)
	print("\nControl Rod Cell:")
	control_cell_by_material(geometry, True)
	print("\nReflector Cell:")
	reflector_cell_by_material(geometry, True)


