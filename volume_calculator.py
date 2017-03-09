# Volume calculator for TREAT lattice cells


import openmc
from math import sqrt

# Global variables
# Extract the geometry from an existing summary
summary = openmc.Summary("summary.h5")
geometry = summary.openmc_geometry


def fuel_cell_by_material(summ):
	"""Calculate the volume of each material a fuel lattice cell
	
	:param summ: instance of openmc.Summary for the TREAT model
	
	:return: fuel_vol, gap_vol, clad_vol, outer_vol
	Volumes in cm^3 for
	"""
	# The main core lattice which the cells appear in
	core_lat = summ.get_lattice_by_id(100)
	[xpitch, ypitch, z] = core_lat.pitch
	# TODO: Ascertain whether `z` is the desired height.
	
	# Useful shorthand
	rt = sqrt(2) / 2    # "root two"
	
	# Create dictionaries for each region
	key_list = ("e", "s", "w", "n", "se", "sw", "nw", "ne")
	n = len(key_list)
	
	fuel_base = 90001
	fuel_list = [None,]*n
	for i in range(n):
		fuel_list[i] = summ.get_surface_by_id(fuel_base + i)
	fuel = dict(zip(key_list, fuel_list))
	
	gap_base = 90011
	gap_list = [None, ] * n
	for i in range(n):
		gap_list[i] = summ.get_surface_by_id(gap_base + i)
	gap = dict(zip(key_list, gap_list))
	
	clad_base = 90021
	clad_list = [None, ] * n
	for i in range(n):
		clad_list[i] = summ.get_surface_by_id(clad_base + i)
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
	cell_area = xpitch * ypitch
	
	# Fuel area: Calculate the middle area and subtract the corners
	# Sides (for middle)
	right = fuel["e"].x0
	left = fuel["w"].x0
	top = fuel["n"].y0
	bottom = fuel["s"].y0
	fuel_middle_area = (right - left) * (top - bottom)
	# Corner
	m = fuel["nw"].coefficients['D'] * rt
	c = sqrt(top ** 2 + left ** 2)
	h = c - m
	fuel_corner_area = h ** 2
	
	# Air gap calculation: same as fuel, but will ultimately subtract fuel
	# Middle
	left = gap["w"].x0
	top = gap["n"].y0
	right = gap["e"].x0
	bottom = gap["s"].y0
	gap_middle_area = (right - left) * (top - bottom)
	# Corner
	m = gap["nw"].coefficients['D'] * rt
	c = sqrt(top ** 2 + left ** 2)
	h = c - m
	gap_corner_area = h ** 2
		
	# Clad: Same deal as the air gap
	# Middle
	left = clad["w"].x0
	top = clad["n"].y0
	right = clad["e"].x0
	bottom = clad["s"].y0
	clad_middle_area = (right - left) * (top - bottom)
	# Corner
	m = clad["nw"].coefficients['D'] * rt
	c = sqrt(top ** 2 + left ** 2)
	h = c - m
	clad_corner_area = h ** 2
	
	# Now that we've got the areas of individual pieces of the puzzle,
	# get those for each of the entire regions
	fuel_area = fuel_middle_area - 4*fuel_corner_area
	gap_area = (gap_middle_area - 4*gap_corner_area) - fuel_area
	clad_area = (clad_middle_area - 4*clad_corner_area) - (fuel_area + gap_area)
	outer_area = cell_area - (fuel_area + gap_area + clad_area)
	
	print("Fuel area: {0:.4} cm^2".format(fuel_area))
	print("Gap area: {0:.4} cm^2".format(gap_area))
	print("Clad area: {0:.4} cm^2".format(clad_area))
	print("Outer cooling area: {0:.4} cm^2".format(outer_area))
	
	# Check
	print("\n")
	print("CELL AREA:", cell_area)
	print("Difference", cell_area - fuel_area - gap_area - clad_area - outer_area)
	
	return fuel_area*z, gap_area*z, clad_area*z, outer_area*z


# TODO: Add getter for reflector cell and control cell


if __name__ == "__main__":
	# Test
	fuel_vol, gap_vol, clad_vol, outer_vol = fuel_cell_by_material(summary)
	


