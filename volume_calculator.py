# Volume calculator for TREAT lattice cells


import openmc
#from openmc import mgxs as mgxs
from math import sqrt

# Global variables
# Extract the geometry from an existing summary
summary = openmc.Summary("summary.h5")
geometry = summary.openmc_geometry


def fuel_cell_materials(summ):
	"""Calculate the volume of each material a fuel lattice cell
	
	:param geom: instance of openmc.Geometry for the TREAT model
	:param TBD:
	
	:return:
	"""
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

	'''
	fuel_base = 90001
	fuel_list = [None, ] * n
	for i in range(n):
		fuel_list[i] = summ.get_surface_by_id(fuel_base + i)
	fuel = dict(zip(key_list, fuel_list))
	'''
	
	
	# Get the surface area (which is all this code does now)
	
	# Fuel area: Calculate the middle area and subtract the corners
	# Sides (for middle)
	right = fuel["e"].x0
	left = fuel["w"].x0
	top = fuel["n"].y0
	bottom = fuel["s"].y0
	middle_area = (right - left)*(top - bottom)
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
	m = fuel["nw"].coefficients['D'] * rt
	c = sqrt(top**2 + left**2)
	h = c - m
	corner_area = h**2
	# By symmetry
	fuel_area = middle_area - 4*corner_area
	print("Fuel area: {0:.4} cm^2".format(fuel_area))
	
	
	
	
	


if __name__ == "__main__":
	# Test
	fuel_cell_materials(summary)


