# Treat mesh
#
# Child of the "Mesh" class in OpenMC, customized to add nuclide capabilities
# for TREAT's specific geometry

import openmc
import numpy
import area_calculator

LAT_ID = 100
FUEL_UNIVERSE = 9  # 99 for active fuel
REFL_UNIVERSE = 4  # 26 for active reflector
CONTROL_UNIVERSE = 5  # 99 (minus hole) for active crd


class Treat_Mesh(openmc.Mesh):
	"""A structured Cartesian mesh in one, two, or three dimensions

	Parameters
	----------
	mesh_id : int
		Unique identifier for the mesh
	name : str
		Name of the mesh
	geometry: openmc.Geometry
		Geometry of the TREAT model

	Attributes
	----------
	id : int
		Unique identifier for the mesh
	name : str
		Name of the mesh
	type : str
		Type of the mesh
	dimension : Iterable of int
		The number of mesh cells in each direction.
	lower_left : Iterable of float
		The lower-left corner of the structured mesh. If only two coordinate are
		given, it is assumed that the mesh is an x-y mesh.
	upper_right : Iterable of float
		The upper-right corner of the structured mesh. If only two coordinate
		are given, it is assumed that the mesh is an x-y mesh.
	width : Iterable of float
		The width of mesh cells in each direction.

	"""
	
	def __init__(self, mesh_id = None, name = '', geometry = None):
		super().__init__(mesh_id, name)
		self.geometry = geometry
		
		# TODO: Add the following to the docstring
		self.surfaces = self.geometry.get_all_surfaces()
		self.zactive = self.surfaces[20010].z0 - self.surfaces[20009].z0
		
		# FIXME: this breaks the code
		# self.lattice = self.geometry.get_all_lattices()[LAT_ID]
		
		self.universes = geometry.get_all_universes()
		self.cells = geometry.get_all_cells()
		self.fuel_cells = self.universes[99].cells
		self.refl_cells = self.universes[26].cells
		self.cont_cells = self.fuel_cells
		for id in (50110, 50210, 50310):
			self.cont_cells[id] = self.cells[id]
		
		# self.nuclides = self._get_nuclides()
	
	# Convert this to private method
	def get_nuclides(self):
		"""Not implemented yet"""
		# TODO: Determine whether to return nuclides for an assembly type or all nuclides
		nuclides = []
		# FIXME: RecursionError
		for cell in tuple(self.cont_cells.values()) + tuple(self.refl_cells.values()):
			# self.cont_cells ontains self.fuel_cells
			for nuclide in cell.get_nuclides():
				if nuclide not in nuclides:
					nuclides.append(nuclide)
		
		return nuclides
	
	# Test: see if the recursion error occurs here too
	# FIXME: RecursionError does indeed occur
	def get_fuel_nuclides(self):
		nuclides = self.universes[99].get_nuclides()
		return nuclides
	
	def get_nuclide_densities(self, assembly):
		raise NotImplementedError("Treat_Mesh.get_nuclide_densities() has not been implemented yet.")


# test
if __name__ == "__main__":
	summ = openmc.Summary("summary.h5")
	geom = summ.geometry
	mesh = Treat_Mesh(geometry = geom)
	# mesh.get_nuclides()
	mesh.get_nuclides()
