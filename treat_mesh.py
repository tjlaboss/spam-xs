# Treat mesh
#
# Child of the "Mesh" class in OpenMC, customized to add nuclide capabilities
# for TREAT's specific geometry

import openmc
import volume_calculator

class Treat_Mesh(openmc.Mesh):
	"""A structured Cartesian mesh in one, two, or three dimensions

	Parameters
	----------
	mesh_id : int
		Unique identifier for the mesh
	name : str
		Name of the mesh

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
		The upper-right corner of the structrued mesh. If only two coordinate
		are given, it is assumed that the mesh is an x-y mesh.
	width : Iterable of float
		The width of mesh cells in each direction.

	"""
	
	
	def get_nuclides(self):
		"""Not implemented yet"""
		raise NotImplementedError("Treat_Mesh.get_nuclides() has not been implemented yet.")