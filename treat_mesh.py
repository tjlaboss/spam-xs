# Treat mesh
#
# Child of the "Mesh" class in OpenMC, customized to add nuclide capabilities
# for TREAT's specific geometry

import openmc
import area_calculator

LAT_ID = 100
FUEL_UNIVERSE = 9       # 99 for active fuel
REFL_UNIVERSE = 4       # 26 for active reflector
CONTROL_UNIVERSE = 5    # 99 (minus hole) for active crd


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
		The upper-right corner of the structrued mesh. If only two coordinate
		are given, it is assumed that the mesh is an x-y mesh.
	width : Iterable of float
		The width of mesh cells in each direction.

	"""
	
	def __init__(self, mesh_id = None, name = '', geometry = None):
		super().__init__(mesh_id, name)
		self.geometry = geometry
	
	
	def get_nuclides(self):
		"""Not implemented yet"""
		raise NotImplementedError("Treat_Mesh.get_nuclides() has not been implemented yet.")
	
	

# test
if __name__ == "__main__":
	summ = openmc.Summary("summary.h5")
	geom = summ.geometry
	mesh = Treat_Mesh(geometry = geom)
	mesh.get_nuclides()



