# Treat mesh
#
# Child of the "Mesh" class in OpenMC, customized to add nuclide capabilities
# for TREAT's specific geometry

import openmc
import numpy
import area_calculator
from copy import deepcopy

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
	mesh_size : float
		Mesh size in units of complete assemblies.
		Must be given in half-integers, e.g. {0.5, 1.0, 1.5, 2.0}

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
	
	def __init__(self, mesh_id = None, name = '', geometry = None, mesh_size = None):
		super().__init__(mesh_id, name)
		self.geometry = geometry
		
		# TODO: Add the following to the docstring
		self.surfaces = self.geometry.get_all_surfaces()
		self.zactive = self.surfaces[20010].z0 - self.surfaces[20009].z0
		
		# FIXME: this breaks the code
		# self.lattice = self.geometry.get_all_lattices()[LAT_ID]
		
		self.universes = geometry.get_all_universes()
		self.materials = geometry.get_all_materials()
		self.cells = geometry.get_all_cells()
		self.fuel_cells = deepcopy(self.universes[99].cells)
		self.refl_cells = deepcopy(self.universes[26].cells)
		self.cont_cells = deepcopy(self.fuel_cells)
		for id in (50210, 50310):
			self.cont_cells[id] = self.cells[id]
		# self.nuclides = self._get_nuclides()
	
	# TODO: Convert this to private method or eliminate altogether
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
	
	def get_refl_nuclides(self):
		nuclides = self.universes[26].get_nuclides()
		return nuclides
	
	def get_cont_nuclides(self):
		nuclides = []
		for cell in self.cont_cells.values():
			for nuc in cell.get_nuclides():
				if nuc not in nuclides:
					nuclides.append(nuc)
		return nuclides
	
	def get_nuclide_densities(self, assembly_type):
		"""Return all nuclides contained in the universe
		
		Parameters
		----------
		assembly_type : str
			{"fuel", "reflector", "control"}

		Returns
		-------
		nuclide_densities : dict
			Dictionary whose keys are nuclide names and values are 3-tuples of
			(nuclide, density percent, density percent type)

		"""
		assembly_type = assembly_type.lower()
		assert assembly_type in ("fuel", "reflector", "refl", "control", "cont"), \
			'assembly_type must be in {"fuel", "reflector", "control"}'
		# nuclides = self.nuclides
		# nuclides = self.get_nuclides()
		nuclide_densities = {}
		
		if assembly_type == "fuel":
			volumes = area_calculator.fuel_cell_by_material(self.geometry)
			fuel_vol, gap_vol, clad_vol, outer_vol = [self.zactive*v for v in volumes]
			for nuc in self.get_fuel_nuclides():
				print(nuc)
				#TODO: The rest
		elif assembly_type in ("reflector", "refl"):
			volumes = area_calculator.reflector_cell_by_material(self.geometry)
			refl_vol, gap_vol, clad_vol, outer_vol = [self.zactive*v for v in volumes]
			for nuc in self.get_refl_nuclides():
				print(nuc)
			# TODO: The rest
		elif assembly_type in ("control", "cont"):
			volumes = area_calculator.control_cell_by_material(self.geometry)
			crd_vol, crd_clad_vol, crd_gap_vol, channel_clad_vol, channel_gap_vol, \
			fuel_vol, gap_vol, clad_vol, outer_vol = [self.zactive*v for v in volumes]
			for nuc in self.get_cont_nuclides():
				print(nuc)
			# TODO: the rest
			
		raise NotImplementedError("Treat_Mesh.get_nuclide_densities() has not been implemented yet.")


# test
if __name__ == "__main__":
	summ = openmc.Summary("summary.h5")
	geom = summ.geometry
	mesh = Treat_Mesh(geometry = geom)
	# mesh.get_nuclides()
	mesh.get_nuclide_densities(assembly_type = "fuel")
