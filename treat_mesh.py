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

MAT_IDS = {"fuel"    : 90000,
           "air"     : 20000,
           "zirc"    : 20005,
           "aluminum": 20008,
           "graphite": 20012}


def merge_nuclide_densities(old_dict, new_dict, vfrac):
	"""Add a dictionary of nuclide densities to an existing dictionary
	by volume fraction.
	
	Parameters
	----------
	old_dict : dictionary
		Existing Dictionary whose keys are nuclide names and values are
			3-tuples of (nuclide, density percent, density percent type)
	
	new_dict : dictionary
		Dictionary to merge into old_dict. Must be of the same format,
			and nuclides must share fraction type.
	
	vfrac : float
		Volume fraction of the material whose nuclides are in new_dict
	
	
	Returns
	-------
	old_dict : dictionary
		The original dictionary updated with the values from new_dict
	
	"""
	for key in new_dict:
		if key in old_dict:
			old_tuple = old_dict[key]
			new_tuple = new_dict[key]
			# Third entry is the percent type
			old_type = old_tuple[2]
			new_type = new_tuple[2]
			errstr = "Density percents must be of the same type. \
			Expected '{}', got '{}'".format(old_type, new_type)
			assert old_type == new_type, errstr
			
			# Update the dictionary with the appropriate fraction of this nuclide
			merged_frac = old_tuple[1] + new_tuple[1]*vfrac
			old_dict[key] = (old_tuple[0], merged_frac, old_type)
		else:
			old_dict[key] = new_dict[key]
	return old_dict


def merge_nuclide_densities_by_cell(cell_dict, vfrac_list, nuclide_densities = None):
	"""Find and merge the nuclide densities for several OpenMC cells

	Parameters
	----------
	cell_dict : collections.OrderedDict
		Ordered dictionary of {cell_id : openmc.Cell} from which each
			cell nuclide densities will be looked up
	
	vfrac_list : list or tuple
		List of the volume fractions corresponding to each Cell.
			Must be the same length as cell_dict.
	
	nuclide_densities : dictionary, optional
		Existing Dictionary whose keys are nuclide names and values are
			3-tuples of (nuclide, density percent, density percent type).
			Nuclide densities from cell_dict will be merged into this.
	
	Returns
	-------
	nuclide_densities: dictionary

	"""
	n = len(cell_dict)
	assert n == len(vfrac_list), \
		"Number of volume fractions given does not equal number of cells."
	
	if nuclide_densities is None:
		nuclide_densities = {}
	i = 0
	
	for id in cell_dict:
		cell_nuc_dens = cell_dict[id].get_nuclide_densities()
		v = vfrac_list[i]
		i += 1
		nuclide_densities = merge_nuclide_densities(nuclide_densities, cell_nuc_dens, v)
	return nuclide_densities


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
		self._surfaces = self.geometry.get_all_surfaces()
		self.zactive = self._surfaces[20010].z0 - self._surfaces[20009].z0
		self._universes = geometry.get_all_universes()
		self._materials = geometry.get_all_materials()
		self._cells = geometry.get_all_cells()
		self.fuel_cells = deepcopy(self._universes[99].cells)
		self.refl_cells = deepcopy(self._universes[26].cells)
		self.cont_cells = deepcopy(self.fuel_cells)
		for id in (50210, 50310):
			self.cont_cells[id] = self._cells[id]
	
	def get_nuclides(self):
		"""Return all of the nuclides in the active region of the core.
		
		Returns
		-------
		nuclides : list of strings
			Names of all the nuclides in the model
		
		"""
		nuclides = []
		for cell in tuple(self.cont_cells.values()) + tuple(self.refl_cells.values()):
			for nuclide in cell.get_nuclides():
				if nuclide not in nuclides:
					nuclides.append(nuclide)
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
		nuclide_densities = {}
		
		if assembly_type == "fuel":
			areas = area_calculator.fuel_cell_by_material(self.geometry)
		elif assembly_type in ("reflector", "refl"):
			areas = area_calculator.reflector_cell_by_material(self.geometry)
		elif assembly_type in ("control", "cont"):
			areas = area_calculator.reflector_cell_by_material(self.geometry)
		else:
			raise NotImplementedError("Unknown assembly type: " + assembly_type)
		
		volumes = [self.zactive*a for a in areas]
		total_volume = sum(volumes)
		vfracs = [v/total_volume for v in volumes]
		nuclide_densities = \
			merge_nuclide_densities_by_cell(self.fuel_cells, vfracs, nuclide_densities)
		return nuclide_densities


# test
if __name__ == "__main__":
	summ = openmc.Summary("summary.h5")
	geom = summ.geometry
	mesh = Treat_Mesh(geometry = geom)
	mesh.get_nuclides()
	fuel_nuc_dens = mesh.get_nuclide_densities(assembly_type = "fuel")
	refl_nuc_dens = mesh.get_nuclide_densities(assembly_type = "refl")
	cont_nuc_dens = mesh.get_nuclide_densities(assembly_type = "cont")
