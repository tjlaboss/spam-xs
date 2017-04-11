#############
# Export to MOC
import openmoc
import openmoc.materialize
import openmc.openmoc_compatible

sp = openmc.StatePoint('statepoint.020.h5')
mesh_lib.load_from_statepoint(sp)

moc_geom = openmc.openmoc_compatible.get_openmoc_geometry(sp.summary.geometry)
materials = openmoc.materialize.load_openmc_mgxs_lib(mesh_lib, moc_geom)
# Generate tracks for OpenMOC
# note: increase num_azim and decrease azim_spacing for actual results (as for TREAT)
track_generator = openmoc.TrackGenerator(moc_geom, num_azim = 32, azim_spacing = 0.1)
track_generator.generateTracks()

# Run OpenMOC
solver = openmoc.CPUSolver(track_generator)
solver.computeEigenvalue()
