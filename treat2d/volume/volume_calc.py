# Volume Calc
#
# Test the volume calculator and compare it to the analytical

import openmc

RUN = True


summ = openmc.Summary("../summary.h5")
geom = summ.geometry

# Extract cells to calculate
lat = geom.get_all_lattices()[100]
universes = geom.get_all_universes()
cells = geom.get_all_cells()
fuel_verse1 = universes[9]
fuel_verse2 = universes[98]
fuel_cell = cells[90011]

if RUN:
	# The calculation itself
	px, py = lat.pitch[0:2]
	ll = (-px/2, -py/2, -0.5)
	ur = (px/2, py/2, +0.5)
	verse_calc = openmc.VolumeCalculation([fuel_verse1, fuel_verse2], int(1E4), ll, ur)
	cell_calc = openmc.VolumeCalculation([fuel_cell], int(1E5), ll, ur)
	
	# Settings file
	settings = openmc.Settings()
	settings.volume_calculations = [verse_calc, cell_calc]
	settings.export_to_xml()
	openmc.calculate_volumes(threads = 12, openmc_exec = "sterlingmc")

# Extract results

