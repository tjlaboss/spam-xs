import pylab


def plot_reaction_rates(eps = 1E-6):
	moc_fission_rates = pylab.loadtxt("moc_data/moc_fission_rates")
	fission_rates = pylab.loadtxt("moc_data/montecarlo_fission_rates")
	errors = pylab.divide(moc_fission_rates - fission_rates, fission_rates/100.0)
	
	# Filter out values that are essentially zero
	# Objective: Ignore zero fission rates in guide tubes with Matplotlib color scheme
	indices = fission_rates <= eps
	fission_rates[indices] = pylab.NaN
	moc_fission_rates[indices] = pylab.NaN
	errors[indices] = pylab.NaN
	
	# Plot OpenMC's fission rates in the left subplot
	fig = pylab.subplot(221)
	pylab.imshow(fission_rates.squeeze(), interpolation = 'none', cmap = 'jet')
	pylab.title('OpenMC Fission Rates')
	
	# Plot OpenMOC's fission rates in the right subplot
	fig2 = pylab.subplot(222)
	pylab.imshow(moc_fission_rates.squeeze(), interpolation = 'none', cmap = 'jet')
	pylab.title('OpenMOC Fission Rates')
	
	fig3 = pylab.subplot(223)
	pct = pylab.imshow(errors.squeeze(), interpolation = 'none', cmap = 'jet')
	pylab.clim(-100, 100)
	pylab.title('Percent error')
	pylab.colorbar(pct)
	
	pylab.tight_layout()
	pylab.show()
	

if __name__ == "__main__":
	plot_reaction_rates()
