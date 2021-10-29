import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import random

def BLACKBODY(x, temp, radius):

	h = 6.63e-27	#erg s
	c = 3e10		#cm/s
	k_B = 1.38e-16	#erg/K
	x = x*1e-8
# 	R = 1e15		#cm

	y = ( 2*h*c**2/ x**5) * np.exp(-1*h*c/(x*k_B*temp)) / (4*np.pi*radius**2)

	return y

path_in = 'snex_spectra/snexdata_target4895/fluxcal_spectra/'
path_out = 'snex_spectra/snexdata_target4895/fluxcal_spectra/bb_extrapolate/'

spectra = [f for f in os.listdir(path_in) if f.endswith('.ascii')]

for spectrum in spectra:
	
	data = np.loadtxt(path_in + spectrum)
	wl = data[:,0]
	flux = data[:,1]
	
	popt, pcov = curve_fit(BLACKBODY, wl, flux, p0=[2.0e4, 1e15])

	errors = np.sqrt(np.diag(pcov))

# 	print('Spectrum:\t\t%s' %spectrum)
# 	print('Temperature:\t\t%f ± %f K' %(popt[0], errors[0]))
# 	print('Radius:\t\t\t%e ± %e cm' %(popt[1], errors[1]))
	
	wl_to_plot = np.linspace(1000, 10000, 10000)
	flux_to_plot = BLACKBODY(wl_to_plot, *popt)

	predicted_flux = BLACKBODY(wl, *popt)
	sigma_flux = predicted_flux - flux
	false_error_cap = np.std(sigma_flux)
	
	# Gauge wavelength sampling
	
	wl_sampling = np.mean(wl[1:] - wl[:-1])
# 	print(wl_sampling)
	
	wl_before = np.arange(1000, wl[0], wl_sampling)
	wl_after = np.arange(wl[-1], 15000, wl_sampling)
	
	noise_before = np.array( [random.uniform(-1*random.random()*false_error_cap, random.random()*false_error_cap) for x in range(len(wl_before))] )
	noise_after = np.array( [random.uniform(-1*random.random()*false_error_cap, random.random()*false_error_cap) for x in range(len(wl_after))] )
	
	flux_before = BLACKBODY(wl_before, *popt) + noise_before
	flux_after = BLACKBODY(wl_after, *popt) + noise_after
	
	total_wl = np.concatenate([wl_before, wl, wl_after])
	total_flux = np.concatenate([flux_before, flux, flux_after])
	
# 	print(total_wl.shape, total_flux.shape)
# 	
# 	plt.plot(total_wl, total_flux)
# 	plt.show()
	
	array_out = np.empty((len(total_wl), 2))
	array_out[:,0] = total_wl
	array_out[:,1] = total_flux
	
	np.savetxt(path_out + 'bbe_' + spectrum, array_out)
	
	
# 	plt.plot(wl, flux, color = 'blue', ls = '-', linewidth = 1, alpha = 0.6)
# 	plt.plot(wl_to_plot, flux_to_plot, color = 'red', ls = '--', linewidth = 2, alpha = 1.0)
# 	plt.fill_between(wl_to_plot, flux_to_plot - np.std(sigma_flux), flux_to_plot + np.std(sigma_flux), alpha = 0.4, color = 'red')
# 	
# 	plt.plot(wl_before, flux_before, color = 'green', ls = '-', linewidth = 1, alpha = 0.6)
# 	plt.plot(wl_after, flux_after, color = 'green', ls = '-', linewidth = 1, alpha = 0.6)
# 
# 	plt.xlabel('Wavelength, $\AA$')
# 	plt.ylabel('Flux')
# 	plt.show()
