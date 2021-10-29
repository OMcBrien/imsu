import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sncosmo

def addTimeSeriesSource(obj_name = 'at2018cow', explosion_epoch = 2458284.8, redshift = 0.014, min_wl = 2000., max_wl = 12000., npoints = 3000):

	import spectres

# 	path_to_spectra = 'custom_sources/%s/source_spectra/' %obj_name
	path_to_spectra = 'snex_spectra/snexdata_target4895/fluxcal_spectra/bb_extrapolate/'

# 	df = pd.read_csv(path_to_spectra + '%s_metadata.csv' %obj_name)
	df = pd.read_csv(path_to_spectra + '%s_metadata.csv' %obj_name)

	new_wl = np.linspace(min_wl, max_wl, npoints)

	fluxes = [np.zeros_like(new_wl)]
	obj_phase = [0.0]


	for index, row in df.iterrows():
	
		try:
			spectrum = np.loadtxt(path_to_spectra + 'bbe_fluxcal_' + row['Ascii file'])
		except:
			print(row['Ascii file'])
		
		wl = spectrum[:,0]
		flux = spectrum[:,1]
		
# 		wl, flux = zip( *sorted( zip( wl, flux ) ) )

		print('### ', np.nanmin(wl), np.nanmax(wl))
		
# 		if  np.nanmin(wl) < min_wl or max_wl > np.nanmax(wl):
# 			print('Spectral range out of bound. Skipping.')
# 			continue
	
		try:
	
			new_flux = spectres.spectres(new_wl, wl, flux, spec_errs=None, fill = 0.0, verbose=True)
			obj_phase.append(row['JD'] - explosion_epoch)
			fluxes.append(new_flux)
			plt.plot(new_wl, new_flux, label = '%.3f' %(row['JD'] - explosion_epoch))
		except:
			
			continue

	plt.legend()
	plt.show()

	fluxes = np.array(fluxes, dtype = float)
	obj_phase = np.array(obj_phase, dtype = float)

	source = sncosmo.TimeSeriesSource(phase = obj_phase, wave = new_wl, flux = fluxes, name = 'fbot_2018cow')

	print(source)

# 	model = sncosmo.Model(source)
# 	model.update({'z': redshift, 't0': obj_phase[0]})
# 	
# 	model.set_source_peakabsmag(-20.0, 'sdssr', 'ab')
# 
# 	tobs = np.linspace(1., 45, 60)
# 	# mags = model.bandmag('desr', 'ab', tobs)
# 	mags = model.bandmag('sdssr', 'ab', tobs)
# 	# mags = model.flux(5, [4000., 4100., 4200.])
# 
# 	plt.plot(tobs, mags)
# 	plt.gca().invert_yaxis()
# 
# 	plt.show()
	
	return source, obj_phase


df = pd.read_csv('photometric_log.csv')


the_source, obj_phase = addTimeSeriesSource()
model = sncosmo.Model(source = the_source)
model.update({'z': 0.014, 't0': 58284.3})


model.set_source_peakabsmag(-19.5, 'sdssr', 'ab')

g_mags = model.bandmag('sdssg', 'ab', df['MJD'])
r_mags = model.bandmag('sdssr', 'ab', df['MJD'])
i_mags = model.bandmag('sdssi', 'ab', df['MJD'])
z_mags = model.bandmag('sdssz', 'ab', df['MJD'])

plt.plot(df['MJD'], df['g'], ls = 'None', marker = 'o', ms = 8, mfc = 'green', mec = 'black', alpha = 0.6)
plt.plot(df['MJD'], df['r'] + 1, ls = 'None', marker = 'o', ms = 8, mfc = 'firebrick', mec = 'black', alpha = 0.6)
plt.plot(df['MJD'], df['i'] + 2, ls = 'None', marker = 'o', ms = 8, mfc = 'goldenrod', mec = 'black', alpha = 0.6)
plt.plot(df['MJD'], df['z'] + 3, ls = 'None', marker = 'o', ms = 8, mfc = 'saddlebrown', mec = 'black', alpha = 0.6)

plt.plot(df['MJD'], g_mags, ls = '--', linewidth = 2, marker = 'None', color = 'green')
plt.plot(df['MJD'], r_mags + 1, ls = '--', linewidth = 2, marker = 'None', color = 'firebrick')
plt.plot(df['MJD'], i_mags + 2, ls = '--', linewidth = 2, marker = 'None', color = 'goldenrod')
plt.plot(df['MJD'], z_mags + 3, ls = '--', linewidth = 2, marker = 'None', color = 'saddlebrown')

for phase in obj_phase:
	plt.axvline(phase + 58284.3, ls = '-', color = 'grey')

plt.xlabel('MJD')
plt.ylabel('Magnitude')

plt.gca().invert_yaxis()
plt.show()

