import sys
import time
import sncosmo
import afterglowpy as grb
import sfdmap
import json
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import population_generator
import extinction

def mJy2ABmag(F_nu, wl):

	F_lam = 3.00e-5*F_nu*1e-3/(wl**2)
	
	AB_mag = -2.5*np.log10(F_lam) - 48.6 - 2.5*np.log10((wl**2)/3e18)
	
	return AB_mag

def registerBandpasses(lightcurve_settings):

	for flt in lightcurve_settings['telescope properties']['filter effective wavelengths'].keys():
		
		print('Registering ', flt)

		flt_data = np.loadtxt('bandpasses/%s' %str(flt))
		flt_band = sncosmo.Bandpass(flt_data[:,0], flt_data[:,1], name = str(flt))
		sncosmo.register(flt_band)


def addTimeSeriesSource(obj_name = 'at2017gfo', explosion_epoch = 2457983.48, redshift = 0.00984, min_wl = 3400., max_wl = 22300., npoints = 40000):

	import spectres

	path_to_spectra = 'custom_sources/%s/source_spectra/' %obj_name

	df = pd.read_csv(path_to_spectra + '%s_metadata.csv' %obj_name)

	new_wl = np.linspace(min_wl, max_wl, npoints)

	fluxes = [np.zeros_like(new_wl)]
	obj_phase = [0.0]


	for index, row in df.iterrows():
	
		try:
			spectrum = np.loadtxt(path_to_spectra + row['Ascii file'])
		except:
			print(row['Ascii file'])
		
		wl = spectrum[:,0]
		flux = spectrum[:,1]
	
		new_flux = spectres.spectres(new_wl, wl, flux, spec_errs=None, fill=None, verbose=True)
	
		obj_phase.append(row['JD'] - explosion_epoch)
	
		fluxes.append(new_flux)
	
# 		plt.plot(new_wl, new_flux, label = '%.3f' %(row['JD'] - explosion_epoch))
# 
# 	plt.legend()
# 	plt.show()

	fluxes = np.array(fluxes, dtype = float)
	obj_phase = np.array(obj_phase, dtype = float)

	source = sncosmo.TimeSeriesSource(phase = obj_phase, wave = new_wl, flux = fluxes, name = 'custom_source')

# 	print(source)
# 
# 	model = sncosmo.Model(source)
# 	model.update({'z': redshift, 't0': obj_phase[0]})
# 	
# 	# model.set_source_peakabsmag(-15.6, 'desr', 'ab')
# 
# 	tobs = np.linspace(0., 45., 50)
# 	# mags = model.bandmag('desr', 'ab', tobs)
# 	mags = model.bandmag('desr', 'ab', tobs)
# 	# mags = model.flux(5, [4000., 4100., 4200.])
# 
# 	print(tobs, mags)
# 
# 	plt.plot(tobs, mags)
# 	plt.gca().invert_yaxis()
# 
# 	plt.show()
	
	return source


def collectTemporalOverlap(query_output, explosion_epoch, duration, column_headings):
	
	time_heading = column_headings['time']

	partial_query_output = query_output.query('{} > {} & {} < {}'.format(time_heading, explosion_epoch, time_heading, explosion_epoch + duration))

# 	print('### partial query, post-temp. filtering')
# 	print(partial_query_output[['jd', 'ra_c', 'dec_c']])

	return partial_query_output

def collectFootprintOverlap(query_output, ra_trans, dec_trans, lightcurve_settings):
	
	ra_heading = lightcurve_settings['telescope properties']['column headings']['ra']
	dec_heading = lightcurve_settings['telescope properties']['column headings']['dec']

	chip_width = lightcurve_settings['telescope properties']['chip width'] / 2.
	chip_height = lightcurve_settings['telescope properties']['chip height'] / 2.
	
	max_allowed_ra = 360.
	min_allowed_ra = 0.

	upper_transient_ra_bound = ra_trans + (chip_width * np.cos(dec_trans*np.pi/180.)) / 2.
	lower_transient_ra_bound = ra_trans - (chip_width * np.cos(dec_trans*np.pi/180.)) / 2.
	
	upper_transient_dec_bound = dec_trans + chip_height / 2.
	lower_transient_dec_bound = dec_trans - chip_height / 2.
	
	# Filter RA first as it wraps (360 --> 0 degrees)
	if upper_transient_ra_bound > max_allowed_ra:
	
		partial_query_output = query_output.query('%s >= %f | %s <= %f' %(ra_heading, lower_transient_ra_bound, ra_heading, (upper_transient_ra_bound - max_allowed_ra)) )
	
	elif lower_transient_ra_bound < min_allowed_ra:
	
		partial_query_output = query_output.query('%s >= %f | %s <= %f' %(ra_heading, (lower_transient_ra_bound + max_allowed_ra), ra_heading, (upper_transient_ra_bound - max_allowed_ra)) )
	
	else:

		partial_query_output = query_output.query('%s >= %f & %s <= %f' %(ra_heading, lower_transient_ra_bound, ra_heading, upper_transient_ra_bound) )
	
	# Now filter by DEC
	partial_query_output = partial_query_output.query('%s >= %f & %s <= %f' %(dec_heading, lower_transient_dec_bound, dec_heading, upper_transient_dec_bound) )

	return partial_query_output

def collectFilterOverlap(query_output, column_headings, flt):
	
	partial_query_output = query_output.query('%s == "%s"' %(column_headings['filter'], flt))

	return partial_query_output


def generateSNLightcurvePopulation(query_output, population_settings, lightcurve_settings, io_settings, source):

	tlc0 = time.time()

	transient_type = lightcurve_settings['population']['transient type']
	
	with open('sources.yaml', 'r') as stream:
		loader = yaml.SafeLoader
		sources = yaml.load(stream, Loader = loader)
	
	source_info = sources[transient_type]
	
	filters_in_use = list(lightcurve_settings['telescope properties']['filter effective wavelengths'].keys())
	filter_keys = lightcurve_settings['telescope properties']['column headings']['filter keys']
	
	dict_flt_keys = dict( zip( filters_in_use, filter_keys ) )
	
	column_headings = lightcurve_settings['telescope properties']['column headings']

	dust = sncosmo.CCM89Dust()
	dustmap = sfdmap.SFDMap(io_settings['sfdmaps path'])
	
	if source:

		model = sncosmo.Model(source = source, effects = [dust, dust], effect_names = ['host', 'mw'], effect_frames = ['rest', 'obs'])
	
	else:
	
		model = sncosmo.Model(source = source_info['source'], effects = [dust, dust], effect_names = ['host', 'mw'], effect_frames = ['rest', 'obs'])
	
	main_lc_dict = {}

	for index, row in population_settings.iterrows():
	
		transient_index_str = row['transient_id']
		t0 = row['explosion epochs']
		distance = row['distances']
		z = row['redshifts']
		ra = row['ras']
		dec = row['decs']

		ebv = dustmap.ebv(ra, dec)

		model.update({'z': z, 't0': t0, 'mwebv': ebv})
		
		try:
		
			model.set_source_peakabsmag(source_info['peak magnitude'], filters_in_use[0], 'ab')
			
		except:
			
			registerBandpasses(lightcurve_settings)
			model.set_source_peakabsmag(source_info['peak magnitude'], filters_in_use[0], 'ab')
		
		partial_query_output = collectTemporalOverlap(query_output, t0, source_info['duration'], column_headings)
		partial_query_output = collectFootprintOverlap(partial_query_output, ra, dec, lightcurve_settings)

		transient_dict = {'explosion epoch': t0,
						  'distance': distance,
						  'redshift': z,
						  'ra': ra,
						  'dec': dec,
						  'ebv': ebv,
						  'lightcurves': {}}
		
		flt_dict = {}
		
		for flt in filters_in_use:
		
			flt_partial_query_output = collectFilterOverlap(partial_query_output, column_headings, dict_flt_keys[flt]).sort_values(lightcurve_settings['telescope properties']['column headings']['time'])
		
			tobs = np.array(flt_partial_query_output[column_headings['time']])
			mag5sig = np.array(flt_partial_query_output[column_headings['limiting mag']])
			
			mags = model.bandmag(flt, 'ab', tobs)
			
			ind = np.where(mags < mag5sig)
			tobs_rec = tobs[ind]
			mags_rec = mags[ind]
			
			lc_dict = {'tobs': list(tobs),
					   'mags': list(mags),
					   'limmag5': list(mag5sig),
					   'tobs_rec': list(tobs_rec),
					   'mags_rec': list(mags_rec)}
			
			flt_dict.update({flt: lc_dict})

		transient_dict['lightcurves'].update(flt_dict)
		
		main_lc_dict.update({transient_index_str: transient_dict})
	
	tlc1 = time.time()
	print('Lightcurves generated in %.4f s' %(tlc1 - tlc0))

	return main_lc_dict

def generateGRBLightcurvePopulation(query_output, population_settings, lightcurve_settings, grb_settings, io_settings):
	
	jet_type_dict = {'TopHat': grb.jet.TopHat,
					 'Gaussian': grb.jet.Gaussian,
					 'PowerLaw': grb.jet.PowerLaw}
	
	filters_in_use = list(lightcurve_settings['telescope properties']['filter effective wavelengths'].keys())
	filter_keys = lightcurve_settings['telescope properties']['column headings']['filter keys']
	
	dict_flt_keys = dict( zip( filters_in_use, filter_keys ) )
	
	column_headings = lightcurve_settings['telescope properties']['column headings']
	
	effective_wavelengths = lightcurve_settings['telescope properties']['filter effective wavelengths']
				   
	c = 3.e18
	seconds_in_day = 8.64e4
	step_sizen = 25
	
	goto_nueffs = {'gotoL': c / 5396.65,
				   'gotoR': c / 6405.05,
				   'gotoG': c / 5355.32, 
				   'gotoB': c / 4600.58}
	
	effective_frequencies = {k: c / v for k, v in effective_wavelengths.items()}
	
	grb_settings['jetType'] = jet_type_dict[grb_settings['jetType']]
	
	main_lc_dict = {}
	
	dustmap = sfdmap.SFDMap(io_settings['sfdmaps path'])

	for index, row in population_settings.iterrows():
	
		transient_index_str = row['transient_id']
		t0 = row['explosion epochs']
		distance = row['distances']
		z = row['redshifts']
		ra = row['ras']
		dec = row['decs']

		grb_settings['z'] = z
		grb_settings['d_L'] = population_generator.redshift2distance(z)['dl_mpc'] * 3.086e24
		
		ebv = dustmap.ebv(ra, dec)
		
		partial_query_output = collectTemporalOverlap(query_output, t0, 350., column_headings)
		partial_query_output = collectFootprintOverlap(partial_query_output, ra, dec, lightcurve_settings)

		transient_dict = {'explosion epoch': t0,
						  'distance': distance,
						  'redshift': z,
						  'ra': ra,
						  'dec': dec,
						  'ebv': ebv,
						  'lightcurves': {}}
		
		flt_dict = {}
		
		for flt in filters_in_use:
		
			flt_partial_query_output = collectFilterOverlap(partial_query_output, column_headings, dict_flt_keys[flt]).sort_values(lightcurve_settings['telescope properties']['column headings']['time'])
		
			tobs = (np.array(flt_partial_query_output[column_headings['time']]) - t0) * seconds_in_day
			mag5sig = np.array(flt_partial_query_output[column_headings['limiting mag']])
			
			if len(tobs) == 0:
			
				lc_dict = {'tobs': [],
				    	   'mags': [],
					       'limmag5': [],
					       'tobs_rec': [],
					   	   'mags_rec': []}
			
				flt_dict.update({flt: lc_dict})

				continue
			
			else:
			
				bandpass = np.loadtxt('bandpasses/%s' %flt)
				wl = bandpass[:,0]
			
				flux_mJy_total = np.zeros_like(tobs)
				nu_single = np.full(len(tobs), c / wl[0])

				flux_old = grb.fluxDensity(tobs, nu_single, **grb_settings)

				for i in range(step_sizen + 1, len(wl), step_sizen):
					
# 					print(i, wl[i])
	
					nu_single[:] = c / wl[i]
	
					flux_new = grb.fluxDensity(tobs, nu_single, **grb_settings)
	
					height = wl[i] - wl[i-step_sizen]
	
					flux_mJy_total += 0.5*(flux_old + flux_new)*height
	
					flux_old = flux_new
			
# 				nu = np.full(len(tobs), c / effective_wavelengths[flt])
# 				flux_mJy = grb.fluxDensity(tobs, nu, **grb_settings)
# 				mags = mJy2ABmag(flux_mJy, effective_wavelengths[flt])
				mags = mJy2ABmag(flux_mJy_total, effective_wavelengths[flt])
				
				A_ex = extinction.fitzpatrick99(np.array([effective_wavelengths[flt]]), 3.1 * ebv)
# 				print(A_ex)

				ind = np.where(mags < mag5sig)
				tobs_rec = tobs[ind]
				mags_rec = mags[ind]
			
				lc_dict = {'tobs': list(tobs / seconds_in_day),
						   'mags': list(mags + A_ex[0]),
						   'limmag5': list(mag5sig),
						   'tobs_rec': list(tobs_rec / seconds_in_day),
						   'mags_rec': list(mags_rec + A_ex[0])}
			
				flt_dict.update({flt: lc_dict})


		transient_dict['lightcurves'].update(flt_dict)
		
		main_lc_dict.update({transient_index_str: transient_dict})
	
	return main_lc_dict

def saveLightcurveData(lightcurve_path, lc_population):
	
	with open(lightcurve_path, 'w') as f:
		json.dump(lc_population, f, indent = 4)

	