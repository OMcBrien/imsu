import sys
import time
import sncosmo
import afterglowpy as grb
import sfdmap
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import population_generator
import extinction

def mJy2ABmag(F_nu, wl):

	F_lam = 3.00e-5*F_nu*1e-3/(wl**2)
	
	AB_mag = -2.5*np.log10(F_lam) - 48.6 - 2.5*np.log10((wl**2)/3e18)
	
	return AB_mag

def registerBandpasses():

	L_data = np.loadtxt('bandpasses/gotoL')
	R_data = np.loadtxt('bandpasses/gotoR')
	G_data = np.loadtxt('bandpasses/gotoG')
	B_data = np.loadtxt('bandpasses/gotoB')

	L_band = sncosmo.Bandpass(L_data[:,0], L_data[:,1], name = 'gotoL')
	R_band = sncosmo.Bandpass(R_data[:,0], R_data[:,1], name = 'gotoR')
	G_band = sncosmo.Bandpass(G_data[:,0], G_data[:,1], name = 'gotoG')
	B_band = sncosmo.Bandpass(B_data[:,0], B_data[:,1], name = 'gotoB')
	
	sncosmo.register(L_band)
	sncosmo.register(R_band)
	sncosmo.register(G_band)
	sncosmo.register(B_band)

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
	
		plt.plot(new_wl, new_flux, label = '%.3f' %(row['JD'] - explosion_epoch))

	plt.legend()
	plt.show()

	fluxes = np.array(fluxes, dtype = float)
	obj_phase = np.array(obj_phase, dtype = float)

	source = sncosmo.TimeSeriesSource(phase = obj_phase, wave = new_wl, flux = fluxes, name = 'kilonova_2017gfo')

# 	print(source)
# 
# 	model = sncosmo.Model(source)
# 	model.update({'z': redshift, 't0': obj_phase[0]})
# 	
# 	# model.set_source_peakabsmag(-15.6, 'desr', 'ab')
# 
# 	tobs = np.linspace(0., 15., 50)
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


def collectTemporalOverlap(query_output, explosion_epoch, duration):

	partial_query_output = query_output.query('jd > {} & jd < {}'.format(explosion_epoch, explosion_epoch + duration))

# 	print('### partial query, post-temp. filtering')
# 	print(partial_query_output[['jd', 'ra_c', 'dec_c']])

	return partial_query_output

def collectFootprintOverlap(query_output, ra_trans, dec_trans):

	chip_height = 2.8 # degrees
	chip_width = 2.1  # degrees
	
	max_allowed_ra = 360.
	min_allowed_ra = 0.

	upper_transient_ra_bound = ra_trans + (chip_width * np.cos(dec_trans*np.pi/180.)) / 2.
	lower_transient_ra_bound = ra_trans - (chip_width * np.cos(dec_trans*np.pi/180.)) / 2.
	
	upper_transient_dec_bound = dec_trans + chip_height / 2.
	lower_transient_dec_bound = dec_trans - chip_height / 2.
	
	# Filter RA first as it wraps (360 --> 0 degrees)
	if upper_transient_ra_bound > max_allowed_ra:
	
		partial_query_output = query_output.query('ra_c >= %f | ra_c <= %f' %(lower_transient_ra_bound, (upper_transient_ra_bound - max_allowed_ra)) )
	
	elif lower_transient_ra_bound < min_allowed_ra:
	
		partial_query_output = query_output.query('ra_c >= %f | ra_c <= %f' %((lower_transient_ra_bound + max_allowed_ra), (upper_transient_ra_bound - max_allowed_ra)) )
	
	else:

		partial_query_output = query_output.query('ra_c >= %f & ra_c <= %f' %(lower_transient_ra_bound, upper_transient_ra_bound) )
	
	# Now filter by DEC
	partial_query_output = partial_query_output.query('dec_c >= %f & dec_c <= %f' %(lower_transient_dec_bound, upper_transient_dec_bound) )

# 	print('### partial query, post-ra/dec filtering')
# 	print(partial_query_output[['jd', 'ra_c', 'dec_c']])

	return partial_query_output

def collectFilterOverlap(query_output, flt):

	flt_dict = {'gotoL': 'L', 'gotoR': 'R', 'gotoG': 'G', 'gotoB': 'B'}
	
	partial_query_output = query_output.query('filter == "%s"' %(flt_dict[flt]))
	
# 	print('### partial query, post-flt filtering')
# 	print(partial_query_output[['jd', 'ra_c', 'dec_c', 'filter']])

	return partial_query_output

def collectFootprintOverlapLoop(query_output, ra_trans, dec_trans):

# 	partial_query_output = query_output.query('ra_sw > {} & ra_nw > {} & ra_se < {} & ra_ne < {}'.format(ra_trans, ra_trans, ra_trans, ra_trans))

	ra_sw = np.array(query_output['ra_sw'])
	ra_se = np.array(query_output['ra_se'])
	dec_nw = np.array(query_output['dec_nw'])
	dec_sw = np.array(query_output['dec_sw'])
	
	for i in range(0, len(ra_sw)):
		
		if dec_sw[i] < dec_trans and dec_nw[i] > dec_trans:
			if ra_sw[i] < ra_trans and ra_se[i] > ra_trans:
				print(ra_sw[i], ra_trans, ra_se[i], dec_sw[i], dec_trans, dec_nw[i])

	
	return None

def generateSNLightcurvePopulation(query_output, population_settings, lightcurve_settings, source):

	tlc0 = time.time()

	transient_type = lightcurve_settings['population']['transient type']
	
	the_source = {'SN Ia': 'nugent-sn1a',
				  'SN Ib/c': 'nugent-sn1bc',
				  'SN II-P': 'nugent-sn2p',
				  'SN II-L': 'nugent-sn2l',
				  'SN IIn': 'nugent-sn2n',
				  'Kilonova': 'kilonova_2017gfo'}
					
	if not transient_type in the_source.keys():
		
		print('Transient type not recognised. Please use one of the following:\n')
		
		for k in the_source.keys():
			print(' *\t', k)
		
		sys.exit()
	
	the_duration = {'SN Ia': 95.,
					'SN Ib/c': 95.,
					'SN II-P': 120.,
					'SN II-L': 140.,
					'SN IIn': 130.,
					'Kilonova': 10.}
					
	the_peak = {'SN Ia': -19.,
				'SN Ib/c': -18.,
				'SN II-P': -17.5,
				'SN II-L': -17.,
				'SN IIn': -17.5,
				'Kilonova': -16.}

	dust = sncosmo.CCM89Dust()
	dustmap = sfdmap.SFDMap('sfdmaps/sfddata-master')
	
	if source:

		model = sncosmo.Model(source = source, effects = [dust, dust], effect_names = ['host', 'mw'], effect_frames = ['rest', 'obs'])
	
	else:
	
		model = sncosmo.Model(source = the_source[transient_type], effects = [dust, dust], effect_names = ['host', 'mw'], effect_frames = ['rest', 'obs'])
	
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
	
		model.set_source_peakabsmag(the_peak[transient_type], 'gotoL', 'ab')
		
		partial_query_output = collectTemporalOverlap(query_output, t0, the_duration[transient_type])
		partial_query_output = collectFootprintOverlap(partial_query_output, ra, dec)

		transient_dict = {'explosion epoch': t0,
						  'distance': distance,
						  'redshift': z,
						  'ra': ra,
						  'dec': dec,
						  'ebv': ebv,
						  'lightcurves': {}}
		
		flt_dict = {}
		
		for flt in lightcurve_settings['filters']:
		
			flt_partial_query_output = collectFilterOverlap(partial_query_output, flt).sort_values(['jd'])
		
			tobs = np.array(flt_partial_query_output['jd'])
			mag5sig = np.array(flt_partial_query_output['limmag5'])
			
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
	
	with open(lightcurve_settings['lightcurve file'] + '.json', 'w') as f:
		json.dump(main_lc_dict, f, indent = 4)

	return main_lc_dict

def generateGRBLightcurvePopulation(query_output, population_settings, lightcurve_settings, grb_settings):
	
	jet_type_dict = {'TopHat': grb.jet.TopHat,
					 'Gaussian': grb.jet.Gaussian,
					 'PowerLaw': grb.jet.PowerLaw}
	
	goto_wleffs = {'gotoL': 5396.65,
				   'gotoR': 6405.05,
				   'gotoG': 5355.32, 
				   'gotoB': 4600.58}
				   
	c = 3.e18
	seconds_in_day = 8.64e4
	step_sizen = 25
	
	goto_nueffs = {'gotoL': c / 5396.65,
				   'gotoR': c / 6405.05,
				   'gotoG': c / 5355.32, 
				   'gotoB': c / 4600.58}
	
	goto_colours = {'gotoL': 'goldenrod',
					'gotoR': 'firebrick',
					'gotoG': 'limegreen', 
					'gotoB': 'royalblue'}
	
	grb_settings['jetType'] = jet_type_dict[grb_settings['jetType']]
	
	main_lc_dict = {}
	
	dustmap = sfdmap.SFDMap('sfdmaps/sfddata-master')

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
		
		partial_query_output = collectTemporalOverlap(query_output, t0, 350.)
		partial_query_output = collectFootprintOverlap(partial_query_output, ra, dec)

		transient_dict = {'explosion epoch': t0,
						  'distance': distance,
						  'redshift': z,
						  'ra': ra,
						  'dec': dec,
						  'ebv': ebv,
						  'lightcurves': {}}
		
		flt_dict = {}
		
		for flt in lightcurve_settings['filters']:
		
			flt_partial_query_output = collectFilterOverlap(partial_query_output, flt).sort_values(['jd'])
		
			tobs = (np.array(flt_partial_query_output['jd']) - t0) * seconds_in_day
			mag5sig = np.array(flt_partial_query_output['limmag5'])
			
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
			
# 				nu = np.full(len(tobs), c / goto_wleffs[flt])
# 				flux_mJy = grb.fluxDensity(tobs, nu, **grb_settings)
# 				mags = mJy2ABmag(flux_mJy, goto_wleffs[flt])
				mags = mJy2ABmag(flux_mJy_total, goto_wleffs[flt])
				
				A_ex = extinction.fitzpatrick99(np.array([goto_wleffs[flt]]), 3.1 * ebv)
# 				print(A_ex)

				ind = np.where(mags < mag5sig)
				tobs_rec = tobs[ind]
				mags_rec = mags[ind]
			
				lc_dict = {'tobs': list(tobs / seconds_in_day),
						   'mags': list(mags + A_ex[0]),
						   'limmag5': list(mag5sig),
						   'tobs_rec': list(tobs_rec),
						   'mags_rec': list(mags_rec)}
			
				flt_dict.update({flt: lc_dict})


		transient_dict['lightcurves'].update(flt_dict)
		
		main_lc_dict.update({transient_index_str: transient_dict})
	
	return main_lc_dict

# if __name__ == '__main__':

	