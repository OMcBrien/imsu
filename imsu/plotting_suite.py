import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def make2DHistogramGOTOPointings(query_output, lightcurve_settings):
	
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex = True, sharey = True, figsize = (15,12))
	
	df_L = query_output.where(query_output['filter'] == 'L').dropna(how = 'all')
	df_R = query_output.where(query_output['filter'] == 'R').dropna(how = 'all')
	df_G = query_output.where(query_output['filter'] == 'G').dropna(how = 'all')
	df_B = query_output.where(query_output['filter'] == 'B').dropna(how = 'all')
	
	panel1 = ax1.hist2d(np.array(df_L['ra_c']), np.array(df_L['dec_c']), bins = [171, 64], cmap = plt.cm.Oranges)
	panel2 = ax2.hist2d(np.array(df_R['ra_c']), np.array(df_R['dec_c']), bins = [171, 64], cmap = plt.cm.Reds)
	panel3 = ax3.hist2d(np.array(df_G['ra_c']), np.array(df_G['dec_c']), bins = [171, 64], cmap = plt.cm.Greens)
	panel4 = ax4.hist2d(np.array(df_B['ra_c']), np.array(df_B['dec_c']), bins = [171, 64], cmap = plt.cm.Blues)
	
	fig.colorbar(panel1[3], ax = ax1)
	fig.colorbar(panel2[3], ax = ax2)
	fig.colorbar(panel3[3], ax = ax3)
	fig.colorbar(panel4[3], ax = ax4)
	
	ax1.set_xlim([0,360])
	ax1.set_ylim([-90,90])
	
	ax1.set_ylabel('Dec, deg.')
	ax3.set_ylabel('Dec, deg.')
	ax3.set_xlabel('RA, deg.')
	ax4.set_xlabel('RA, deg.')
	
	plt.tight_layout()
	plt.show()

def plotLimitingMagnitudesFromQuery(query_output, lightcurve_settings):

	goto_colours = {'gotoL': 'goldenrod', 'gotoR': 'firebrick', 'gotoG': 'limegreen', 'gotoB': 'royalblue'}
	
	for flt in lightcurve_settings['filters']:
		
		query_output_toplot = query_output.where(query_output['filter'] == flt[-1])
		query_output_toplot.dropna()
		
		plt.plot(query_output_toplot['jd'], query_output_toplot['limmag5'], marker = 'v', ms = 8, mec = 'black', mew = 0.5, mfc = goto_colours[flt], ls = 'None', alpha = 0.5, label = flt)
	
	plt.xlabel('JD')
	plt.ylabel('5$\sigma$ limiting mag.')
	
	plt.legend(loc = 'best', frameon = False, ncol = 2)
	
	plt.gca().invert_yaxis()
		
	plt.tight_layout()
	plt.show()

def makeRedshiftHistogram(population_settings, survey_settings):

	df_rec = population_settings.where(population_settings['detected'] == True)
	df_nrec = population_settings.where(population_settings['detected'] == False)
	
	the_bins = np.linspace(survey_settings['depth']['min redshift'], survey_settings['depth']['max redshift'], 15)
	
	plt.hist([df_rec['redshifts'], df_nrec['redshifts']],
				bins = the_bins,
				histtype = 'barstacked',
				stacked = True,
				edgecolor = 'black',
				color = ['limegreen', 'lightcoral'],
				linewidth = 2,
				label = ['Recovered', 'Non-recovered'])
	
	plt.legend(loc = 'upper left', frameon = False)
	
	plt.xlabel('Redshift, $z$')
	plt.xticks(the_bins, rotation = 315)
	
	plt.tight_layout()
	plt.show()

def makeExplosionEpochHistogram(population_settings, survey_settings):

	df_rec = population_settings.where(population_settings['detected'] == True)
	df_nrec = population_settings.where(population_settings['detected'] == False)
	
	the_bins = np.linspace(survey_settings['time observing']['survey begin'], survey_settings['time observing']['survey end'], 15)
	
	plt.hist([df_rec['explosion epochs'], df_nrec['explosion epochs']],
				bins = the_bins,
				histtype = 'barstacked',
				stacked = True,
				edgecolor = 'black',
				color = ['limegreen', 'lightcoral'],
				linewidth = 2,
				label = ['Recovered', 'Non-recovered'])
	
	plt.legend(loc = 'upper left', frameon = False)
	
	plt.xlabel('Explosion epoch, JD')
	plt.xticks(the_bins, rotation = 315)
	
	plt.tight_layout()
	plt.show()

def make2DHistogramRecoveredPointings(population_settings, survey_settings):
	
	df_rec = population_settings.where(population_settings['detected'] == True).dropna(how = 'all')
	df_nrec = population_settings.where(population_settings['detected'] == False).dropna(how = 'all')
	
	plt.hist2d(df_rec['ras'], df_rec['decs'], bins = [40, 20], cmap = plt.cm.Greens)
	plt.colorbar()
	plt.title('Recovered')
	plt.xlabel('RA, deg.')
	plt.ylabel('Dec, deg.')
	plt.xlim([0., 360.])
	plt.ylim([-90., 90.])
	
	plt.tight_layout()
	plt.show()
	
	plt.hist2d(df_nrec['ras'], df_nrec['decs'], bins = [40, 20], cmap = plt.cm.Reds)
	plt.colorbar()
	plt.title('Non-recovered')
	plt.xlabel('RA, deg.')
	plt.ylabel('Dec, deg.')
	plt.xlim([0., 360.])
	plt.ylim([-90., 90.])
	
	plt.tight_layout()
	plt.show()

def plotTransientLightcurves(lc_population, lightcurve_settings):

	goto_colours = {'gotoL': 'goldenrod', 'gotoR': 'firebrick', 'gotoG': 'limegreen', 'gotoB': 'royalblue'}
	counter = 0
	
	for transient, properties in lc_population.items():
		
		lc_properties = properties['lightcurves']
		
		for flt in lightcurve_settings['filters']:
			
			tobs = lc_properties[flt]['tobs']
			mags = lc_properties[flt]['mags']
			limmag5 = lc_properties[flt]['limmag5']
			tobs_rec = lc_properties[flt]['tobs_rec']
			mags_rec = lc_properties[flt]['mags_rec']
            
			plt.plot(tobs, mags, ls = 'None', marker = 'o', ms = 8, mec = 'k', label = '{}'.format(flt), mfc = goto_colours[flt], linewidth = 2)
			plt.plot(tobs, limmag5, ls = 'None', marker = 'v', ms = 8, mec = 'k', mfc = goto_colours[flt])
			plt.plot(tobs_rec, mags_rec, ls = 'None', marker = 'o', ms = 14, mfc = 'None', mec = 'red', mew = 2)
		
		plt.legend(loc = 'best', frameon = False, ncol = 2)

		plt.xlabel('Phase, days')
		plt.ylabel('Apparent magnitude')
		plt.title(transient)
	
		plt.tight_layout()
		plt.gca().invert_yaxis()
		plt.show()
		
		counter += 1
		
		if counter >= 25:
			break
	
	
	
	
	
	
	
	
	