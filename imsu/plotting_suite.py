import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def makeRedshiftHistogram(population_settings, lightcurve_settings):

	the_bins = np.linspace(lightcurve_settings['depth']['min redshift'], lightcurve_settings['depth']['max redshift'], 15)
	
	for flt in lightcurve_settings['telescope properties']['filter effective wavelengths'].keys():
	
		try:
		
			df_rec = population_settings.where(population_settings['detected_' + flt] == True)
			df_nrec = population_settings.where(population_settings['detected_' + flt] == False)
		
		except RuntimeWarning:
		
			continue
	
		plt.hist([df_rec['redshifts'], df_nrec['redshifts']],
					bins = the_bins,
					histtype = 'barstacked',
					stacked = True,
					edgecolor = 'black',
					color = ['limegreen', 'lightcoral'],
					linewidth = 2,
					label = ['Recovered', 'Non-recovered'])
	
		plt.legend(loc = 'upper left', frameon = False)
		
		plt.title(flt)
		
		plt.xlabel('Redshift, $z$')
		plt.xticks(the_bins, rotation = 315)
	
		plt.tight_layout()
		plt.show()

def makeExplosionEpochHistogram(population_settings, the_survey, lightcurve_settings):

	the_bins = np.linspace(the_survey.begin_jd, the_survey.end_jd, 15)

	for flt in lightcurve_settings['telescope properties']['filter effective wavelengths'].keys():
	
		try:
		
			df_rec = population_settings.where(population_settings['detected_' + flt] == True)
			df_nrec = population_settings.where(population_settings['detected_' + flt] == False)
		
		except RuntimeWarning:
			
			continue
	
		plt.hist([df_rec['explosion epochs'], df_nrec['explosion epochs']],
					bins = the_bins,
					histtype = 'barstacked',
					stacked = True,
					edgecolor = 'black',
					color = ['limegreen', 'lightcoral'],
					linewidth = 2,
					label = ['Recovered', 'Non-recovered'])
	
		plt.legend(loc = 'upper left', frameon = False)
		
		plt.title(flt)
		
		plt.xlabel('Explosion epoch, JD')
		plt.xticks(the_bins, rotation = 315)
	
		plt.tight_layout()
		plt.show()

def make2DHistogramRecoveredPointings(population_settings, lightcurve_settings):

	for flt in lightcurve_settings['telescope properties']['filter effective wavelengths'].keys():
		
		try:
		
			df_rec = population_settings.where(population_settings['detected_' + flt] == True).dropna(how = 'all')
			df_nrec = population_settings.where(population_settings['detected_' + flt] == False).dropna(how = 'all')
		
		except RuntimeWarning:
			
			continue
	
		plt.hist2d(df_rec['ras'], df_rec['decs'], bins = [40, 20], cmap = plt.cm.Greens)
		plt.colorbar()
		plt.title(flt, ' recovered')
		plt.xlabel('RA, deg.')
		plt.ylabel('Dec, deg.')
		plt.xlim([0., 360.])
		plt.ylim([-90., 90.])
	
		plt.tight_layout()
		plt.show()
	
		plt.hist2d(df_nrec['ras'], df_nrec['decs'], bins = [40, 20], cmap = plt.cm.Reds)
		plt.colorbar()
		plt.title(flt, ' non-recovered')
		plt.xlabel('RA, deg.')
		plt.ylabel('Dec, deg.')
		plt.xlim([0., 360.])
		plt.ylim([-90., 90.])
	
		plt.tight_layout()
		plt.show()

	
	
	
	
	
	
	
	
	