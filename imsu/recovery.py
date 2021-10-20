import sys
import json
import numpy as np

def nightlyCount(tobs):

	night_count = np.unique(np.floor(tobs))

	return list(night_count)

def performMinimumDetectionCount(lc_population, recovery_settings):
	
	for transient, properties in lc_population.items():
	
		min_det_count_dict = {'min det count check': True}

		for flt, lightcurves in properties['lightcurves'].items():

			if len(lightcurves['mags_rec']) < recovery_settings['min detection counts'][flt]:
				
				min_det_count_dict = {'min det count check': False}
				break
	
		properties.update(min_det_count_dict)

	return lc_population

def performMinimumNightlyDetectionCount(lc_population, recovery_settings):
	
	for transient, properties in lc_population.items():
		
		min_nightly_det_count_dict = {'min nightly det count check': True}

		for flt, lightcurves in properties['lightcurves'].items():
		
			tobs_rec = lightcurves['tobs_rec']
			
			tobs_rec_night_count = nightlyCount(tobs_rec)
			
			if len(tobs_rec_night_count) < recovery_settings['min nightly detection counts'][flt]:
				
				min_nightly_det_count_dict = {'min nightly det count check': False}
				break
				
		properties.update(min_nightly_det_count_dict)
	
	return lc_population

def countRecoveries(lc_population, recovery_settings, population_settings):

	detection_count = 0
	
	detection_flag = [False for item in population_settings['transient_id']]

	for index, (transient, properties) in enumerate(lc_population.items()):
	
# 		print(properties['min det count check'], properties['min nightly det count check'])
	
		if properties['min det count check'] and properties['min nightly det count check']:
		
			detection_count += 1
			detection_flag[index] = True
	
# 	print(detection_count)
	population_settings['detected'] = detection_flag
	print(population_settings)
	
	return detection_count

def saveRecoveryData(lc_population, population_settings, recovery_lc_path, recovery_pop_path):

	with open(recovery_lc_path, 'w') as f:
		json.dump(lc_population, f, indent = 4)
	
	population_settings.to_csv(recovery_pop_path, index = False)



