import sys
import json
import numpy as np

def nightlyCount(tobs):

	night_count = np.unique(np.floor(tobs))

	return list(night_count)
	
def initialiseRecoveryProcedure(lc_population):

	for transient, properties in lc_population.items():
	
		properties.update({'recovery': {}})	

	return lc_population

def performMinimumDetectionCount(lc_population, recovery_settings):
	
	for transient, properties in lc_population.items():
	
		min_det_count_dict = {'min detection counts': {}}

		for flt, lightcurves in properties['lightcurves'].items():

			if len(lightcurves['mags_rec']) < recovery_settings['min detection counts'][flt]:
				
				min_det_count_dict['min detection counts'].update({flt: False})
# 				break
			else:

				min_det_count_dict['min detection counts'].update({flt: True})
			
		properties['recovery'].update(min_det_count_dict)

	return lc_population

def performMinimumNightlyDetectionCount(lc_population, recovery_settings):
	
	for transient, properties in lc_population.items():
		
		min_nightly_det_count_dict = {'min nightly detection counts': {}}

		for flt, lightcurves in properties['lightcurves'].items():
		
			tobs_rec = lightcurves['tobs_rec']
			
			tobs_rec_night_count = nightlyCount(tobs_rec)
			
			if len(tobs_rec_night_count) < recovery_settings['min nightly detection counts'][flt]:
				
				min_nightly_det_count_dict['min nightly detection counts'].update({flt: False})
# 				break
			else:
			
				min_nightly_det_count_dict['min nightly detection counts'].update({flt: True})
				
		properties['recovery'].update(min_nightly_det_count_dict)
	
	return lc_population

def markRecoveriesByFilter(lc_population, recovery_settings, population_settings):

	for flt in recovery_settings['min detection counts'].keys():

		detection_count = 0	
		detection_flag = [False for item in population_settings['transient_id']]

		for index, (transient, properties) in enumerate(lc_population.items()):
	
			if properties['recovery']['min detection counts'][flt] and properties['recovery']['min nightly detection counts'][flt]:
		
				detection_count += 1
				detection_flag[index] = True
	
	# 	print(detection_count)
		population_settings['detected_' + flt] = detection_flag
	
	return population_settings
	
def calculateEfficiencies(recovery_settings, population_settings):
	
	filters_in_use = recovery_settings['min detection counts'].keys()
	
	for flt in filters_in_use:
		
		detected_status = population_settings['detected_' + flt]
		
		print('Efficiency in %s:\t%.3f' %(flt, np.sum(detected_status) * 100 / len(detected_status)))

def saveRecoveryData(lc_population, population_settings, recovery_lc_path, recovery_pop_path):

	with open(recovery_lc_path, 'w') as f:
		json.dump(lc_population, f, indent = 4)
	
	population_settings.to_csv(recovery_pop_path, index = False)



