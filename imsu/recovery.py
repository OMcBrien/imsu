import sys
import json
import numpy as np

def nightlyCount(tobs):

	night_count = np.unique(np.floor(tobs))

	return list(night_count)

def performMinimumDetectionCount(lc_population, recovery_settings):
	
# 	print(recovery_settings)
	
	for transient, properties in lc_population.items():
	
# 		print('\n### {} ###'.format(transient))
		min_det_count_dict = {'min det count check': True}

		for flt, lightcurves in properties['lightcurves'].items():
# 			print('Filter = ', flt)
# 			print('must pass min det. count ==> ', recovery_settings['min detection counts'][flt])
# 			print('min det. count = ', len(lightcurves['mags_rec']))
			
			if len(lightcurves['mags_rec']) < recovery_settings['min detection counts'][flt]:
				
# 				print('non recovery')
				min_det_count_dict = {'min det count check': False}
				break
			
# 			print('###')
	
		properties.update(min_det_count_dict)
# 		print(properties.keys())
		
# 	print(lc_population)

	return lc_population

def performMinimumNightlyDetectionCount(lc_population, recovery_settings):
	
	for transient, properties in lc_population.items():
		
# 		print('\n### {} ###'.format(transient))
		min_nightly_det_count_dict = {'min nightly det count check': True}

		for flt, lightcurves in properties['lightcurves'].items():
# 			print('Filter = ', flt)
# 			print('must pass min nightly det. count ==> ', recovery_settings['min nightly detection counts'][flt])
		
			tobs_rec = lightcurves['tobs_rec']
			
			tobs_rec_night_count = nightlyCount(tobs_rec)
			
# 			print('nightly det. count = ', len(tobs_rec_night_count))
			
			if len(tobs_rec_night_count) < recovery_settings['min nightly detection counts'][flt]:
				
# 				print('non recovery')
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

def saveRecoveryData(lc_population, lightcurve_settings, population_settings):

	with open(lightcurve_settings['lightcurve file'] + '_rec.json', 'w') as f:
		json.dump(lc_population, f, indent = 4)
	
	population_settings.to_csv(lightcurve_settings['population']['population file'] + '_rec.csv', index = False)



