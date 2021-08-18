import common_tools, db_queries, population_generator, lc_generator, recovery, plotting_suite
import sys
import os
import json
import pandas as pd

# =============================

def main():

	(all_settings,
	program_settings,
	query_settings,
	lightcurve_settings,
	grb_settings,
	recovery_settings) = common_tools.readSettings()
	
	common_tools.setPlottingParameters()
	
	if program_settings['Query Database']:

		db_conn = db_queries.establishConnection()
		query_output = db_queries.queryDB(db_conn, query_settings)

	else:
		
		if not os.path.exists(query_settings['query file'] + '.csv'):
			sys.exit('No query results exist. Enable query execution.')
		else:
			query_output = pd.read_csv(query_settings['query file'] + '.csv', dtype = {'dbevent': str, 'tilename': str})
	
	the_survey = common_tools.BoxSurvey(query_output)
	
	print(the_survey.begin_jd)
	print(the_survey.end_jd)
	print(the_survey.lower_ra_bound)
	print(the_survey.upper_ra_bound)
	print(the_survey.lower_dec_bound)
	print(the_survey.upper_dec_bound)
	
	if program_settings['Generate Transients']:
		
		population_settings = population_generator.generatePopulation(the_survey, lightcurve_settings)
		
		population_settings = pd.DataFrame(population_settings)
		population_settings.to_csv(lightcurve_settings['population']['population file'] + '.csv', index = False)

	else:

		if not os.path.exists(lightcurve_settings['population']['population file'] + '.csv'):
			sys.exit('No population file exists. Enable population generator.')
		else:
			population_settings = pd.read_csv(lightcurve_settings['population']['population file'] + '.csv')

	
	if program_settings['Generate Lightcurves']:
	
		if lightcurve_settings['population']['transient type'] == 'GRB':
			
			lc_population = lc_generator.generateGRBLightcurvePopulation(query_output, population_settings, lightcurve_settings, grb_settings)
		
		elif lightcurve_settings['population']['transient type'] == 'Kilonova':
			
			lc_generator.registerBandpasses()
			source = lc_generator.addTimeSeriesSource(obj_name = 'at2017gfo', explosion_epoch = 2457983.48, redshift = 0.00984, min_wl = 3400., max_wl = 22300., npoints = 40000)
			lc_population = lc_generator.generateSNLightcurvePopulation(query_output, population_settings, lightcurve_settings, source)
		
		else:

			lc_generator.registerBandpasses()
			lc_population = lc_generator.generateSNLightcurvePopulation(query_output, population_settings, lightcurve_settings, source = None)
	
	else:
	
		if not os.path.exists(lightcurve_settings['lightcurve file'] + '.json'):
			sys.exit('No lightcurve population has been generated. Enable lightcurve generator.')
		else:
			with open(lightcurve_settings['lightcurve file'] + '.json', 'r') as f:
				lc_population = json.load(f)

	
	if program_settings['Perform Recovery']:
			
		lc_population = recovery.performMinimumDetectionCount(lc_population, recovery_settings)
		lc_population = recovery.performMinimumNightlyDetectionCount(lc_population, recovery_settings)
			
		recovered_count = recovery.countRecoveries(lc_population, recovery_settings, population_settings)
		print('Efficiency = %.3f per cent' %(recovered_count * 100 / lightcurve_settings['population']['number to inject'] ) )
			
		recovery.saveRecoveryData(lc_population, lightcurve_settings, population_settings)
		
	else:
		
		if not os.path.exists(lightcurve_settings['lightcurve file'] + '_rec.json'):
			sys.exit('Lightcurve recovery has not been performed. Enable recovery module.')
		else:
			with open(lightcurve_settings['lightcurve file'] + '_rec.json', 'r') as f:
				lc_population = json.load(f)
			
			population_settings = pd.read_csv(lightcurve_settings['population']['population file'] + '_rec.csv')
	
	if program_settings['Make Plots']:
		
		common_tools.setPlottingParameters()
		
		plotting_suite.make2DHistogramGOTOPointings(query_output, lightcurve_settings)
		plotting_suite.plotLimitingMagnitudesFromQuery(query_output, lightcurve_settings)
		
		plotting_suite.makeRedshiftHistogram(population_settings, lightcurve_settings)
		plotting_suite.makeExplosionEpochHistogram(population_settings, the_survey)
		
		plotting_suite.make2DHistogramRecoveredPointings(population_settings)

		plotting_suite.plotTransientLightcurves(lc_population, lightcurve_settings)
		
		
	


if __name__ == '__main__':
	main()
	