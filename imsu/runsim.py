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
	recovery_settings,
	io_settings) = common_tools.readSettings()
	
	common_tools.setPlottingParameters()
	
	output_path = os.path.join('results', io_settings['results directory'])
	
	if not os.path.exists('results'):
		os.mkdir('results')
		
	if not os.path.exists(output_path):
		os.mkdir(output_path)
		
	
	# Query database or read previous query from file
	query_path = os.path.join(output_path, query_settings['query file']) + '.csv'
	
	if program_settings['Query Database']:

		db_conn = db_queries.establishConnection()
		query_output = db_queries.queryDB(db_conn, query_settings)
		query_output.to_csv(query_path, index = False)

	else:
	
		if not os.path.exists(query_path):
			sys.exit('No query results exist. Enable query execution or provide history of pointings.')
		else:
			query_output = pd.read_csv(query_path)
	
	
	# Set special treatment for GRB simulations
	if lightcurve_settings['population']['transient type'] == 'GRB':
		is_grb = True
	else:
		is_grb = False
	
	# Establish Survey() class for given transient type
	if is_grb:
		the_survey = common_tools.EllipseSurvey(query_output, lightcurve_settings['telescope properties']['column headings'], lightcurve_settings['ellipse geometry'])
	else:
		the_survey = common_tools.BoxSurvey(query_output, lightcurve_settings['telescope properties']['column headings'])

	# Begin generating transient population (RA, Dec, z, etc.)	
	population_path = os.path.join(output_path, lightcurve_settings['population']['population file']) + '.csv'
	
	if program_settings['Generate Transients']:
	
		if is_grb:
			
			population_settings = population_generator.generateGRBPopulation(the_survey, lightcurve_settings)
		
			population_settings = pd.DataFrame(population_settings)
			population_settings.to_csv(population_path, index = False)			
		
		else:
		
			population_settings = population_generator.generatePopulation(the_survey, lightcurve_settings)
		
			population_settings = pd.DataFrame(population_settings)
			population_settings.to_csv(population_path, index = False)

	else:

		if not os.path.exists(population_path):
			sys.exit('No population file exists. Enable population generator.')
		else:
			population_settings = pd.read_csv(population_path)

	# Generate lightcurves from population
	lightcurve_path = os.path.join(output_path, lightcurve_settings['lightcurve file']) + '.csv'
	
	if program_settings['Generate Lightcurves']:
		
		lc_generator.registerBandpasses(lightcurve_settings)
	
		if lightcurve_settings['population']['transient type'] == 'GRB':
			
			lc_population = lc_generator.generateGRBLightcurvePopulation(query_output, population_settings, lightcurve_settings, grb_settings, io_settings)
			
			with open(lightcurve_path, 'w') as f:
				json.dump(lc_population, f, indent = 4)
		
		elif lightcurve_settings['population']['transient type'] == 'Kilonova':
			
			source = lc_generator.addTimeSeriesSource(obj_name = 'at2017gfo', explosion_epoch = 2457983.48, redshift = 0.00984, min_wl = 3400., max_wl = 22300., npoints = 40000)
			lc_population = lc_generator.generateSNLightcurvePopulation(query_output, population_settings, lightcurve_settings, io_settings, source)
			
			lc_generator.saveLightcurveData(lightcurve_path, lc_population)
		
		else:

			lc_population = lc_generator.generateSNLightcurvePopulation(query_output, population_settings, lightcurve_settings, io_settings, source = None)
			
			lc_generator.saveLightcurveData(lightcurve_path, lc_population)
	
	else:
	
		if not os.path.exists(lightcurve_path):
			sys.exit('No lightcurve population has been generated. Enable lightcurve generator.')
		else:
			with open(lightcurve_path, 'r') as f:
				lc_population = json.load(f)

	# Run recovery tests
	recovery_lc_path = os.path.join(output_path, lightcurve_settings['lightcurve file']) + '_rec.json'
	recovery_pop_path = os.path.join(output_path, lightcurve_settings['lightcurve file']) + '_rec.csv'
	
	
	if program_settings['Perform Recovery']:
			
		lc_population = recovery.performMinimumDetectionCount(lc_population, recovery_settings)
		lc_population = recovery.performMinimumNightlyDetectionCount(lc_population, recovery_settings)
			
		recovered_count = recovery.countRecoveries(lc_population, recovery_settings, population_settings)
		print('Efficiency = %.3f per cent' %(recovered_count * 100 / lightcurve_settings['population']['number to inject'] ) )
			
		recovery.saveRecoveryData(lc_population, population_settings, recovery_lc_path, recovery_pop_path)
		
	else:
		
		if not os.path.exists(lightcurve_settings['lightcurve file'] + '_rec.json'):
			sys.exit('Lightcurve recovery has not been performed. Enable recovery module.')
		else:
			with open(lightcurve_settings['lightcurve file'] + '_rec.json', 'r') as f:
				lc_population = json.load(f)
			
			population_settings = pd.read_csv(lightcurve_settings['population']['population file'] + '_rec.csv')
	
	if program_settings['Make Plots']:
		
		common_tools.setPlottingParameters()
		
# 		plotting_suite.make2DHistogramGOTOPointings(query_output, lightcurve_settings)
# 		plotting_suite.plotLimitingMagnitudesFromQuery(query_output, lightcurve_settings)
		
		plotting_suite.makeRedshiftHistogram(population_settings, lightcurve_settings)
		plotting_suite.makeExplosionEpochHistogram(population_settings, the_survey)
		
		plotting_suite.make2DHistogramRecoveredPointings(population_settings)

# 		plotting_suite.plotTransientLightcurves(lc_population, lightcurve_settings)
# 		plotting_suite.plotRecoveredLightcurves(lc_population, lightcurve_settings, population_settings)
		
		
	


if __name__ == '__main__':
	main()
	