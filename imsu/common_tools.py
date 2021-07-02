import matplotlib
import math
import re
import sys

def readSettings(path_to_settings = 'settings.yaml'):

	import yaml
	
	with open(path_to_settings, 'r') as stream:
		
		loader = yaml.SafeLoader
		loader.add_implicit_resolver(u'tag:yaml.org,2002:float',
									re.compile(u'''^(?:
								    [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
								    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
								    |\\.[0-9_]+(?:[eE][-+][0-9]+)?
								    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
								    |[-+]?\\.(?:inf|Inf|INF)
								    |\\.(?:nan|NaN|NAN))$''', re.X),
								    list(u'-+0123456789.'))
	
		all_settings = yaml.load(stream, Loader = loader)

	
	program_settings = all_settings['Program Settings']
	query_settings = all_settings['Query Settings']
	survey_settings = all_settings['Survey Settings']
	lightcurve_settings = all_settings['Lightcurve Settings']
	grb_settings = all_settings['GRB Settings']
	recovery_settings = all_settings['Recovery Settings']
	
	return (all_settings,
			program_settings,
			query_settings,
			survey_settings,
			lightcurve_settings,
			grb_settings,
			recovery_settings)

def setPlottingParameters():

	"""
	Plot Metadata
	"""

	SMALL_SIZE = 15
	MEDIUM_SIZE = 20
	BIGGER_SIZE = 25


	matplotlib.rcParams['figure.figsize'] = (10, 8)

	matplotlib.rcParams['font.size'] = SMALL_SIZE		
	matplotlib.rcParams['axes.titlesize'] = BIGGER_SIZE  
	matplotlib.rcParams['axes.labelsize'] = MEDIUM_SIZE 
	matplotlib.rcParams['xtick.labelsize'] = SMALL_SIZE 
	matplotlib.rcParams['ytick.labelsize'] = SMALL_SIZE 
	matplotlib.rcParams['legend.fontsize'] = MEDIUM_SIZE
	matplotlib.rcParams['figure.titlesize'] = BIGGER_SIZE

	matplotlib.rcParams["font.family"] = "serif"
	matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'


	"""
	Script begin
	"""