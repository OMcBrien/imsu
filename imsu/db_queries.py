import gotocat as gc
import time
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def establishConnection(phase = 4):

	db_conn = gc.GOTOdb(phase = phase)

	if not type(db_conn.conn) == gc.psycopg2.extensions.connection:
		sys.exit("Database connection not created, check your connection to {}".format(d4.host))

	return db_conn

def queryDB(db_conn, query_settings):

	t0 = time.time()

	if not query_settings['dbevent']:
	
		survey_begin = query_settings['time observing']['survey begin']
		survey_end = query_settings['time observing']['survey end']
	
		query_output = db_conn.query("SELECT * FROM image WHERE image.jd BETWEEN {} and {}".format(survey_begin, survey_end))
		query_output.to_csv('{}.csv'.format(query_settings['query file']), index = False)

	else:

		query_output = db_conn.query("SELECT * FROM image WHERE image.dbevent = '{}'".format(query_settings['dbevent']))
		
	t1 = time.time()
	
	print('Query ran in %.4f s' %(t1 - t0))

	return query_output

# Going to need to move this snippet elsewhere
def transformRightAscension(ra, origin = 0):

	# Fixing plot axes for sky projection

	RA = np.remainder(ra + 360 - origin, 360)
	ind = RA > 180
	RA[ind] -= 360
	RA = -RA
	# That's better
	
	return RA

def visualiseFootprints(query_output):

	goto_colours = {'L': 'goldenrod', 'R': 'firebrick', 'G': 'limegreen', 'B': 'royalblue'}

	fig = plt.figure()
	ax = plt.subplot(111)
	ax.grid(True)

	ra_se = np.array(query_output['ra_se'])
	ra_sw = np.array(query_output['ra_sw'])
	dec_se = np.array(query_output['dec_se'])
	dec_nw = np.array(query_output['dec_nw'])
	dec_sw = np.array(query_output['dec_sw'])
	dec_c = np.array(query_output['dec_c'])
	ra_c = np.array(query_output['ra_c'])
	flt = np.array(query_output['filter'])
	

	footprint_width = ra_sw - ra_se
	footprint_height = dec_nw - dec_sw
	
	ind = np.where(footprint_width > 350.)
	footprint_width[ind] = -1. * (footprint_width[ind] - 360.)
	
	ax.scatter(ra_c, dec_c, s = 45, marker = 'x', c = 'red')

	for i in range(0, len(ra_se)):
		
		new_footprint = Rectangle( (ra_se[i], dec_se[i]), width = footprint_width[i], height = footprint_height[i], alpha = 0.25, facecolor = goto_colours[flt[i]], edgecolor = goto_colours[flt[i]], linewidth = 2)		
		ax.add_patch(new_footprint)

	
	ax.set_xlim([0,360])
	ax.set_ylim([-90,90])
	
	ax.set_xlabel('RA (deg.)')
	ax.set_ylabel('Dec (deg.)')
	
	plt.tight_layout()
	plt.show()

	return None
	




