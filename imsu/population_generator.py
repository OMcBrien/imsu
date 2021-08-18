import sncosmo
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import time

## LAST MODIFIED : April 15, 2013
## CREATED : April 15, 2013
## AUTHOR : DRYX
def redshift2distance(z):
    """Convert a redshift to various distance units

    **Key Arguments:**
        - ``z`` -- the redshift to be converted

    **Return:**
        - ``results`` -- result dictionary including
            - ``dcmr_mpc`` -- co-moving radius distance
            - ``da_mpc`` -- angular distance
            - ``da_scale`` -- angular distance scale
            - ``dl_mpc`` -- luminosity distance (usually use this one)
            - ``dmod`` -- distance modulus (determined from luminosity distance)
    """
    ################ > IMPORTS ################
    ## STANDARD LIB ##
    ## THIRD PARTY ##
    ## LOCAL APPLICATION ##

    ################ >ACTION(S) ################
    # Cosmological Parameters (to be changed if required)
    WM = 0.3           # Omega_matter
    WV = 0.7           # Omega_vacuum
    H0 = 70.0           # Hubble constant (km s-1 Mpc-1)

    # Other variables
    h = H0/100.0
    WR = 4.165E-5/(h*h)     # Omega_radiation
    WK = 1.0-WM-WV-WR       # Omega_curvature = 1 - Omega(Total)
    c = 299792.458          # speed of light (km/s)

    # Arbitrarily set the values of these variables to zero just so we can define them.
    DCMR = 0.0             # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0          # comoving radial distance in units of Mpc
    DA = 0.0                # angular size distance in units of c/H0
    DA_Mpc = 0.0            # angular size distance in units of Mpc
    DA_scale = 0.0          # scale at angular size distance in units of Kpc / arcsec
    DL = 0.0                # luminosity distance in units of c/H0
    DL_Mpc = 0.0            # luminosity distance in units of Mpc
    DMOD = 0.0              # Distance modulus determined from luminosity distance
    a = 0.0                 # 1/(1+z), the scale factor of the Universe

    az = 1.0/(1.0+z)        # 1/(1+z), for the given redshift

    # Compute the integral over a=1/(1+z) from az to 1 in n steps
    n = 1000
    for i in range(n):
        a = az+(1.0-az)*(i+0.5)/n
        adot = math.sqrt(WK+ (WM/a) + (WR/(math.pow(a,2))) +(WV*math.pow(a,2)))
        DCMR = DCMR + 1.0/(a*adot)

    DCMR = (1.0-az)*DCMR/n           # comoving radial distance in units of c/H0
    DCMR_Mpc = (c/H0)*DCMR           # comoving radial distance in units of Mpc

    # Tangental comoving radial distance
    x = math.sqrt(abs(WK))*DCMR
    if x > 0.1:
      if WK > 0.0:
         ratio = 0.5*(math.exp(x)-math.exp(-x))/x
      else:
         ratio = math.sin(x)/x
    else:
      y = math.pow(x,2)
      if WK < 0.0:
         y=-y
      ratio = 1 + y/6.0 + math.pow(y,2)/120.0

    DA = az*ratio*DCMR               #angular size distance in units of c/H0
    DA_Mpc = (c/H0)*DA               #angular size distance in units of Mpc
    DA_scale = DA_Mpc/206.264806     #scale at angular size distance in units of Kpc / arcsec
    DL = DA/math.pow(az,2)                #luminosity distance in units of c/H0
    DL_Mpc = (c/H0)*DL               #luminosity distance in units of Mpc
    DMOD = 5*math.log10(DL_Mpc*1e6)-5     #Distance modulus determined from luminosity distance


    results = \
    {
      "dcmr_mpc": DCMR_Mpc,
      "da_mpc": DA_Mpc,
      "da_scale": DA_scale,
      "dl_mpc": DL_Mpc,
      "dmod": DMOD,
      "z" : z
    }

    return results


## LAST MODIFIED : April 15, 2013
## CREATED : April 15, 2013
## AUTHOR : DRYX
def distance2redshift(
        DL_Mpc):
    """Convert a luminosity distance to a redshift

    **Key Arguments:**
        - ``DL_Mpc`` -- luminosity distance (Mpc)

    **Return:**
        - ``redshift`` -- the calculated redshift
    """
    ################ > IMPORTS ################
    ## STANDARD LIB ##
    import math
    ## THIRD PARTY ##
    ## LOCAL APPLICATION ##

    lowerLimit = 0.
    upperLimit = 1.
    redshift = upperLimit - lowerLimit
    distGuess = float(redshift2distance(redshift)['dl_mpc'])

    distDiff = DL_Mpc - distGuess

    while math.fabs(distDiff) > 0.01:
        if distGuess < DL_Mpc:
            lowerLimit = redshift
            redshift = lowerLimit + (upperLimit - lowerLimit)/2.
            distGuess = float(redshift2distance(redshift)['dl_mpc'])
        elif distGuess > DL_Mpc:
            upperLimit = redshift
            redshift = lowerLimit + (upperLimit - lowerLimit)/2.
            distGuess = float(redshift2distance(redshift)['dl_mpc'])
        distDiff = DL_Mpc - distGuess

    redshift = float("%5.4f" % (redshift,))

    return redshift
    
######

def setTransientIDs(nobjects):

	return ['transient_%s' %(str(x)) for x in range(0, nobjects)]

def setPopulationDistances(nobjects, lightcurve_settings):

	t0 = time.time()

	nevents = lightcurve_settings['population']['number to inject']
	nshells = 75
	
	if lightcurve_settings['depth']['min redshift'] == 0.0:
		min_redshift = 0.000001
	else:
		min_redshift = lightcurve_settings['depth']['min redshift']
	
	max_redshift = lightcurve_settings['depth']['max redshift']

	redshift = np.linspace(min_redshift, max_redshift, nshells)

	redshiftl = redshift[:-1]
	redshiftu = redshift[1:]
	redshiftm = (redshiftl + redshiftu) / 2.

	distancem = [redshift2distance(x)['dl_mpc'] for x in redshiftm]
	distancel = [redshift2distance(x)['dl_mpc'] for x in redshiftl]
	distanceu = [redshift2distance(x)['dl_mpc'] for x in redshiftu]


	volumem = (4. / 3.)* np.pi * (np.array(distancem))**3
	volumel = (4. / 3.)* np.pi * (np.array(distancel))**3
	volumeu = (4. / 3.)* np.pi * (np.array(distanceu))**3
	volume_total = volumeu[-1]

	wredshift = (volumeu - volumel) / volume_total
	indices = np.arange(0, nshells - 1, 1)

	redshifts = []

	for i in range(0, nevents):

		the_index = random.choices(indices, weights = wredshift, k = 1)
		redshifts.append(random.uniform(redshiftl[the_index[0]], redshiftu[the_index[0]]))
	
	distances = [redshift2distance(x)['dl_mpc'] for x in redshifts]	
	
	t1 = time.time()
	print('Distances/redshifts set in %.3f s' %(t1-t0))

	return list(distances), list(redshifts)

def setPopulationExplosionEpochs(nobjects, the_survey, lightcurve_settings):

	nobjects = lightcurve_settings['population']['number to inject']

	survey_begin = the_survey.begin_jd
	survey_end = the_survey.end_jd

	explosion_epochs = [random.uniform( survey_begin, survey_end) for nelements in range(nobjects)]
	
# 	population_settings.update({'explosion epochs': explosion_epochs})

	return list(explosion_epochs)
	
def setPopulationSkyCoords(nobjects, the_survey, lightcurve_settings):

	min_ra = the_survey.lower_ra_bound
	max_ra = the_survey.upper_ra_bound
	min_dec = the_survey.lower_dec_bound
	max_dec = the_survey.upper_dec_bound
	
	ras = [random.uniform( min_ra, max_ra) for nelements in range(nobjects)]
	
# 	print('Declination setting is done with a triangular distribution. Please fix.')
# 	decs = [random.triangular( min_dec, max_dec, 0.0) for nelements in range(nobjects)]

	decs = []

	dec = np.linspace(min_dec, max_dec, 50)

	decl = dec[:-1]
	decu = dec[1:]

	decm = (decl + decu) / 2.
	wdec = np.cos(decm * np.pi/180.)

	indices = np.arange(0, 49, 1)

	for i in range(0, nobjects):

		the_index = random.choices(indices, weights = wdec, k = 1)
	# 	print(the_index)
		decs.append(random.uniform(decl[the_index[0]], decu[the_index[0]]))
	
# 	plt.hist(decs, bins = 50)
# 	plt.show()
	
# 	population_settings.update({'ras': ras, 'decs': decs})
	
	return list(ras), list(decs)

def generatePopulation(the_survey, lightcurve_settings):

	nobjects = lightcurve_settings['population']['number to inject']
	
	transient_ids = setTransientIDs(nobjects)
	distances, redshifts = setPopulationDistances(nobjects, lightcurve_settings)
	explosion_epochs = setPopulationExplosionEpochs(nobjects, the_survey, lightcurve_settings)
	ras, decs = setPopulationSkyCoords(nobjects, the_survey, lightcurve_settings)
	
	return {'transient_id': transient_ids, 'distances': distances, 'redshifts': redshifts, 'explosion epochs': explosion_epochs, 'ras': ras, 'decs': decs}
	
	

