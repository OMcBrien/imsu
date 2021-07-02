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

def setPopulationDistances(survey_settings, lightcurve_settings):

	t0 = time.time()

	nobjects = lightcurve_settings['population']['number to inject']

	min_redshift = np.array(survey_settings['depth']['min redshift'])
	max_redshift = np.array(survey_settings['depth']['max redshift'])
	
	if min_redshift == 0.0:
		
		min_distance = 0.0
		min_volume_enclosed = 0.0
		
	else:
	
		min_distance = redshift2distance(min_redshift)['dl_mpc']
		min_volume_enclosed = (4./3.) * np.pi * min_distance**3
	
	max_distance = redshift2distance(max_redshift)['dl_mpc']
	max_volume_enclosed = (4./3.) * np.pi * max_distance**3
	
	volumes_enclosed = [random.uniform( min_volume_enclosed, max_volume_enclosed) for nelements in range(nobjects)]
	
	distances = ((3./4.) * (1./np.pi) * np.array(volumes_enclosed))**(1./3.)
	
# 	plt.hist(distances, bins = int(np.sqrt(nobjects)))
# 	plt.xlabel('Distance (Mpc)')
# 	plt.show()
	
	redshifts = []
	
	t1 = time.time()
	print('Distance distribution performed in %.4f s' %(t1 - t0))
	
	for distance in distances:
	
		redshifts.append(distance2redshift(distance))
	
	t2 = time.time()
	print('Distance to redshift conversion in %.4f s' %(t2 - t1))
	
# 	plt.hist(redshifts, bins = int(np.sqrt(nobjects)))
# 	plt.xlabel('Redshifts, $z$')
# 	plt.show()
	
	transient_id = ['transient_%s' %(str(x)) for x in range(0, nobjects)]
	
	population_settings = {'transient_id': transient_id, 'distances': list(distances), 'redshifts': redshifts}

	return population_settings

def setPopulationExplosionEpochs(survey_settings, lightcurve_settings, population_settings):

	nobjects = lightcurve_settings['population']['number to inject']

	survey_begin = survey_settings['time observing']['survey begin']
	survey_end = survey_settings['time observing']['survey end']

	explosion_epochs = [random.uniform( survey_begin, survey_end) for nelements in range(nobjects)]
	
	population_settings.update({'explosion epochs': explosion_epochs})

	return population_settings
	
def setPopulationSkyCoords(survey_settings, lightcurve_settings, population_settings):

	nobjects = lightcurve_settings['population']['number to inject']

	min_ra = survey_settings['sky visible']['min ra']
	max_ra = survey_settings['sky visible']['max ra']
	min_dec = survey_settings['sky visible']['min dec']
	max_dec = survey_settings['sky visible']['max dec']
	
	ras = [random.uniform( min_ra, max_ra) for nelements in range(nobjects)]
	
# 	plt.hist(ras, bins = int(np.sqrt(nobjects)))
# 	plt.xlabel('RA (deg.)')
# 	plt.show()
	
	print('Declination setting is done with a triangular distribution. Please fix.')
	decs = [random.triangular( min_dec, max_dec, 0.0) for nelements in range(nobjects)]
	
# 	plt.hist(decs, bins = int(np.sqrt(nobjects)))
# 	plt.xlabel('Dec (deg.)')
# 	plt.show()
	
	population_settings.update({'ras': ras, 'decs': decs})
	
	return population_settings

