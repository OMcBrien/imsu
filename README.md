# imsu

This is (hopefully) a generalised efficiency simulation code for astronomical surveys to use. It takes a list of pointings made by a telescope and returns a measurement of the efficiency with which that telescope can detect a particular astrophysical transients, like supernovae and kilonovae. I developed a [prototype version](https://github.com/OMcBrien/imitation-survey) of this for my PhD, specifically with the ATLAS project in mind ([Tonry et al., (2018)](https://ui.adsabs.harvard.edu/abs/2018PASP..130f4505T/abstract)), but rewrote it to be compatible with other robotic surveys, like GOTO ([Steeghs et al., (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv211005539S/abstract)). The principle of operation of both codes is the same, but this version is a good deal faster and has been extended to allow for GRB lightcurve injection and kilonova lightcurve synthesis.

A note on the naming of this cose: I called it **imsu** as it is short for *imitation survey*, but at this moment in time the code doesn't provide proper imitation survey capabilities. It is instead, like the PhD prototype, an efficiency calculator for different transient phenomena. 

## Basis of operation

The code is divided into a few subroutines that handle its various stages of execution. The only user input is supposed to be through the settings file, `<settings.yaml>`. This documentation will explain in detail how the code runs and how to interpret the results. To give an overview, the code will generate a population of transient objects of the user's choosing, distribute them through a volume of space (including across the sky and over redshift/physical distance) and over a period of time considered by the telescope under consideration, given the exposure history provided. For each transient, a lightcurve is generated using the corresponding exposure times of the frames that contain the transient which, from the limiting magnitudes of the exposure, enables it to identify which transients are formally recovered by the survey. 

## Installation

Unfortunately, I've not packaged this yet such that it can be installed via pip. For now, simply clone the repo somewhere on your machine. change directory into `imsu/` and run the main script with `python runsim.py`. You'll need a few dependencies for it to operate fully. The 'standard' packages to install are `NumPy`, `pandas`, `matplotlib` and `json` if you don't already have them. Beyond that, please install the following:

1. `pyyaml`, which can be pip installed or obtained from [here](https://pypi.org/project/PyYAML/).
2. `SNCosmo`, which is used for lightcurve generation of supernovae and kilonovae-type objects. Read the docs [here](https://sncosmo.readthedocs.io/en/stable/). 
3. `AfterglowPy`, which is used for generating GRB lightcurves. There isn't much documentation for this package, but you can read some example on it [here](https://afterglowpy.readthedocs.io/en/latest/).
4. `sfdmap`, which is used with `SNCosmo` to apply extinction corrections to generated lightcurves. It is pip installable or can be installed from source [here](https://github.com/kbarbary/sfdmap).
5. `extinction`, another extinction correction package this time used with `AfterglowPy`. It's also pip installable, but if you want the documentation, look [here](https://extinction.readthedocs.io/en/latest/).

With the exception of `sfdmap`, these packages should all work out of the box. `sfdmap` has an additional step in its install of downloading the [Schlegel, Finkbeiner and Davis (1998)](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract) dustmap files themselves. These are then queried during the lightcurve generation for E(B-V) values. The dustmaps are pretty big once untarred and are contained in folder named `sfddata-master`. Remember the path to this folder as you'll need it in a moment.

Lastly, if you're affiliated with the GOTO project, you can also install [`gotocat`](https://github.com/GOTO-OBS/gotocat) and query the GOTO database directly for exposure histories. Otherwise, you're own exposure history will need to be provided.

## Contents of code

When you change into the main `imsu/` directory, you can see its contents include the Python code needed for execution and three additional subdirectories. The first is `bandpasses/` and it contains text files with the transmission functions of most of the major filters used by optical telescopes. It also has transmission functions for some of the non-standard filter systems, like ATLAS' *co* bands and the GOTO Baader filters *LRGB*. The second subdirectory is `custom_sources/` and this is where the spectra for unique objects are held to make spectral time series that can be used with `SNCosmo`. At the moment, only AT2017gfo is included, but if there are any additional objects that need added, email me and I'll try to add them. The third is `query_files` and its purpose will be explained in the next section.

In terms of code, the efficiency calculator is split into several modules to handle different parts of execution. To provide an idea of how the code works, I'll provide a walkthrough of its operation.

### Specifying the settings in `settings.yaml`

Of the two YAML files included, `settings.yaml` is the most important. It governs the program controls and input parameters and is intended to be the main way users interact with the efficiency simulation. At the top of the settings file are a few boolean toggles that read:

```
Program Settings:
    Query Database: False
    Generate Transients: False
    Generate Lightcurves: False
    Perform Recovery: False
    Make Plots: False
```

These are used to turn on the various stages of execution. The code has been written in such a way that you can enable these one by one (and actually disable the previous sections) without needed to re-run large sections of code. This can be useful if you generate a population of transients and want to see how different lightcurves behave with those population settings. Note that for GOTO users, the `Query Database` option uses `gotocat` to query the pointing databases using the information provideded under `Query Settings`. You can select all pointings associated with a particular `dbevent`, if you're interested in following a specific GRB for example, or all pointings between two dates, specified under `time observing`. Your query results will be saved to a file named specified by `query file` in the directory `query_files/` and this same file is read if you disabled the `Query Database` option. **For all other users, please supply your own file with a list of exposures to the `query_files/` folder and pass its name to `query file`. You don't need to include the `query_files/` path or the extension, but this should be `.csv`.**

When you first run the code, a main `results/` directory will be created, if it does not already exist, along with a subdirectory unique to your execution. You can set the name of this subdirectory under `IO Settings > results directory`. At this point, all results files excluding the query file will be saved to that subdurectory.

### Generating the transient population with `population_generator.py`

Setting `Program Settings > Generate Transients` enables the population generator. The population generator sets the sky coordinates, distances (as both physical distances and redshifts) and explosion epochs of your transients. This is a separate step to the lightcurve generation, so you can reuse the same transient population with a different transient type when you enable the lightcurve generator. That said, the settings needed for the population generator can be specified from the `Lightcurve Settings`, in addition to the settings needed for the lightcurve generator. The first parameters to set are the `Lightcurve Settings > population` values.

```
Lightcurve Settings:
    population:
        transient type: SN Ia \ SN Ib/c \ SN II-L \ SN II-P \ SN IIn \ Kilonova \ GRB
        number to inject: 1000
        population file: 'name_of_population_file'
```

The transient type must take a specific value. For the standard supernova types, use: 

1. `SN Ia` for Type Ia
2. `SN Ib/c` for Type Ib/c
3. `SN II-L` for Type II-L
4. `SN II-P` for Type II-P
5. `SN IIn` for Type IIn

Additionally, you may set this to `Kilonova` or `GRB` to use the custom AT2017gfo source or switch the lightcurve generator to `AfterglowPy`.

The population generator attempts to distribute events uniformly in time and isotropically over a volume of space that is partly defined by the exposures contained in the query file. For a given list of pointings, assuming a non-GRB transient type injection, the region of sky that transients are distributed over is set by drawing a box from the minimum RA and Dec values to the maximum RA and Dec values. Within that, values of RA are set assuming a uniform distribution and values of Dec are weighted by a cos(Dec) function. If a GRB transient type is used, you specify the geometry of the area in which the transients are distributed. Under `Lightcurve Settings > ellipse geometry` you can specify an elliptical region of sky to place events. It requies the central RA and Dec and values for the semi-major and smei-minor axes. Additionally, you can supply a rotation angle (measued clockwise from the horizontal on a cartesian axis). To distribute events over distance, the code weights their placement such that the number of events located out to a redshift z is propertional to z cubed. You can specify the minimum and maximum redshifts for this distribution under `Lightcurve Settings > depth`.

Once these values are set, the population file containing this data is saved to the specified results directory to a `.csv` file with the name passed to `Lightcurve Settings > population > population file`. If a populationfile already exists in results directory, such as from a previous execution, you can disable the population generator and the file specified here will be read instead.

### Lightcurve generation with `lc_generator.py`

Enabling `Program Settings > Generate Lightcurves` will begin synthesising lightcurves of the transient type sepcified under `Lightcurve Settings > population > transient type`. Note that you must have a population file made for this to work. Two different packages are used for lightcurve generation - `SNCosmo` for optical transients and `AfterglowPy` for GRBs. Both routines work by collecting the exposures from the query file that are temporally coincident with the specific transient under consideration and, from them, selecting the exposures that would contain the transient given the geometry of the detector and the transients location on the sky. At this point, the code determines how bright the transient would be at these epochs using the lightcurve generator and, if they are brighter than the limiting magnitude of the corresponding exposures, they are saved as 'recovered' exposures. These recovered exposures are passed on to the recovery suite at the next stage of execution. For more information on how the lightcurve generators work, see below.

To collect the temporally coincident and spatially overlapping exposures from the exposure list, the code needs to know which columns to use for the exposure times, central RA and Dec, filters and which filters are used (i.e. while `SNCosmo` and `bandpasses/` might recognise Sloan *g* as `sdssg`, for a real exposure this information is probably only saved as `g`). These can be specified under `Lightcurve Settings > telescope properties > column headings`. Likewise, the code also needs to know the dimensions of the chip to identify whether a given exposure contains a transient on the sky. The horizontal and vertical dimensions, in degrees, are set under `Lightcurve Settings > telescope properties > chip width` and `chip height`.

#### SNCosmo

`SNCosmo` is a useful Python package that allows you to synthesise lightcurves of transient phenomena by performing synthetic photometry on a spectral time series of a given transient type. The in-built sources include standard supernova types as well the the SALT-II models, though these are not time series sources so are incompatible with the efficiency calculator. The synthetic photometry capabilities allow you to extract magnitudes in most of the standard filter systems (e.g. SDSS, DES, GUNN etc.) or add your own filter profiles. Additionally, if you have a number of spectra of a particular transient, you can add them as their own time series source. I have done with, under `custom_sources/AT2017gfo` for the first kilonova, AT2017gfo, and can add more upon request.

The filters `SNCosmo` supports by default are viewable in the [documentation](https://sncosmo.readthedocs.io/en/stable/bandpass-list.html). You can select which filters to use by specifying them under `Lightcurve Settings > telescope properties > filter effective wavelengths`. If using in-built filters (e.g. SDSS *griz*), you only need to set this like so:

```
    telescope properties:
        filter effective wavelengths:
            sdssg: 
            sdssr: 
            sdssi: 
            sdssz:
```

The effective wavelengths of these filters are known by SNCosmo already. If using a non-supported filter, check to see if it's available under `bandpasses/`. You must also supply the effective wavelengths, like so:

```
    telescope properties:
        filter effective wavelengths:
            gotoL: 5396.65
            gotoR: 6405.05
            gotoG: 5355.32
            gotoB: 4600.58
```

The units for effective wavelength are Angstrom. The generated lightcurves are saved in a JSON file with the name passed to `Lightcurve Setting > lightcurve file`. No extension is needed here. As before, if you already have a lightcurve file generated, you can disable the lightcurve generator and it will read this file instead for the next stage of execution.

#### AfterglowPy

`AfterglowPy` is another Python package used to generate lightcurves, but this time for GRB afterglows. It can return lightcurves across the electromagnetic spectrum, including the optical, but requires a more delicate selection of parameters. For more information, consult [Ryan et al., (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...896..166R/abstract). These parameters are set under `GRB Settings` and will be utilised by the GRB lightcurve generator if `Lightcurve Settings > population > transient type` is set to `GRB`. The default parameters are for GRB170817A, the gamma accompaniment to AT2017gfo.

Bear in mind that to produce accurate optical lightcurves, the code has to integrate across the width of the transmission function for each filter requested. This can take some time and limits the use cases to filters available under `bandpasses/`.

### Recovery of transients with `recovery.py`

The recovery suit is rather simplistic. Because the lightcurve generator is able to flag which transients are not only visible in each exposure but also which are bright enough to be detected, the recovery suite is used to provide constrains on when a survey would formally consider the events recovered. If a given transient is only viewed above the detection limit in one exposures, it might not be pursued for follow-up by the survey. Hence, within the recovery suite, you may specify the total number of detections that must be made in each filter, as well as on how many nights detections must be made in each filter, for a transient even to count toward the efficiency measurement.

Enable the recovery suite under `Program Settings > Perform Recovery` and refer to the `Recovery Settings`. The `min detection counts` lets you specify that minimum number of total detections that must be made in each filter. Please set the filters to the same as under `Lightcurve Settings > telescope properties > filter effective wavelengths`. Likewise, the minimum number of nights an event must be detected on in each filter can be specified under `min nightly detection counts`. If a given transient passes these specified criteria, the lightcurve file is update to reflect that when returning the final efficiency calculation. This efficiency is the ratio of the number of events generated that pass these criteria to the total number of events generated.
