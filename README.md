![Python](https://img.shields.io/badge/Python->=3.10-blue)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# Binary ecLipsE Prediction Software (BLEPS)
A simple code for getting fast predictions on upcoming or recent eclipses or transits for eclipsing binary stars. 

BLEPS is primarily designed for ground-based observers who wish to observe eclipsing binaries from smaller observatories or with amateur equipment, e.g. to measure the flux ratio in several photometric bands.

BLEPS takes basic information about your targets and calculates all feasible eclipses for your location, considering daylight and altitude limits.
Output options include:
- A .csv file containing all of the info you need to prepare observations
- A .csv file which you can import to Google Calendar
- An astropy table object (if scripting with BLEPS)

#### CAUTION: All times are calculated in UTC. Take care to consider timezones when planning your observations!

## Setting up
BLEPS uses a configuration file to set up the calculations. It should be named `config.yaml` and kept in the same folder as the `bleps.py` script.
In here, you can point to your list of target info, specify your location / altitude / date constraints, and choose how you want to save the output.
You can also specify what observing strategy you want to use to hunt for visible eclipses.

### Observing strategy
Those familiar with observing long-period eclipsing binaries are familiar with the problem of eclipse durations being longer than a single night. 
In this case, if your system has (near)-total eclipses, you can choose to split the eclipse in half, so at least you can get a measure of the flux ratio 
by fixing the geometry of the system with e.g. a prior fit to TESS light curves.
See [Miller et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.5129M/abstract) for an example of this strategy in practice.

There are three options for observing strategy in BLEPS:
- `'split'`: Consider each half of an eclipse separately in feasibility checks (ingress/egress).
- `'full'`: The whole eclipse must be pass feasibility checks in order to be returned.
- `'any'`: Returns info if the system is observable at any point during the eclipse.

### Target file
BLEPS takes basic information on your target(s) from a .csv file which *must* have the following named columns:

- 'Name' - str - name of system(s)
- 'RA' - str - right ascension of target in hour angle format 'XX XX XX.XX'
- 'Dec' - str - declination of target in degree format '+/-XX XX XX.X'
- 'T0' - float - time of primary eclipse minimum in HJD or BJD
- 'Period' - float - orbital period of system in days
- 'Width1' - float - width of primary eclipse in phase units
- 'Phase2' - float - phase of secondary eclipse in phase units
- 'Width2' - float - width of secondary eclipse in phase units

You can specify the path to this file in the `config.yaml` file.

## Google Calendar integration
If you run BLEPS with the `gcal` option set to `True`, then you should generate a .csv file that can be imported into Google Calendar. 
I suggest making a new calendar for this, setting the timezone to UTC to prevent unexpected time conversions. 
You can import a calendar from the Settings menu - [see this guide from Google for more info](https://support.google.com/calendar/answer/37118?hl=en&co=GENIE.Platform%3DDesktop).
