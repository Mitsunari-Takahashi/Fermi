#!/usr/bin/env python

import sys
from astropy.io import fits
from astropy.time import Time
#from astropy.time import TimeDelta
#from astropy.coordinates import SkyCoord  # High-level coordinates
#from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
#from astropy.coordinates import Angle, Latitude, Longitude  # Angles
#import astropy.units as u
from skyfield.api import load
from skyfield.api import utc
from array import array
#import math
#import numpy as np
#import datetime
from pMETandMJD import *

par = sys.argv
print par
metTI = float(par[1])

planets = load('de421.bsp')
earth, jupiter = planets['EARTH'], planets['JUPITER BARYCENTER']
ts = load.timescale()
mjdTI = ConvertMetToMjd(metTI)
tTI = Time(mjdTI, format='mjd')
utcTI = ts.utc(tTI.to_datetime(timezone=utc))
amJupiter = earth.at(utcTI).observe(jupiter)
raJ, decJ, distJ = amJupiter.radec()
degRaJ, degDecJ = raJ._degrees, decJ._degrees
bJ, lJ, distGJ = amJupiter.galactic_latlon()
degLJ, degBJ = lJ._degrees, bJ._degrees
#radRaJ, radDecJ = math.radians(degRaJ), math.radians(degDecJ)
aOutCel = [degRaJ, degDecJ, distJ]
print "CEL:", aOutCel
aOutGal = [degLJ, degBJ, distGJ]
print "GAL:", aOutGal
