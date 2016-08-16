#!/usr/bin/env python
"""Calculate angular distance between an sky position and every source in a catalogue
"""
import sys
from astropy.io import fits
#from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
#from astropy.coordinates import Angle, Latitude, Longitude  # Angles
#import astropy.units as u
from array import array
import math
from math import sin, cos, asin, acos, pi
import numpy as np


def vecterSkyPosition(aPosiRad):
    """Return a 3D vector (NumPy array) corresponding to the given (RA, DEC)
vecterSkyPosition(RA in radians, DEC in radians)
"""
    vecPosi = np.array([cos(aPosiRad[0])*sin(pi/2.-aPosiRad[1]), sin(aPosiRad[0])*sin(pi/2.-aPosiRad[1]), sin(pi/2.-aPosiRad[1])])
    return vecterSkyPosition


def calcAngDist_Catalogue(aPosiDeg, fitsCatalogue):
    """Calculate angular distance between a sky position and evergy source in a catalogue. Return an array of angular distance in radians.
calcAngDist_Catalogue([RA in deg, DEC in deg], 'Path of your catalogue')
"""
    HDULIST = fits.open(fitsCatalogue)
    TBDATA = HDULIST[1].data
    npaPosiDeg = np.array(aPosiDeg)
    A_POSI_RAD = npaPosiDeg * math.radians(1)
    VEC_POSI = vecterSkyPosition(A_POSI_RAD)
    NSRC = len(TBDATA['Source_Name'])
    aAngDist = np.array([])
    for iSrc in range(NSRC):
        vecSrc = vecterSkyPosition([TBDATA['RAJ2000'][iSrc], TBDATA['DECJ2000'][iSrc]])
        angRad = acos(np.dot(VEC_POSI, vecSrc))
        aAngDist.append(angRad)
    return aAngDist
