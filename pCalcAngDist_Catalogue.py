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


def vecterSkyPosition(aPosiDeg):
    """Return a 3D vector (NumPy array) corresponding to the given (RA, DEC)
vecterSkyPosition(RA in radians, DEC in degrees)
"""
    npaPosiDeg = np.array(aPosiDeg)
    aPosiRad = npaPosiDeg * math.radians(1)
    cosPhi = cos(aPosiRad[0])
    cosTh = cos(pi/2.-aPosiRad[1])
    sinPhi = sin(aPosiRad[0])
    sinTh = sin(pi/2.-aPosiRad[1])
    vecPosi = np.array([cosPhi*sinTh, sinPhi*sinTh, cosTh])
    for iel in range(len(vecPosi)):
        if abs(vecPosi[iel])<1e-15:
            vecPosi[iel] = 0.
    return vecPosi


def calcAngDist_Catalogue(aPosiDeg, fitsCatalogue, nSrc=0):
    """Calculate angular distance between a sky position and evergy source in a catalogue. Return an array of angular distance in radians.
calcAngDist_Catalogue([RA in deg, DEC in deg], 'Path of your catalogue')
"""
    HDULIST = fits.open(fitsCatalogue)
    TBDATA = HDULIST[1].data
    npaPosiDeg = np.array(aPosiDeg)
    VEC_POSI = vecterSkyPosition(aPosiDeg)
    print "Given sky position vector:", VEC_POSI
    if nSrc==0:
        nSrc = len(TBDATA['Source_Name'])
    aAngDist = np.array([])
    for iSrc in range(nSrc):
        vecSrc = vecterSkyPosition([TBDATA['RAJ2000'][iSrc], TBDATA['DEJ2000'][iSrc]])
        prodDot = np.dot(VEC_POSI, vecSrc)
        angRad = acos(prodDot)
        aAngDist = np.append(aAngDist, angRad)
        if math.degrees(angRad)<15:
            print TBDATA['Source_Name'][iSrc], "(", TBDATA['ASSOC'][iSrc], ") :", math.degrees(angRad), "deg, Flux50=", TBDATA['Flux50'][iSrc]
    return [TBDATA['Source_Name'][:nSrc], aAngDist]
