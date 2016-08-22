#!/usr/bin/env python
"""Calculate angular distance between an sky position and every source in a catalogue
"""
import sys
from astropy.io import fits
from array import array
import math
from math import sin, cos, asin, acos, pi
import numpy as np


def suppress_vector_too_small_compo(vecPosi):
    vecSupp = vecPosi.copy()
    for iel in range(len(vecPosi)):
        if abs(vecPosi[iel])<1e-15:
            vecSupp[iel] = 0.
    return vecSupp


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
    vecRslt = suppress_vector_too_small_compo(vecPosi)
    return vecRslt


def calcAngDist_Catalogue(aPosiDeg, fitsCatalogue, threshold=0, ang_print=0, nSrc=0):
    """Calculate angular distance between a sky position and evergy source in a catalogue. Return an array of dictionaries which have "Source_Name", "ASSOC", "ANG_DIST" in degrees, "Flux50".
calcAngDist_Catalogue([RA in deg, DEC in deg], 'Path of your catalogue')
"""
    HDULIST = fits.open(fitsCatalogue)
    TBDATA = HDULIST[1].data
    npaPosiDeg = np.array(aPosiDeg)
    VEC_POSI = vecterSkyPosition(aPosiDeg)
    #print "Given sky position vector:", VEC_POSI
    aDictResults = []
    if nSrc==0:
        nSrc = len(TBDATA['Source_Name'])
    for iSrc in range(nSrc):
        vecSrc = vecterSkyPosition([TBDATA['RAJ2000'][iSrc], TBDATA['DEJ2000'][iSrc]])
        angRad = get_ang_dist_vectors(VEC_POSI, vecSrc)
        if TBDATA['Flux50'][iSrc]>=threshold:
            aDictResults.append({"Source_Name":TBDATA['Source_Name'][iSrc], "ASSOC":TBDATA['ASSOC'][iSrc], "ANG_DIST":math.degrees(angRad), "Flux50":TBDATA['Flux50'][iSrc]})
            if math.degrees(angRad)<ang_print:
                print aDictResults[-1]
    return aDictResults


def get_ang_dist_vectors(vec0, vec1):
    """Return the angular distance between two vectors (should be normalized) in radians.
"""
    prodDot = np.dot(vec0, vec1)
    return acos(prodDot)
    
