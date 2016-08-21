#!/usr/bin/env python

import sys
import math
from math import pi, degrees, radians
import ROOT
from ROOT import TTree
import numpy as np
from array import array
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
import healpy as hp
from healpy import pixelfunc as hppf
from pConvertSkyCoords import *
from pCalcAngDist_Catalogue import *


def find_galoff_healpxs(nhpside=16, threshold_flux50=0, path_catalogue_excl='/nfs/farm/g/glast/u/mtakahas/data/catalogue/gll_psch_v09.fit'):
    """Returns a list of HEALPiX pixles which satisfy your conditions as OFF region of Galactic ridge. The coordinate system is RADEC.
"""
    NPIX = hppf.nside2npix(nhpside)
    MAX_PIX_DEG = degrees(hppf.max_pixrad(nhpside))
    ar_pix_incl = []
    for ipix in range(NPIX):
        theta, phi = hppf.pix2ang(nhpside, ipix)
        radec_rad = [phi, pi/2.-theta] #[pi/2.-theta, phi]
        radec_deg = [degrees(radec_rad[0]), degrees(radec_rad[1])]
        gal_deg = convert_radec_to_galactic(radec_deg)
        if abs(gal_deg[1])>=50 and 90<=gal_deg[0]<270:
            ar_ang_dist = calcAngDist_Catalogue(radec_deg, path_catalogue_excl, threshold_flux50)
            for di in ar_ang_dist:
                if di['ANG_DIST']<=MAX_PIX_DEG:
                    break
            else:
                ar_pix_incl.append(ipix)
    return ar_pix_incl

    
def find_pointon_healpxs(ra_poi, dec_poi, ang_cut, nhpside=512):
    NPIX = hppf.nside2npix(nhpside)
    MAX_PIX_DEG = degrees(hppf.max_pixrad(nhpside))
    RADEC_POI = [radians(ra_poi), radians(dec_poi)]
#    print RADEC_POI
    VEC_POI = hppf.ang2vec(pi/2.-RADEC_POI[1], RADEC_POI[0])
    ar_pix_incl = []
    for ipix in range(NPIX):
        vec = hppf.pix2vec(nhpside, ipix)
        ang_dist = get_ang_dist_vectors(VEC_POI, vec)
        if degrees(ang_dist)<ang_cut:
            ar_pix_incl.append(ipix)
    return ar_pix_incl
    

    
    
    
