#!/usr/bin/env python

import sys
import math
from math import pi, degrees, radians
import ROOT
from ROOT import TTree
import numpy as np
import numpy.ma as ma
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


def find_galoff_healpxs(nhpside=16, path_catalogue_excl='/nfs/farm/g/glast/u/mtakahas/data/catalogue/gll_psch_v13.fit', dict_flux_cutangle={2.5E-10:4., 5E-9:10.}):
#def find_galoff_healpxs(nhpside=16, threshold_flux50=0, path_catalogue_excl='/nfs/farm/g/glast/u/mtakahas/data/catalogue/gll_psch_v13.fit', dict_flux_cutangle={1E-10:4., 1E-9:10.}):
    """Returns a list of HEALPiX pixles which satisfy your conditions as OFF region of Galactic ridge. The coordinate system is RADEC.
"""
    NPIX = hppf.nside2npix(nhpside)
    MAX_PIX_DEG = np.degrees(hppf.max_pixrad(nhpside))
    ar_pix_incl = []
    pixs = np.array(range(NPIX))
    theta, phi = hppf.pix2ang(nhpside, pixs)
    radec_rad = [phi, np.pi/2.-theta]
    scoord = SkyCoord(radec_rad[0], radec_rad[1], unit="rad")
    qualif_highb = (np.abs(scoord.galactic.b)>=50*u.deg) * (90*u.deg<=scoord.galactic.l) * (scoord.galactic.l<270*u.deg)
    scoord_highb = scoord[qualif_highb]
    pixs_highb = pixs[qualif_highb]

    qualif_multi = qualif_highb

    dict_separation_min = {}
    dict_qualif_angsep = {}

    for flux_threshold, cutangle in dict_flux_cutangle.items():
        print flux_threshold
        dict_separation_min[flux_threshold] = calcAngDist_pixarray_Catalogue(scoord_highb, path_catalogue_excl, threshold=flux_threshold)
        dict_qualif_angsep[flux_threshold] = dict_separation_min[flux_threshold]>(cutangle+np.degrees(hppf.nside2resol(nhpside)))
        print dict_qualif_angsep[flux_threshold]

    qualif_nosrc = np.prod(np.array(dict_qualif_angsep.values()), axis=0)
    print qualif_nosrc
    return pixs_highb[qualif_nosrc>0]


    # for ipix in range(NPIX):
    #     theta, phi = hppf.pix2ang(nhpside, ipix)
    #     radec_rad = [phi, pi/2.-theta] #[pi/2.-theta, phi]
    #     radec_deg = [degrees(radec_rad[0]), degrees(radec_rad[1])]
    #     gal_deg = convert_radec_to_galactic(radec_deg)
    #     if abs(gal_deg[1])>=50 and 90<=gal_deg[0]<270:
    #         ar_ang_dist = calcAngDist_Catalogue(radec_deg, path_catalogue_excl, threshold_flux50)
    #         for di in ar_ang_dist:
    #             if di['ANG_DIST']<=MAX_PIX_DEG:
    #                 break
    #         else:
    #             ar_pix_incl.append(ipix)
    #return ar_pix_incl

    
def find_pointon_healpxs(ra_poi, dec_poi, ang_cut, nhpside=512):
    NPIX = hppf.nside2npix(nhpside)
#    MAX_PIX_DEG = degrees(hppf.max_pixrad(nhpside))
    RADEC_POI = [radians(ra_poi), radians(dec_poi)]
    VEC_POI = hppf.ang2vec(pi/2.-RADEC_POI[1], RADEC_POI[0])
    ar_pix_incl = []
    ang_cut_rad = radians(ang_cut)
    for ipix in range(NPIX):
        vec = hppf.pix2vec(nhpside, ipix)
        ang_dist = get_ang_dist_vectors(VEC_POI, vec)
        if ang_dist<ang_cut_rad:
            ar_pix_incl.append(ipix)
    return ar_pix_incl
    

# def find_pointon_healpxs(ra_poi, dec_poi, ang_cut, nhpside=512):
#     NPIX = hppf.nside2npix(nhpside)
#     coo_poi = SkyCoord(ra_poi, dec_poi, unit="deg")
#     ar_pix_incl = []
#     for ipix in range(NPIX):
#         ang_pix = hppf.pix2ang(nhpside, ipix)
#         coo_pix = SkyCoord(ang_pix[1], pi/2.-ang_pix[0], unit="deg")
#         ang_pix2poi = coo_poi.separation(coo_pix)
#         deg_pix2poi = float(ang_pix2poi.to_string(unit=u.deg, decimal=True))
#         if degrees(deg_pix2poi)<ang_cut:
#             ar_pix_incl.append(ipix)
#     return ar_pix_incl
    
    
    
