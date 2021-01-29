#!/usr/bin/env python

import sys
import ROOT
#from ROOT import TTree
#from ROOT import TChain
import numpy as np
import datetime
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
#import healpy as hp
from healpy import pixelfunc as hppf
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pMETandMJD import *
from pAnalysisConfig import *
from pFindHEALPix import *
import commands


def make_livetime_histogram(aHtgLt, nRegion, pathFileScAll, metStart, metStop, aFileToI, aCoordsPix_array, aAreaPix_array, origin_time=0, metExStart=0, metExStop=0):
    """Look over spacecraft files and make a histogram of (solid angle [sr] * time interval [sec]) on MET vs. Zenith angle vs. Cos(Inclination)"""
    fmwStart = ConvertMetToFMW(metStart)
    fmwStop = ConvertMetToFMW(metStop)
    print '=========================================='
    print 'Going to look MET', metStart, '-', metStop
    print "Fermi Mission Week:", fmwStart, "-", fmwStop
    print '=========================================='
    print 'MET',metExStart, '-', metExStop, 'is excluded.'
    cmd = "ls {0}".format(pathFileScAll)
    ret = commands.getoutput(cmd)
    aPathFileScAll = ret.split("\n")
    for iFileSc in range(len(aPathFileScAll)):
        strPathFileSc = aPathFileScAll[iFileSc]
        #fmwFile = int(strPathFileSc[-27:-24])
        # Make file list
        #if fmwFile>=int(fmwStart-1) and fmwFile<=int(fmwStop+1) :
        aFileToI.append(aPathFileScAll[iFileSc])
    timeStart = datetime.datetime.now() # For progress bar
    for fileToI in aFileToI:
        hdulistSC = fits.open(fileToI)
        tbdataSC = hdulistSC[1].data
        aSTART, aSTOP = tbdataSC.field('START'), tbdataSC.field('STOP')
        aRA_ZENITH = tbdataSC.field('RA_ZENITH')
        aDEC_ZENITH = tbdataSC.field('DEC_ZENITH')
        aRA_SCZ = tbdataSC.field('RA_SCZ')
        aRA_SCX = tbdataSC.field('RA_SCX')
        aDEC_SCZ = tbdataSC.field('DEC_SCZ')
        aDEC_SCX = tbdataSC.field('DEC_SCX')
        aLIVETIME = tbdataSC.field('LIVETIME')
        aDATA_QUAL = tbdataSC.field('DATA_QUAL')
        aLAT_CONFIG = tbdataSC.field('LAT_CONFIG')
        nTI = len(aSTART)
        print "  ", fileToI, "(", nTI, "intervals in total)"

        dict_cut = {}
        dict_cut['time'] = np.array(( (aSTART>=metStart) * (aSTART<metStop) ) * ( ((metExStart==0) * (metExStop==0)) + (aSTOP<metExStart) + (aSTART>=metExStop)), dtype=bool)
        dict_cut['qual'] = aDATA_QUAL>0
        dict_cut['config'] = aLAT_CONFIG==1
        print '* Cuts'
        for k, v in dict_cut.items():
            print '{k}: {v}'.format(k=k, v=v[:100])
        cut_combined_interval =  np.array(dict_cut['time'] * dict_cut['qual'] * dict_cut['config'], dtype=bool) #np.prod(dict_cut.values())
        print '* Combined cut'
        print cut_combined_interval[:100]

        # coordsSCZ_interval = SkyCoord(aRA_SCZ, aDEC_SCZ, unit="deg")
        # coordsZenith_interval = SkyCoord(aRA_ZENITH, aDEC_ZENITH, unit="deg")
        # tplot_interval = aSTART-origin_time
        #vecSCX = np.array(hppf.ang2vec(math.pi/2.-math.radians(aDEC_SCX), math.radians(aRA_SCX)))
        #vecSCZ = np.array(hppf.ang2vec(math.pi/2.-math.radians(aDEC_SCZ), math.radians(aRA_SCZ)))

        for iR in range(nRegion):
            if (isinstance(aAreaPix_array[iR], int) or isinstance(aAreaPix_array[iR], float)) and aAreaPix_array[iR]==0: # Point-like region
                angSCZ = coordsSCZ.separation(aCoordsPix_array[iR])
                radSCZ = float(angSCZ.to_string(unit=u.rad, decimal=True))
                angZenith = coordsZenith.separation(aCoordsPix_array[iR])
                degZenith = float(angZenith.to_string(unit=u.deg, decimal=True))
                aHtgLt[iR].Fill(cos(radSCZ), degZenith, tplot, tti)
            else:  # Extended region
                npixel = len(aAreaPix_array[iR])
                coordsSCZ_interval = SkyCoord(np.array([aRA_SCZ]).T, np.array([aDEC_SCZ]).T, unit="deg")
                print 'Telescope boresight has been reconstructed.'
                coordsZenith_interval = SkyCoord(np.array([aRA_ZENITH]).T, np.array([aDEC_ZENITH]).T, unit="deg")
                print 'Zenith direction has been reconstructed.'
                #coordsPix_pixel = aCoordsPix_array[iR]

                angSCZ_interval_pixel = coordsSCZ_interval.separation(aCoordsPix_array[iR])
                cosSCZ_interval_pixel = np.cos(angSCZ_interval_pixel)
                print 'Inclination for each pixel has been calculated.'

                angZenith_interval_pixel = coordsZenith_interval.separation(aCoordsPix_array[iR])
                #degZenith = float(angZenith.to_string(unit=u.deg, decimal=True))
                print 'Zenith angle for each pixel has been calculated.'

                tplot_interval = aSTART-origin_time

                time_valid = cut_combined_interval * aLIVETIME #* aAreaPix_array[iR]
                print time_valid[:100]
                livetime_interval_pixel = np.dot(np.array([time_valid]).T, np.array([aAreaPix_array[iR]]))
                print 'Livetime for each pixel has been calculated.'

                for iTI in range(nTI):
                    for jpix in range(len(aCoordsPix_array[iR])):
#                         if iTI<10 and livetime_interval_pixel[iTI][jpix]>0:
#                             print '''Inclination: {cth}
# Zenith: {za}
# Time: {tp}
# Livetime: {lv}
# '''.format(cth=cosSCZ_interval_pixel[iTI][jpix], za=float(angZenith_interval_pixel[iTI][jpix].to_string(unit=u.deg, decimal=True)), tp=tplot_interval[iTI], lv=livetime_interval_pixel[iTI][jpix])
                        aHtgLt[iR].Fill(cosSCZ_interval_pixel[iTI][jpix], float(angZenith_interval_pixel[iTI][jpix].to_string(unit=u.deg, decimal=True)), tplot_interval[iTI], livetime_interval_pixel[iTI][jpix])
                    if iTI%100==0:
                        print '{0:2.0f}% have been filled.'.format(iTI*100/nTI)
                        sys.stdout.flush()

        # nTI_included = 0
        # for iTI in range(nTI):
        #     if ( aSTART[iTI]>=metStart and aSTART[iTI]<metStop ) and ( (metExStart==0 and metExStop==0) or aSTOP[iTI]<metExStart or aSTART[iTI]>=metExStop):
        #         nTI_included = nTI_included+1
        #         if not aDATA_QUAL[iTI]>0:
        #             print 'Bad time interval', aSTART[iTI], '-', aSTOP[iTI], ':', aDATA_QUAL[iTI]
        #             continue
        #         if not aLAT_CONFIG[iTI]==1:
        #             print 'LAT config:', aSTART[iTI], '-', aSTOP[iTI], ':', aLAT_CONFIG[iTI]
        #             continue
        #         tti = aLIVETIME[iTI]
        #         coordsSCZ = SkyCoord(aRA_SCZ[iTI], aDEC_SCZ[iTI], unit="deg")
        #         coordsZenith = SkyCoord(aRA_ZENITH[iTI], aDEC_ZENITH[iTI], unit="deg")
        #         tplot = aSTART[iTI]-origin_time
        #         vecSCX = np.array(hppf.ang2vec(math.pi/2.-math.radians(aDEC_SCX[iTI]), math.radians(aRA_SCX[iTI])))
        #         vecSCZ = np.array(hppf.ang2vec(math.pi/2.-math.radians(aDEC_SCZ[iTI]), math.radians(aRA_SCZ[iTI])))
        #         for iR in range(nRegion):
        #             if aAreaPix_array[iR]==0:
        #                 angSCZ = coordsSCZ.separation(aCoordsPix_array[iR])
        #                 radSCZ = float(angSCZ.to_string(unit=u.rad, decimal=True))
        #                 angZenith = coordsZenith.separation(aCoordsPix_array[iR])
        #                 degZenith = float(angZenith.to_string(unit=u.deg, decimal=True))
        #                 aHtgLt[iR].Fill(cos(radSCZ), degZenith, tplot, tti)
        #             else:                        
        #                 for (jpix, coordsPix) in enumerate(aCoordsPix_array[iR]):
        #                     angSCZ = coordsSCZ.separation(coordsPix)
        #                     radSCZ = float(angSCZ.to_string(unit=u.rad, decimal=True))
        #                     angZenith = coordsZenith.separation(coordsPix)
        #                     degZenith = float(angZenith.to_string(unit=u.deg, decimal=True))
        #                     aHtgLt[iR].Fill(cos(radSCZ), degZenith, tplot, tti*aAreaPix_array[iR][jpix])
        #         if iTI%20==0:
        #             print iTI, aSTART[iTI], aRA_SCZ[iTI], aDEC_SCZ[iTI], tbdataSC.field('LAT_MODE')[iTI], aDATA_QUAL[iTI]
#                sys.stdout.flush()
        print 'Histogram filling has been finished.' #print nTI_included, 'intervals are included.'
