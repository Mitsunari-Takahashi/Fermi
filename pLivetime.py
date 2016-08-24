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


def make_livetime_histogram(aHtgLt, nRegion, pathFileScAll, metStart, metStop, aFileToI, aCoordsPix_array, aAreaPix_array, origin_time=0):
    """Look over spacecraft files and make a histogram of (solid angle [sr] * time interval [sec]) on MET vs. Zenith angle vs. Cos(Inclination)"""
    fmwStart = ConvertMetToFMW(metStart)
    fmwStop = ConvertMetToFMW(metStop)
    
    print "Fermi Mission Week:", fmwStart, "-", fmwStop

    cmd = "ls {0}".format(pathFileScAll)
    ret = commands.getoutput(cmd)
    aPathFileScAll = ret.split("\n")
    for iFileSc in range(len(aPathFileScAll)):
        strPathFileSc = aPathFileScAll[iFileSc]
        fmwFile = int(strPathFileSc[-27:-24])
        # Make file list
        if fmwFile>=int(fmwStart-1) and fmwFile<=int(fmwStop+1) :
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
        nTI = len(aSTART)
        print "  ", fileToI, "(", nTI, "intervals )"
        for iTI in range(nTI):
            if aSTART[iTI]>=metStart and aSTOP[iTI]<metStop:
                tti = aLIVETIME[iTI]
                coordsSCZ = SkyCoord(aRA_SCZ[iTI], aDEC_SCZ[iTI], unit="deg")
                coordsZenith = SkyCoord(aRA_ZENITH[iTI], aDEC_ZENITH[iTI], unit="deg")
                tplot = (aSTART[iTI]+aSTOP[iTI])/2.-origin_time
                vecSCX = np.array(hppf.ang2vec(math.pi/2.-math.radians(aDEC_SCX[iTI]), math.radians(aRA_SCX[iTI])))
                vecSCZ = np.array(hppf.ang2vec(math.pi/2.-math.radians(aDEC_SCZ[iTI]), math.radians(aRA_SCZ[iTI])))
                vecSCY = np.cross(vecSCZ, vecSCX)
                aAngSCY = hppf.vec2ang(vecSCY)
                coordsSCY = SkyCoord(math.degrees(aAngSCY[1]), math.degrees(math.pi/2.-aAngSCY[0]), unit="deg")
                for iR in range(nRegion):
                    for (jpix, coordsPix) in enumerate(aCoordsPix_array[iR]):
                        angSCZ = coordsSCZ.separation(coordsPix)
                        radSCZ = float(angSCZ.to_string(unit=u.rad, decimal=True))
                        angSCY = coordsSCY.separation(coordsPix)
                        radSCY = float(angSCY.to_string(unit=u.rad, decimal=True))
                        angZenith = coordsZenith.separation(coordsPix)
                        degZenith = float(angZenith.to_string(unit=u.deg, decimal=True))
                     #   aHtgLt[iR].Fill(cos(radSCZ), degZenith, (aSTART[iTI]+aSTOP[iTI])/2.-origin_time, tti*aAreaPix_array[iR][jpix])
                        aHtgLt[iR].Fill(cos(radSCZ), degZenith, tplot, tti*aAreaPix_array[iR][jpix])
                print iTI, aSTART[iTI], aRA_SCZ[iTI], aDEC_SCZ[iTI], tbdataSC.field('LAT_MODE')[iTI]#math.degrees(aAngSCY[1]), math.degrees(math.pi/2.-aAngSCY[0]), degZenith, math.degrees(radSCY)
            #if iTI%1==0:
             #   print iTI
                #rate = int((aSTOP[iTI]-metStart)/(metStop-metStart)*100.+0.5)
              #  rate = int(iTI/nTI*100.+0.5)
               # if rate>0:
                #    nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                 #   if iTI%2==0:
                  #      meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
                   # else:
                    #    meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                #else:
                 #   meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                #sys.stdout.write(meter
                #sys.stdout.write("\r{0}".format(tplot))
                #sys.stdout.write("coordsPix: {0}\n".format(coordsPix))
                #sys.stdout.write("angSCZ: {0}\n".format(angSCZ))
                #sys.stdout.write("angZenith: {0}\n\n".format(angZenith))
                sys.stdout.flush()
