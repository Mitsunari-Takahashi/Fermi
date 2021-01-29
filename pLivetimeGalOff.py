#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
import numpy as np
import xml.etree.ElementTree as ET
import datetime
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
import healpy as hp
from healpy import pixelfunc as hppf
#import commands
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pAnalysisConfig import *
from pFindHEALPix import *
from pLivetime import *


par = sys.argv
nameFileOut = par[2]
#nameFileSuffix = par[2]

# Spacecraft data
#pathFileScAll = "/disk/gamma/cta/store/takhsm/FermiData/spacecraft/FERMI_POINTING_FINAL_???_20?????_20?????_??.fits"
pathFileScAll = par[1]
#"/nfs/farm/g/glast/u/mtakahas/data/spacecraft/lat_spacecraft_weekly_w*.fits" 
#"/disk/gamma/cta/store/takhsm/FermiData/spacecraft/mtakahas-AstroServer-00009-ft2-30s.fits"

print "===================="
metStart = float(par[3])
metStop = float(par[4])
if len(par)>5:
    pathCatalogue = par[5]
else:
    pathCatalogue = "/nfs/farm/g/glast/u/mtakahas/data/catalogue/gll_psch_v13.fit" #"/disk/gamma/cta/store/takhsm/FermiData/catalogue/gll_psch_v09.fit"
print "Time domain:", metStart, "-", metStop
# ON/OFF regions
NHPSIDE_OFF = 16
aHpxOFF = [find_galoff_healpxs(NHPSIDE_OFF, 0, pathCatalogue)]
print aHpxOFF
aCoordsPix_array = [[]]
aAreaPix_array = [[]]
aStrRegion = ["GalacticOFF"]
for npix in aHpxOFF[0]:
    aAngPix = hppf.pix2ang(NHPSIDE_OFF, npix)
    aCoordsPix_array[0].append(SkyCoord(aAngPix[1], pi/2-aAngPix[0], unit="rad"))
    aAreaPix_array[0].append(hppf.nside2pixarea(NHPSIDE_OFF))
    
# Output objects
aFileToI = []
fileRoot = ROOT.TFile(nameFileOut, "update")
#fileRoot = ROOT.TFile("Livetime_GalacticOff_{0}.root".format(nameFileSuffix), "update")
aHtgLt = []
for hRegion in range(len(aStrRegion)):
    aHtgLt.append(ROOT.TH3D("htgLt_{0}".format(aStrRegion[hRegion]), "Livetime over {0} [sec sr];Cos(Inclination angle);Zenith angle [deg];Time[sec]".format(aStrRegion[hRegion]), 100, -1.0, 1.0, 180, 0, 180, max(10, int(metStop-metStart)/54000), metStart, metStop+1))
make_livetime_histogram(aHtgLt, len(aStrRegion), pathFileScAll, metStart, metStop, aFileToI, aCoordsPix_array, aAreaPix_array)
aHtgLt_projYX = []
fileRoot.cd()
for jR in range(len(aStrRegion)):
    aHtgLt[jR].Write()
    aHtgLt_projYX.append(aHtgLt[jR].Project3D("yx"))
    aHtgLt_projYX[jR].Write()
