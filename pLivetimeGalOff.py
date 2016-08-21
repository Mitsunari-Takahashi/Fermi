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
#from pColor import *
from pAnalysisConfig import *
from pFindHEALPix import *
from pLivetime import *


par = sys.argv
nameFileSuffix = par[1]

# Data
#if len(par)>4:
#    pathList = par[4]
#else:
#    pathList = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/PublicTableGRBs.xml" #/Users/Mitsunari/FermiAnalysis/catalogue/PublicTableGRBs.xml"
#fileList = ET.parse(pathList)
#rtXml = fileList.getroot()


# Spacecraft data
pathFileScAll = "/disk/gamma/cta/store/takhsm/FermiData/spacecraft/FERMI_POINTING_FINAL_???_20?????_20?????_??.fits"

print "===================="
# Photon data
#listFileIn = par[5:]
#print listFileIn
#listNameGrb = ['160509374']

#for nameFileIn in listFileIn:
#for nameGrb in listNameGrb:
    #------ Source data -----
    #indexGrbName = nameFileIn.rindex('GRB') + 3
    #indexGrbNameEnd = indexGrbName + 9
    #nameGrb = nameFileIn[indexGrbName:indexGrbNameEnd]
    # for grb in rtXml: 
    #     if grb[0].text==nameGrb: 
    #         raSrc = float(grb[5].text) 
    #         decSrc = float(grb[6].text) 
    #         trigger_time = float(grb[2].text) 
    #         if float(par[3])>0 and float(par[2])<=0:
    #             metStart = trigger_time+float(par[2])
    #             metStop = trigger_time+float(par[3])
    #             tPro = float(par[2])
    #             tPost = float(par[3])
    #         elif float(par[2])>0 and float(par[3])>float(par[2]):
    #             metStart = float(par[2])
    #             metStop = float(par[3])
    #             tPro = metStart - trigger_time
    #             tPost = metStop - trigger_time
    #         else:
    #             print "Time domain is not correct."
    #         if grb[7].text == "--":
    #             err_rad = 0.
    #         else:
    #             err_rad = float(grb[7].text) 
    # print ""
    # print "==============="
    # print "GRB", nameGrb
    # print "==============="
    # print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 
metStart = float(par[2])
metStop = float(par[3])
if len(par)>4:
    pathCatalogue = par[4]
else:
    pathCatalogue = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/gll_psch_v09.fit"
print "Time domain:", metStart, "-", metStop
# ON/OFF regions
NHPSIDE_OFF = 16
aHpxOFF = [find_galoff_healpxs(NHPSIDE_OFF, 0, pathCatalogue)]
aCoordsPix_array = [[]]
aAreaPix_array = [[]]
aStrRegion = ["GalacticOFF"]
for npix in aHpxOFF[0]:
    aAngPix = hppf.pix2ang(NHPSIDE_OFF, npix)
    aCoordsPix_array[0].append(SkyCoord(pi/2.-aAngPix[1], aAngPix[0], unit="deg"))
    aAreaPix_array[0].append(hppf.nside2pixarea(NHPSIDE_OFF, npix))
    
# Output objects
aFileToI = []
fileRoot = ROOT.TFile("Livetime_GalacticOff_{0}.root".format(nameFileSuffix), "update")
aHtgLt = []
for hRegion in range(len(aStrRegion)):
    aHtgLt.append(ROOT.TH3D("htgLt_{0}".format(aStrRegion[hRegion]), "Livetime over {0} [sec sr];Cos(Inclination angle);Zenith angle [deg];Time[sec]".format(aStrRegion[hRegion]), 40, 0.2, 1.0, 180, 0, 180, max(10, int(metStop-metStart)/54000), metStart, metStop+1))
make_livetime_histogram(aHtgLt, len(aStrRegion), pathFileScAll, metStart, metStop, aFileToI, aCoordsPix_array, aAreaPix_array)
aHtgLt_projYX = []
fileRoot.cd()
for jR in range(len(aStrRegion)):
    aHtgLt[jR].Write()
    aHtgLt_projYX.append(aHtgLt[jR].Project3D("yx"))
    aHtgLt_projYX[jR].Write()
