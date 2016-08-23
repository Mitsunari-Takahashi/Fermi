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
if len(par)>4:
    pathList = par[4]
else:
    pathList = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/PublicTableGRBs.xml" #/Users/Mitsunari/FermiAnalysis/catalogue/PublicTableGRBs.xml"
fileList = ET.parse(pathList)
rtXml = fileList.getroot()


# Spacecraft data
pathFileScAll = "/disk/gamma/cta/store/takhsm/FermiData/spacecraft/FERMI_POINTING_FINAL_414_20?????_20?????_??.fits"

print "===================="
# Photon data
#listFileIn = par[5:]
#print listFileIn
listNameGrb = ['160509374']

#for nameFileIn in listFileIn:
for nameGrb in listNameGrb:
    #------ Source data -----
    #indexGrbName = nameFileIn.rindex('GRB') + 3
    #indexGrbNameEnd = indexGrbName + 9
    #nameGrb = nameFileIn[indexGrbName:indexGrbNameEnd]
    for grb in rtXml: 
        if grb[0].text==nameGrb: 
            raSrc = float(grb[5].text) 
            decSrc = float(grb[6].text) 
            trigger_time = float(grb[2].text) 
            if float(par[3])>0 and float(par[2])<=0:
                metStart = trigger_time+float(par[2])
                metStop = trigger_time+float(par[3])
                tPro = float(par[2])
                tPost = float(par[3])
            elif float(par[2])>0 and float(par[3])>float(par[2]):
                metStart = float(par[2])
                metStop = float(par[3])
                tPro = metStart - trigger_time
                tPost = metStop - trigger_time
            else:
                print "Time domain is not correct."
            if grb[7].text == "--":
                err_rad = 0.
            else:
                err_rad = float(grb[7].text) 
    print ""
    print "==============="
    print "GRB", nameGrb
    print "==============="
    print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 
    print "Time domain:", metStart, "-", metStop

    # ON/OFF regions
    nOff = 0;
#    NHPSIDE_ON = 512
    NHPSIDE_ON = 64

    aHpx_array = [find_pointon_healpxs(raSrc, decSrc, 1.0, nhpside=NHPSIDE_ON)]
    aCoordsRegion = [SkyCoord(raSrc, decSrc, unit="deg")]
    aCoordsPix_array = []
    aAreaPix_array = []
    aStrRegion = []
    for (iRegion, coordsRegion) in enumerate(aCoordsRegion):
        if iRegion==0:
            aStrRegion.append("ON")
        else:
            aStrRegion.append("OFF{0}".format(iRegion))
        aCoordsPix_array.append([])
        aAreaPix_array.append([])
        for npix in aHpx_array[iRegion]:
            aAngPix = hppf.pix2ang(NHPSIDE_ON, npix)
            aCoordsPix_array[-1].append(SkyCoord(aAngPix[1], pi/2.-aAngPix[0], unit="rad"))
            aAreaPix_array[-1].append(hppf.nside2pixarea(NHPSIDE_ON, npix))
        print aCoordsPix_array[-1]
    # Output objects
    #fmw = ConvertMetToFMW(trigger_time)
    aFileToI = []
    fileRoot = ROOT.TFile("Livetime_GRB{0}{1}.root".format(nameGrb, nameFileSuffix), "update")
    aHtgLt = []
    for hRegion in range(nOff+1):
        aHtgLt.append(ROOT.TH3D("htgLt_GRB{0}_{1}".format(nameGrb, hRegion), "Livetime [sec sr];Cos(Inclination angle);Zenith angle [deg];Time from the GRB trigger [sec]".format(aStrRegion[hRegion], nameGrb), 40, 0.2, 1.0, 180, 0, 180, max(10, int(tPost-tPro)/54000), metStart, metStop))#tPro, tPost))
    make_livetime_histogram(aHtgLt, nOff+1,pathFileScAll, metStart, metStop, aFileToI, aCoordsPix_array, aAreaPix_array)#, trigger_time)
    aHtgLt_projYX = []
    fileRoot.cd()
    for jR in range(nOff+1):
        aHtgLt[jR].Write()
        aHtgLt_projYX.append(aHtgLt[jR].Project3D("yx"))
        aHtgLt_projYX[jR].Write()
