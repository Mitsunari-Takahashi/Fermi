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
import commands
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pMETandMJD import *
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

#from pCutBDT import cutBDT
from pAnalysisConfig import *

par = sys.argv
nameFileSuffix = par[1]
tPro = float(par[3])
tPost = float(par[4])

# Performance file
pathFilePerf = par[2]
filePerf = ROOT.TFile(pathFilePerf, 'READ')
htgAcc = filePerf.Get("acc_cth_hist")
print htgAcc.GetName(), "is found."

# Data
pathList = "/Users/Mitsunari/FermiAnalysis/catalogue/PublicTableGRBs.xml"
fileList = ET.parse(pathList)
rtXml = fileList.getroot()

# OFF regions
nOff = 4;
degOffOffset = 14.0;
aStrRegion = ["ON"]
for gRegion in range(nOff):
    aStrRegion.append("OFF{0}".format(gRegion))
# Spacecraft data
pathFileScAll = "/Volumes/KARYU/FermiData/scdata/FERMI_POINTING_*_???_???????_???????_??.fits"
cmd = "ls {0}".format(pathFileScAll)
ret = commands.getoutput(cmd)
aPathFileScAll = ret.split("\n")

print "===================="
# Photon data
#pathFileIn = par[5]
#cmd = "ls {0}".format(pathFileIn)
#ret = commands.getoutput(cmd)
#listFileIn = ret.split("\n")
listFileIn = par[5:]
print listFileIn

aHtgExpSum = []
aHtgLtSum = []
for hRegion in range(nOff+1):
    aHtgExpSum.append(ROOT.TH2D('htgExpSum{0}'.format(hRegion), '{0} summed exposure histogram'.format(aStrRegion[hRegion]),7, 4.35, 5.75, 180, 0, 180))
    aHtgLtSum.append(ROOT.TH1D("htgLtSum{0}".format(hRegion), "{0} summed livetime".format(aStrRegion[hRegion]), 180, 0, 180))

for nameFileIn in listFileIn:
    #------ Source data -----
    indexGrbName = nameFileIn.rindex('GRB') + 3
    indexGrbNameEnd = indexGrbName + 9
    nameGrb = nameFileIn[indexGrbName:indexGrbNameEnd]
    for grb in rtXml: #for iGrb in range(trList.GetEntries()):
        if grb[0].text==nameGrb: #if trList.name==int(nameGrb):
            raSrc = float(grb[5].text) #trList.ra
            decSrc = float(grb[6].text) #trList.dec
            trigger_time = float(grb[2].text) #ConvertMjdToMet(trList.trigger_time)
            metStart = trigger_time + tPro
            metStop = trigger_time + tPost
            if grb[7].text == "--":
                err_rad = 0.
            else:
                err_rad = float(grb[7].text) #trList.error_radius
    print ""
    print "==============="
    print "GRB", nameGrb
    print "==============="
    print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 
    print "Time domain:", metStart, "-", metStop
    fmw = ConvertMetToFMW(trigger_time)
    print "Fermi Mission Week:", fmw
    coordsGrb = SkyCoord(raSrc, decSrc, unit="deg")
    aFileToI = ["", "", ""]
    fileRoot = ROOT.TFile("Exposure_GRB{0}{1}.root".format(nameGrb, nameFileSuffix), "update")
    aHtgLt = []
    aHtgExp = []
    for iRegion in range(nOff+1):
        aHtgLt.append(ROOT.TH1D("htgLt{0}".format(iRegion), "{0} Livetime".format(aStrRegion[iRegion]), 180, 0, 180))
        aHtgExp.append(ROOT.TH2D("htgExp{0}".format(iRegion), "{0} Exposure".format(aStrRegion[iRegion]), 7, 4.35, 5.75, 180, 0, 180))
    for iFileSc in range(len(aPathFileScAll)):
        strPathFileSc = aPathFileScAll[iFileSc]
        #indexFileWeek = aPathFileScAll[iFileSc].rindex('FERMI_POINTING_FINAL') + 21
        #indexFileWeekEnd = indexFileWeek + 3
        #fmwFile = int(strPathFileSc[indexFileWeek:indexFileWeekEnd])
        fmwFile = int(strPathFileSc[-27:-24])
        if fmwFile==int(fmw)-1:
            aFileToI[0] = aPathFileScAll[iFileSc]
        elif fmwFile==int(fmw):
            aFileToI[1] = aPathFileScAll[iFileSc]
        elif fmwFile==int(fmw)+1:
            aFileToI[2] = aPathFileScAll[iFileSc]
    for fileToI in aFileToI:
        hdulistSC = fits.open(fileToI)
        tbdataSC = hdulistSC[1].data
        aSTART, aSTOP = tbdataSC.field('START'), tbdataSC.field('STOP')
        aRA_ZENITH = tbdataSC.field('RA_ZENITH')
        aDEC_ZENITH = tbdataSC.field('DEC_ZENITH')
        aRA_SCZ, aRA_SCX = tbdataSC.field('RA_SCZ'), tbdataSC.field('RA_SCX')
        aDEC_SCZ, aDEC_SCX = tbdataSC.field('DEC_SCZ'), tbdataSC.field('DEC_SCX')
        aLIVETIME = tbdataSC.field('LIVETIME')
        nTI = len(aSTART)
        print "  ", fileToI, "(", nTI, "intervals )"
        for iTI in range(nTI):
            if aSTART[iTI]>=metStart and aSTOP[iTI]<metStop:
                tti = aLIVETIME[iTI]
#                 elif aSTART[iTI]<metStart and aSTOP[iTI]<metStop:
#                     tti = aSTOP[iTI] - metStart
#                     print "Time interval No.", iTI, "LIVETIME:", tti, "sec"
#                 elif aSTART[iTI]>=metStart and aSTOP[iTI]>=metStop:
#                     tti = metStop - aSTART[iTI]
#                     print "Time interval No.", iTI, "LIVETIME:", tti, "sec"
                coordsSCZ = SkyCoord(aRA_SCZ[iTI], aDEC_SCZ[iTI], unit="deg")            
                angSCZ = coordsSCZ.separation(coordsGrb)
                radSCZ = float(angSCZ.to_string(unit=u.rad, decimal=True))
                coordsZenith = SkyCoord(aRA_ZENITH[iTI], aDEC_ZENITH[iTI], unit="deg")
                angZenith = coordsZenith.separation(coordsGrb)
                degZenith = float(angZenith.to_string(unit=u.deg, decimal=True))
                for iR in range(nOff+1):
                    aHtgLt[iR].Fill(degZenith, tti)
                    for iE in range(aHtgExp[iR].GetNbinsX()):
                        acc = htgAcc.GetBinContent(htgAcc.GetXaxis().FindBin(aHtgExp[iR].GetXaxis().GetBinCenter(iE+1)), htgAcc.GetYaxis().FindBin(cos(radSCZ)))
                        exp = acc * tti / 2. / math.pi / htgAcc.GetYaxis().GetBinWidth(1)
                        aHtgExp[iR].Fill(aHtgExp[iR].GetXaxis().GetBinCenter(iE+1), degZenith, exp)
    fileRoot.cd()
    for jR in range(nOff+1):
        aHtgExp[jR].Write()
        aHtgLt[jR].Write()
        aHtgExpSum[jR].Add(aHtgExp[jR])
        aHtgLtSum[jR].Add(aHtgLt[jR])
fileRootSum = ROOT.TFile('SummedExposure.root', 'UPDATE')
fileRootSum.cd()
for kR in range(nOff+1):
    aHtgExpSum[kR].Write()
    aHtgLtSum[kR].Write()
