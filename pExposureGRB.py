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
nOff = 0;
aStrRegion = ["ON"]
for gRegion in range(nOff):
    aStrRegion.append("OFF{0}".format(gRegion))
# Spacecraft data
#pathFileScAll = "/Volumes/KARYU/FermiData/scdata/FERMI_POINTING_*_???_???????_???????_??.fits"
pathFileScAll = "/Volumes/KARYU/FermiData/scdata/lat_spacecraft_weekly_w???_p202_v???.fits"
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

            if float(par[4])>0 and float(par[3])<=0:
                metStart = trigger_time+float(par[3])
                metStop = trigger_time+float(par[4])
                tPro = float(par[3])
                tPost = float(par[4])
            elif float(par[3])>0 and float(par[4])>float(par[3]):
                metStart = float(par[3])
                metStop = float(par[4])
                tPro = metStart - trigger_time
                tPost = metStop - trigger_time
            else:
                print "Time domain is not correct."
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
    

    # Output objects
    fmw = ConvertMetToFMW(trigger_time)
    fmwStart = ConvertMetToFMW(metStart)
    fmwStop = ConvertMetToFMW(metStop)
    print "Fermi Mission Week:", fmwStart, "-", fmwStop
    coordsGrb = SkyCoord(raSrc, decSrc, unit="deg")
    aFileToI = []
    fileRoot = ROOT.TFile("Exposure_GRB{0}{1}.root".format(nameGrb, nameFileSuffix), "update")
    aHtgLt = []
    aHtgExp = []
    aHtgExp_ZE = []
    aHtgExp_ZT = []
    aHtgExpZ90 = []
    for hRegion in range(nOff+1):
        #aHtgLt.append(ROOT.TH1D("htgLt{0}".format(iRegion), "{0} Livetime".format(aStrRegion[iRegion]), 180, 0, 180))
        #aHtgExp.append(ROOT.TH2D("htgExp{0}".format(iRegion), "{0} Exposure".format(aStrRegion[iRegion]), 7, 4.35, 5.75, 180, 0, 180))
        aHtgLt.append(ROOT.TH2D("htgLt_GRB{0}_{1}".format(nameGrb, hRegion), "Livetime vs. zenith for GRB{1} (Spatially {0});Cos(Inclination angle);Zenith angle [deg];Live time [s]".format(aStrRegion[hRegion], nameGrb), 40, 0.2, 1.0, 180, 0, 180))
        aHtgExp.append(ROOT.TH3D('htgExp_GRB{0}_{1}'.format(nameGrb, hRegion), 'Exposure on zenith vs. inclination vs. energy for GRB{1} (Spatially {0});log10(Energy [MeV]);Cos(Inclination angle);Zenith angle [deg]'.format(aStrRegion[hRegion],nameGrb),7, 4.35, 5.75, 40, 0.2, 1.0, 180, 0, 180))
#        aHtgExp_ZE.append(ROOT.TH2D('htgExp_ZE_GRB{0}_{1}'.format(nameGrb, hRegion), 'Exposure on zenith vs. energy for GRB{1} (Spatially {0});log10(Energy [MeV]);Zenith angle [deg]'.format(aStrRegion[hRegion],nameGrb),7, 4.35, 5.75, 180, 0, 180))
        aHtgExp_ZT.append(ROOT.TH2D('htgExp_ZT_GRB{0}_{1}'.format(nameGrb, hRegion), 'Exposure on zenith vs. time for GRB{1} (Spatially {0});Time from the GRB trigger [sec];Zenith angle [deg]'.format(aStrRegion[hRegion],nameGrb), max(10, int(tPost-tPro)/54000), tPro, tPost, 180, 0, 180))
 
    for iFileSc in range(len(aPathFileScAll)):
        strPathFileSc = aPathFileScAll[iFileSc]
        #indexFileWeek = aPathFileScAll[iFileSc].rindex('FERMI_POINTING_FINAL') + 21
        #indexFileWeekEnd = indexFileWeek + 3
        #fmwFile = int(strPathFileSc[indexFileWeek:indexFileWeekEnd])
        #fmwFile = int(strPathFileSc[-27:-24])
        fmwFile = int(strPathFileSc[-18:-15])

        # Make file list
        if fmwFile>=int(fmwStart-1) and fmwFile<=int(fmwStop+1) :
            aFileToI.append(aPathFileScAll[iFileSc])
#         if fmwFile==int(fmw)-1:
#             aFileToI[0] = aPathFileScAll[iFileSc]
#         elif fmwFile==int(fmw):
#             aFileToI[1] = aPathFileScAll[iFileSc]
#         elif fmwFile==int(fmw)+1:
#             aFileToI[2] = aPathFileScAll[iFileSc]
    timeStart = datetime.datetime.now()
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
                    aHtgLt[iR].Fill(cos(radSCZ), degZenith, tti)
                    for iE in range(aHtgExp[iR].GetNbinsX()):
                        acc = htgAcc.GetBinContent(htgAcc.GetXaxis().FindBin(aHtgExp[iR].GetXaxis().GetBinCenter(iE+1)), htgAcc.GetYaxis().FindBin(cos(radSCZ)))
                        exp = acc * tti / 2. / math.pi / htgAcc.GetYaxis().GetBinWidth(1)
                        aHtgExp[iR].Fill(aHtgExp[iR].GetXaxis().GetBinCenter(iE+1), cos(radSCZ), degZenith, exp)
                        aHtgExp_ZT[iR].Fill((aSTART[iTI]+aSTOP[iTI])/2.-trigger_time, degZenith, exp)
            if iTI%(nTI/200)==0:
                rate = int((aSTOP[iTI]-metStart)/(metStop-metStart)*100.+0.5)
                if rate>0:
                    nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                    meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
                else:
                    meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                sys.stdout.write(meter)
                sys.stdout.flush()

    fileRoot.cd()
    for jR in range(nOff+1):
        aHtgLt[jR].Write()
        aHtgExp[jR].Write()
        aHtgExp_ZT[jR].Write()
        # aHtgExp_ZE.append(aHtgExp[jR].ProjectionX('htgExp_ZE_GRB{0}_{1}'.format(nameGrb, jR), 1, aHtgExp_ZE[jR].GetYaxis().FindBin(90.-0.1)))
        # aHtgExp_ZE[jR].SetTitle('Exposure with inclination vs. energy for GRB{1} (Spatially {0})'.format(aStrRegion[jR],nameGrb))
        # aHtgExp_ZE[jR].SetXTitle('log10(Energy [MeV])')
        # aHtgExp_ZE[jR].SetYTitle('Exposure [a.u.]')
        # aHtgExp_ZE[jR].Write()
        # aHtgExpZ90.append(aHtgExp[jR].ProjectionX('htgExpZ90_GRB{0}_{1}'.format(nameGrb, jR), 1, aHtgExp_ZE[jR].GetYaxis().FindBin(90.-0.1)))
        # aHtgExpZ90[jR].SetTitle('Exposure with zenith<90 vs. energy for GRB{1} (Spatially {0})'.format(aStrRegion[jR],nameGrb))
        # aHtgExpZ90[jR].SetXTitle('log10(Energy [MeV])')
        # aHtgExpZ90[jR].SetYTitle('Exposure [a.u.]')
        # aHtgExpZ90[jR].Write()
