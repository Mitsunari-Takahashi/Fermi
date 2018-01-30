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
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
import pandas
from pAnalysisConfig import *
from pFindHEALPix import *
from pLivetime import *


par = sys.argv
nameFileSuffix = par[1]

# Data
#if len(par)>6:
#    pathList = par[6]
#else:
pathList = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/PublicTableGRBs.xml" #/Users/Mitsunari/FermiAnalysis/catalogue/PublicTableGRBs.xml"
fileList = ET.parse(pathList)
rtXml = fileList.getroot()


# Spacecraft data
#pathFileScAll = "/disk/gamma/cta/store/takhsm/FermiData/spacecraft/FERMI_POINTING_FINAL_414_20?????_20?????_??.fits"
pathFileScAll = "/disk/gamma/cta/store/takhsm/FermiData/spacecraft/mtakahas-AstroServer-00011-ft2-30s.fits"

print "===================="
# Photon data
#listFileIn = par[5:]
#print listFileIn
listNameGrb = par[9:] #['130427324']
print listNameGrb
#catalogue
path_catlogue_csv = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/GBM-BusrtCatalogue_20170623.csv"
csv = pandas.read_csv(path_catlogue_csv)
num_lines = sum(1 for line in open(path_catlogue_csv))
rclass  =par[6]

# PSF histogram
path_file_perf = '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_{0}_perf.root'.format(rclass) #'/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root'
file_perf = ROOT.TFile(path_file_perf, 'READ')
print file_perf.GetName(), 'is opened.'
htg2_psf = file_perf.Get('psf_cth_q68_hist')
FIXED_PSF_ENERGY = float(par[4])
if FIXED_PSF_ENERGY<=0:
    print 'Scaling is variable for energy'
else:
    print 'Scaling is fixed for energy at', FIXED_PSF_ENERGY
FIXED_PSF_INCLIN = float(par[5])
if FIXED_PSF_INCLIN<=0:
    print 'Scaling is variable for inclination'
else:
    print 'Scaling is fixed for inclination at', FIXED_PSF_INCLIN

METSTARTLAT = 239557417. # 2008-08-04 15:43:36.000 UTC
METSTOPLAT = 536457605. # 2018Jan02 00:00

#for nameFileIn in listFileIn:
for nameGrb in listNameGrb:
    #------ Source data -----
    #indexGrbName = nameFileIn.rindex('GRB') + 3
    #indexGrbNameEnd = indexGrbName + 9
    #nameGrb = nameFileIn[indexGrbName:indexGrbNameEnd]

    for iGrb in range(num_lines-1):
        if int(nameGrb) == int(csv.ix[iGrb,'name']):
            raSrc = float(csv.ix[iGrb,'ra'])
            decSrc = float(csv.ix[iGrb,'dec']) 
            trigger_time = ConvertMjdToMet(float(csv.ix[iGrb,'trigger_time']))
            err_rad = float(csv.ix[iGrb,'error_radius'])
    # for grb in rtXml: 
    #     if grb[0].text==nameGrb: 
    #         raSrc = float(grb[5].text) 
    #         decSrc = float(grb[6].text) 
    #         trigger_time = float(grb[2].text) 
            if float(par[3])>0 and float(par[2])<METSTARTLAT and float(par[2])<float(par[3]):
                metStart = trigger_time+float(par[2])
                metStop = trigger_time+float(par[3])
                tPro = float(par[2])
                tPost = float(par[3])
            elif float(par[2])>=METSTARTLAT and float(par[3])>float(par[2]):
                metStart = float(par[2])
                metStop = float(par[3])
                tPro = metStart - trigger_time
                tPost = metStop - trigger_time
            else:
                print "Time domain is not correct."
            if float(par[7])==0 and float(par[8])==0:
                print 'No excluded time.'
                metExStart = 0
                metExStop = 0
            if float(par[7])<METSTARTLAT and float(par[8])<METSTARTLAT:
                metExStart = trigger_time-float(par[7])
                metExStop = trigger_time+float(par[8])
            elif float(par[7])>=METSTARTLAT and float(par[8])<float(par[7]):
                metExStart = float(par[7])
                metExStop = float(par[8])
            else:
                print "Excluded time domain is not correct."

            if metStart<METSTARTLAT:
                print "Your start time is outside of the histogram range!!!"
            if metStop>METSTOPLAT:
                print "Your stop time is outside of the histogram range!!!"

    print ""
    print "==============="
    print "GRB", nameGrb
    print "==============="
    print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 
    print "Time domain:", metStart, "-", metStop
    print "Excluded time:", metExStart, "-", metExStop

    # ON/OFF regions
    nOff = 0;
#    NHPSIDE_ON = 512
    NHPSIDE_ON = 32
    ANG_CUT = 5.
    ANG_CUT_RAD = radians(5.)

    aHpx_array = [find_pointon_healpxs(raSrc, decSrc, ANG_CUT, nhpside=NHPSIDE_ON)]
    aCoordsRegion = [SkyCoord(raSrc, decSrc, unit="deg")]
    aCoordsPix_array = []
    aAreaPix_array = []
    aAreaPix_sum = []
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
            area_pix = hppf.nside2pixarea(NHPSIDE_ON)
            aAreaPix_array[-1].append(area_pix)
        aAreaPix_sum.append(sum(aAreaPix_array[-1]))
        print aCoordsPix_array[-1]
#        print 'Solid angle:', aAreaPix_array[-1]
        print 'Solid angle =', aAreaPix_sum[-1], 'sr'
    # Output objects
    #fmw = ConvertMetToFMW(trigger_time)
    aFileToI = []
    if nameFileSuffix!="":
        nameFileSuffix = "_" + nameFileSuffix
    fileRoot = ROOT.TFile("Livetime_GRB{0}_T{1}-{2}{3}.root".format(nameGrb, par[2], par[3], nameFileSuffix), "update")
    aHtgLt = []
    NBIN_CTH = 50
    EDGE_CTH_LOW = 0.0
    EDGE_CTH_UP = 1.0
    NBIN_ZEN = 180
    EDGE_ZEN_LOW = 0
    EDGE_ZEN_UP = 180
    NBIN_ENE = htg2_psf.GetNbinsX()
    EDGE_ENE_LOW =  htg2_psf.GetXaxis().GetBinLowEdge(1)
    EDGE_ENE_UP =  htg2_psf.GetXaxis().GetBinUpEdge(htg2_psf.GetNbinsX())
    for hRegion in range(nOff+1):
        aHtgLt.append(ROOT.TH3D("htgLt_{0}".format(hRegion), "Livetime [sec sr];Cos(Inclination angle);Zenith angle [deg];Time from the GRB trigger [sec]".format(aStrRegion[hRegion], nameGrb), NBIN_CTH, EDGE_CTH_LOW, EDGE_CTH_UP, NBIN_ZEN, EDGE_ZEN_LOW, EDGE_ZEN_UP, max(10, int(METSTOPLAT-METSTARTLAT)/54000), METSTARTLAT, METSTOPLAT))#tPro, tPost))
    make_livetime_histogram(aHtgLt, nOff+1,pathFileScAll, metStart, metStop, aFileToI, aCoordsPix_array, aAreaPix_array, trigger_time, metExStart, metExStop)
    aHtgLt_projYX = []
    aHtgLt_scaled = []
    print 'Making output products...'
    fileRoot.cd()
    for jR in range(nOff+1):
        print 'Region', jR
        aHtgLt[jR].Write()
        aHtgLt_projYX.append(aHtgLt[jR].Project3D("yx"))
        aHtgLt_projYX[jR].Write()
        aHtgLt_scaled.append(ROOT.TH3D('{0}_scaled'.format(aHtgLt[jR].GetName()), '{0} scaled;{1};{2};log_{{10}}Energy [MeV]'.format(aHtgLt[jR].GetTitle(), aHtgLt_projYX[jR].GetXaxis().GetTitle(), aHtgLt_projYX[jR].GetYaxis().GetTitle()), NBIN_CTH, EDGE_CTH_LOW, EDGE_CTH_UP, NBIN_ZEN, EDGE_ZEN_LOW, EDGE_ZEN_UP, NBIN_ENE, EDGE_ENE_LOW, EDGE_ENE_UP))
        for iz in range(1, 1+NBIN_ENE):
            print '  Energy {0} - {1}'.format(aHtgLt_scaled[-1].GetZaxis().GetBinLowEdge(iz), aHtgLt_scaled[-1].GetZaxis().GetBinUpEdge(iz))
            if FIXED_PSF_ENERGY<=0:
                nbin_ene_psf = htg2_psf.GetXaxis().FindBin(aHtgLt_scaled[-1].GetZaxis().GetBinCenter(iz))
            else:
                nbin_ene_psf = htg2_psf.GetXaxis().FindBin(FIXED_PSF_ENERGY)
            for ix in range(1, 1+NBIN_CTH):
                if FIXED_PSF_INCLIN<=0:
                    nbin_inc_psf = htg2_psf.GetYaxis().FindBin(aHtgLt_scaled[-1].GetXaxis().GetBinCenter(ix))
                else:
                    nbin_inc_psf = htg2_psf.GetYaxis().FindBin(FIXED_PSF_INCLIN)
                psf_cut_rad = radians(htg2_psf.GetBinContent(nbin_ene_psf, nbin_inc_psf) + err_rad)
                area_ratio = 2.*pi*(1.0-cos(psf_cut_rad)) / aAreaPix_sum[jR]
                print '    Inclination {0} - {1} : Scaling factor = {2}'.format(aHtgLt_scaled[-1].GetXaxis().GetBinLowEdge(ix), aHtgLt_scaled[-1].GetXaxis().GetBinUpEdge(ix), area_ratio)
                for iy in range(1, 1+NBIN_ZEN):
                    aHtgLt_scaled[-1].SetBinContent(ix, iy, iz, aHtgLt_projYX[jR].GetBinContent(ix, iy)*area_ratio)
        aHtgLt_scaled[-1].Write()
    print 'Livetime calculation finished.'
