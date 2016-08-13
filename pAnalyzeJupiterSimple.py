#!/usr/bin/env python

import sys
from astropy.io import fits
from astropy.table import Table, Column
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from skyfield.api import load
from skyfield.api import utc
from array import array
import math
import numpy as np
import datetime
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
# Original modules
from pColor import *
from pMETandMJD import *
par = sys.argv
print par

planets = load('de421.bsp')
earth, jupiter, sun, moon = planets['EARTH'], planets['JUPITER BARYCENTER'], planets['SUN'], planets['MOON']
ts = load.timescale()

aOffset = [0., +6.]
aStrData = ['ON', 'OFF']
# SPACECRAFT DATA
aPathFitsSC = par[1:] #pathFitsPH.replace('_PH', '_SC')
print aPathFitsSC
for pathFitsSC in aPathFitsSC:
    print pathFitsSC
    tbSC = Table.read(pathFitsSC, hdu=1)
    nTI = len(tbSC)
    print nTI, "time intervals"

    for iData in range(len(aOffset)):
# PHOTON DATA
        pathFitsPH = pathFitsSC.replace('_SC00.fits', '{0}_PH00_gti.fits'.format(aStrData[iData]))
        tbPH = Table.read(pathFitsPH, hdu=1)
        tbPH.sort('TIME')
        nEvt = len(tbPH)

# Output
        atOffOffset = aOffset[iData] #[0.] #month The 1st element is the ON region
        aPathFileOut = pathFitsPH.replace('.fits', '.root') #.format(aStrRegion[-1]))
        aFileOut = ROOT.TFile(aPathFileOut, 'RECREATE')
# TTree of Events
        aTrEvt = ROOT.TTree('trEvt', 'EVENTS')
        aaAngDist = np.zeros(1, dtype=float)
        aaAngSepSun = np.zeros(1, dtype=float)
        aaAngSepMoon = np.zeros(1, dtype=float)
        aaEnergy = np.zeros(1, dtype=float)
        aaTime = np.zeros(1, dtype=float)
        aaRa = np.zeros(1, dtype=float)
        aaDec = np.zeros(1, dtype=float)
        aaL = np.zeros(1, dtype=float)
        aaB = np.zeros(1, dtype=float)
        aTrEvt.Branch('Energy',aaEnergy,'Energy/D')
        aTrEvt.Branch('Time',aaTime,'Time/D')
        aTrEvt.Branch('RA',aaRa,'RA/D')
        aTrEvt.Branch('DEC',aaDec,'DEC/D')
        aTrEvt.Branch('L',aaL,'L/D')
        aTrEvt.Branch('B',aaB,'B/D')
        aTrEvt.Branch('AngDist',aaAngDist,'AngDist/D')
        aTrEvt.Branch('AngularSeparationSun',aaAngSepSun,'AngularSeparatoinSun/D')
        aTrEvt.Branch('AngularSeparationMoon',aaAngSepMoon,'AngularSeparatoinMoon/D')

        nSecBin = 1
        htgGammaCount = ROOT.TH1D('htgGammaCount', 'Gamma-ray count curve', int((tbSC['STOP'][nTI-1]-tbSC['START'][0])/nSecBin), tbSC['START'][0], tbSC['STOP'][nTI-1])
        indexWeight = 2.0
        htgGammaCountWeight = ROOT.TH1D('htgGammaCountWeight', 'Gamma-ray count curve (Weight index = {0})'.format(indexWeight), int((tbSC['STOP'][nTI-1]-tbSC['START'][0])/nSecBin), tbSC['START'][0], tbSC['STOP'][nTI-1])

        dtOff = TimeDelta(atOffOffset*365.25/12., format='jd')
        print "Offset:", atOffOffset, "month"
        print dtOff
        timeStart = datetime.datetime.now()
        print timeStart
        
        iTI = 0
        for iEvt in range(nEvt):
            while tbSC['STOP'][iTI]<=tbPH['TIME'][iEvt]:
                iTI = iTI + 1
            if tbPH['TIME'][iEvt]>=tbSC['START'][iTI] and tbPH['TIME'][iEvt]<tbSC['STOP'][iTI]:
                coordsEvt = SkyCoord(tbPH['RA'][iEvt], tbPH['DEC'][iEvt], unit="deg")
                mjd = ConvertMetToMjd(tbPH['TIME'][iEvt])
                aaTime[0] = tbPH['TIME'][iEvt]
                aaEnergy[0] = tbPH['ENERGY'][iEvt]
                aaL[0] = tbPH['L'][iEvt]  
                aaB[0] =  tbPH['B'][iEvt]
                aaRa[0] = tbPH['RA'][iEvt]  
                aaDec[0] = tbPH['DEC'][iEvt]
                tOn = Time(mjd, format='mjd')
                tOff = Time(mjd, format='mjd') + dtOff #+ aDtRegion[iRegion]
                utcOn = ts.utc(tOn.to_datetime(timezone=utc))
                utcOff = ts.utc(tOff.to_datetime(timezone=utc))
                strTopoN = '{0} N'.format(tbSC['LAT_GEO'][iTI])
                strTopoE = '{0} E'.format(tbSC['LON_GEO'][iTI])
                sc = earth.topos(strTopoN, strTopoE)
                amJupiter = sc.at(utcOff).observe(jupiter)
                bJ, lJ, distGJ = amJupiter.galactic_latlon()
                degLJ, degBJ = lJ._degrees, bJ._degrees
                raJ, decJ, distJ = amJupiter.radec()
                degRaJ, degDecJ = raJ._degrees, decJ._degrees
                radRaJ, radDecJ = math.radians(degRaJ), math.radians(degDecJ)
                coordsJ = SkyCoord(degRaJ, degDecJ, unit="deg")
                angEvt = coordsEvt.separation(coordsJ)
                degEvt = float(angEvt.to_string(unit=u.deg, decimal=True))
                aaAngDist[0] = degEvt
                amSun = sc.at(utcOn).observe(sun)#amSun = earth.at(utcTI).observe(sun)
                raS, decS, distS = amSun.radec()
                degRaS, degDecS = raS._degrees, decS._degrees
                radRaS, radDecS = math.radians(degRaS), math.radians(degDecS)
                coordsS = SkyCoord(degRaS, degDecS, unit="deg")
                amMoon = sc.at(utcOn).observe(moon) #amMoon = earth.at(utcTI).observe(moon)
                raM, decM, distM = amMoon.radec()
                degRaM, degDecM = raM._degrees, decM._degrees
                radRaM, radDecM = math.radians(degRaM), math.radians(degDecM)
                coordsM = SkyCoord(degRaM, degDecM, unit="deg")
            #if degEvt < 2.0:
                angS = coordsS.separation(coordsEvt)
                degS = float(angS.to_string(unit=u.deg, decimal=True))
                aaAngSepSun[0] = degS
                angM = coordsM.separation(coordsEvt)
                degM = float(angM.to_string(unit=u.deg, decimal=True))
                aaAngSepMoon[0] = degM
                if degEvt < 5.0 and abs(degBJ)>=10 and aaAngSepSun[0]>=10 and aaAngSepMoon[0]>=10:
                    htgGammaCount.Fill(aaTime[0])
                    htgGammaCountWeight.Fill(aaTime[0], pow(aaEnergy[0]/1000., indexWeight))
                aTrEvt.Fill()
    
            if iEvt%(nEvt/200)==0:
                rate = int((iEvt*100.)/nEvt+0.5)
                if rate>0:
                    nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                    meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
            else:
                meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                sys.stdout.write(meter)
            sys.stdout.flush()
    
        aFileOut.cd()
        aTrEvt.Write()
        htgGammaCountWeight.Write()
        htgGammaCount.Write()
        htgMultipleEvt = ROOT.TH1D("htgMultipleEvt", "Multiple events in {0} sec".format(nSecBin), 10, 0, 10)
        for iBin in range(htgGammaCount.GetNbinsX()):
            nMulti = htgGammaCount.GetBinContent(iBin+1)
            htgMultipleEvt.Fill(nMulti)
            if nMulti>1:
                print nMulti, "event at MET", htgGammaCount.GetBinCenter(iBin)
        htgMultipleEvt.Write()
print ""        
print "Finished!"
