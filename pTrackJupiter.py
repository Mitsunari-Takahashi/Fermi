#!/usr/bin/env python

import sys
from astropy.io import fits
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

# SPACECRAFT DATA
pathFitsSC = par[1] #pathFitsPH.replace('_PH', '_SC')
hdulistSC = fits.open(pathFitsSC)
tbdataSC = hdulistSC[1].data
aSTART, aSTOP = tbdataSC.field('START'), tbdataSC.field('STOP')
aLAT_GEO, aLON_GEO, aRAD_GEO = tbdataSC.field('LAT_GEO'), tbdataSC.field('LON_GEO'), tbdataSC.field('RAD_GEO')
aRA_SCZ, aRA_SCX = tbdataSC.field('RA_SCZ'), tbdataSC.field('RA_SCX')
aDEC_SCZ, aDEC_SCX = tbdataSC.field('DEC_SCZ'), tbdataSC.field('DEC_SCX')
aLIVETIME = tbdataSC.field('LIVETIME')
nTI = len(aSTART)
print nTI, "time intervals"

# Output
aOffOffset = [0., -4., -3., -2., -1., 2., 3., 4.] #[month] The 1st is the ON region
nRegion = len(aOffOffset) # 1 ON region + other OFF resions
aStrRegion = []
aDtRegion = []
aLtRegion = []
pathFileOut = pathFitsSC.replace('.fits', '_TrackJupiter.root')
fileOut = ROOT.TFile(pathFileOut, 'UPDATE')
# TTree of Light curve
aTrLC = []
aaStartTi = []
aaStopTi = []
aaLivetimeTi = []
aaRaTi = []
aaDecTi = []
aaGalLTi = []
aaGalBTi = []
aaDistTi = []
aaAngSepSunTi = []
aaAngSepMoonTi = []
aStrRegion.append('ON')
#aPathFileOut.append(pathFitsPH.replace('.fits', '_ON.root'))
print "Offset for ON/OFF regions:", aOffOffset, "month"
timeStart = datetime.datetime.now()
print timeStart

for jRegion in range(nRegion):
    if jRegion>0:
        aStrRegion.append('OFF{0}'.format(jRegion))
    aDtRegion.append(TimeDelta(aOffOffset[jRegion]*365.25/12, format='jd'))
    aLtRegion.append(0.0)
    # LC tree
    aTrLC.append(ROOT.TTree('trTrackJupiter_{0}'.format(aStrRegion[-1]), 'Track of Jupiter {0} month'.format(aOffOffset[jRegion])))
    aaStartTi.append(np.zeros(1, dtype=float))
    aaStopTi.append(np.zeros(1, dtype=float))
    aaLivetimeTi.append(np.zeros(1, dtype=float))
    aaRaTi.append(np.zeros(1, dtype=float))
    aaDecTi.append(np.zeros(1, dtype=float))
    aaGalLTi.append(np.zeros(1, dtype=float))
    aaGalBTi.append(np.zeros(1, dtype=float))
    aaDistTi.append(np.zeros(1, dtype=float))
    aaAngSepSunTi.append(np.zeros(1, dtype=float))
    aaAngSepMoonTi.append(np.zeros(1, dtype=float))
    aTrLC[-1].Branch('Start',aaStartTi[-1],'Start/D')
    aTrLC[-1].Branch('Stop',aaStopTi[-1],'Stop/D')
    aTrLC[-1].Branch('Livetime',aaLivetimeTi[-1],'Livetime/D')
    aTrLC[-1].Branch('Ra',aaRaTi[-1],'Ra/D')
    aTrLC[-1].Branch('Dec',aaDecTi[-1],'Dec/D')
    aTrLC[-1].Branch('GalL',aaGalLTi[-1],'GalL/D')
    aTrLC[-1].Branch('GalB',aaGalBTi[-1],'GalB/D')
    aTrLC[-1].Branch('Distance',aaDistTi[-1],'Distance/D')
    aTrLC[-1].Branch('AngularSeparationSun',aaAngSepSunTi[-1],'AngularSeparatoinSun/D')
    aTrLC[-1].Branch('AngularSeparationMoon',aaAngSepMoonTi[-1],'AngularSeparatoinMoon/D')

for iTI in range(nTI):
    metTI = (aSTART[iTI] + aSTOP[iTI])/2.0
    mjdTI = ConvertMetToMjd(metTI)
    for iRegion in range(nRegion):
        #print "  ", aStrRegion[iRegion]
        aaStartTi[iRegion][0] = aSTART[iTI]
        aaStopTi[iRegion][0] = aSTOP[iTI]
        aaLivetimeTi[iRegion][0] = aLIVETIME[iTI]
        #print aaTime[iRegion][0]
        tTI = Time(mjdTI, format='mjd') + aDtRegion[iRegion]
        utcTI = ts.utc(tTI.to_datetime(timezone=utc))
        #print "    ", utcTI.utc_datetime()
        strTopoN = '{0} N'.format(aLAT_GEO[iTI])
        strTopoE = '{0} E'.format(aLON_GEO[iTI])
        sc = earth.topos(strTopoN, strTopoE)
        #amJupiter = earth.at(utcTI).observe(jupiter)
        amJupiter = sc.at(utcTI).observe(jupiter)
        bJ, lJ, distGJ = amJupiter.galactic_latlon()
        degLJ, degBJ = lJ._degrees, bJ._degrees
        raJ, decJ, distJ = amJupiter.radec()
        degRaJ, degDecJ = raJ._degrees, decJ._degrees
        #print "     Target:", degRaJ, degDecJ
        aaGalLTi[iRegion][0] = degLJ
        aaGalBTi[iRegion][0] = degBJ
        aaRaTi[iRegion][0] = degRaJ
        aaDecTi[iRegion][0] = degDecJ
        aaDistTi[iRegion][0] = distJ.au #.to(u.au)
        radRaJ, radDecJ = math.radians(degRaJ), math.radians(degDecJ)
        coordsJ = SkyCoord(degRaJ, degDecJ, unit="deg")
        if iRegion==0:
            #amSun = earth.at(utcTI).observe(sun)
            amSun = sc.at(utcTI).observe(sun)
            raS, decS, distS = amSun.radec()
            degRaS, degDecS = raS._degrees, decS._degrees
            radRaS, radDecS = math.radians(degRaS), math.radians(degDecS)
            coordsS = SkyCoord(degRaS, degDecS, unit="deg")
            #amMoon = earth.at(utcTI).observe(moon)
            amMoon = sc.at(utcTI).observe(moon)
            raM, decM, distM = amMoon.radec()
            degRaM, degDecM = raM._degrees, decM._degrees
            radRaM, radDecM = math.radians(degRaM), math.radians(degDecM)
            coordsM = SkyCoord(degRaM, degDecM, unit="deg")
        angS = coordsS.separation(coordsJ)
        degS = float(angS.to_string(unit=u.deg, decimal=True))
        angM = coordsM.separation(coordsJ)
        degM = float(angM.to_string(unit=u.deg, decimal=True))
        #coordsSCZ = SkyCoord(aRA_SCZ[iTI], aDEC_SCZ[iTI], unit="deg")
        #angSCZ = coordsSCZ.separation(coordsJ)
        #degSCZ = float(angSCZ.to_string(unit=u.deg, decimal=True))
        #aaCthTi[iRegion][0] = math.cos(math.radians(degSCZ))
        aaAngSepSunTi[iRegion][0] = degS
        aaAngSepMoonTi[iRegion][0] = degM
        aTrLC[iRegion].Fill()

    if iTI%(nTI/200)==0:
        rate = int((iTI*100.)/nTI+0.5)
        if rate>0:
            nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
            meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
        else:
            meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
        sys.stdout.write(meter)
        sys.stdout.flush()
    
fileOut.cd()
for kRegion in range(nRegion):
    aTrLC[kRegion].Write()
print ""        
print "Finished!"
