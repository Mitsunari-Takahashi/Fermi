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

# PHOTON DATA
pathFitsPH = par[1]
hdulistPH = fits.open(pathFitsPH)
tbdataPH = hdulistPH[1].data
aENERGY = tbdataPH.field('ENERGY')
aRA = tbdataPH.field('RA')
aDEC = tbdataPH.field('DEC')
aL = tbdataPH.field('L')
aB = tbdataPH.field('B')
aTIME = tbdataPH.field('TIME')
nEvt = len(aENERGY)

# SPACECRAFT DATA
pathFitsSC = pathFitsPH.replace('_PH', '_SC')
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
nRegion = 2 # 1 ON region + other OFF resions
aStrRegion = []
aDtRegion = []
aLtRegion = []
aPathFileOut = []
aFileOut = []
# TTree of Light curve
aTrLC = []
aaTime = []
aaGamCountTi = []
aaGamCountTiWeight = []
aaLivetimeTi = []
aaCthTi = []
aaGalLTi = []
aaGalBTi = []
aaAngSepSunTi = []
aaAngSepMoonTi = []
# TTree of Events
aTrEvt = []
aaTimeInIt = []
aaEnergy = []
aaAngDist = []
aStrRegion.append('ON')
#aPathFileOut.append(pathFitsPH.replace('.fits', '_ON.root'))
tOffOffset = 3. #month
print "Offset for OFF regoin:", tOffOffset, "month"
timeStart = datetime.datetime.now()
print timeStart

for jRegion in range(nRegion):
    if jRegion>0:
        aStrRegion.append('OFF{0}'.format(jRegion))
    aDtRegion.append(TimeDelta(-jRegion*tOffOffset*365.25/12, format='jd'))
    aLtRegion.append(0.0)
    aPathFileOut.append(pathFitsPH.replace('.fits', '_{0}.root'.format(aStrRegion[-1])))
    aFileOut.append(ROOT.TFile(aPathFileOut[-1], 'RECREATE'))
    # LC tree
    aTrLC.append(ROOT.TTree('trLC_{0}'.format(aStrRegion[-1]), '{0} Light Curve'.format(aStrRegion[-1])))
    aaTime.append(np.zeros(1, dtype=float))
    aaGamCountTi.append(np.zeros(1, dtype=int))
    aaGamCountTiWeight.append(np.zeros(1, dtype=float))
    aaLivetimeTi.append(np.zeros(1, dtype=float))
    aaGalLTi.append(np.zeros(1, dtype=float))
    aaGalBTi.append(np.zeros(1, dtype=float))
    aaAngSepSunTi.append(np.zeros(1, dtype=float))
    aaAngSepMoonTi.append(np.zeros(1, dtype=float))
    aaCthTi.append(np.zeros(1, dtype=float))
    aTrLC[-1].Branch('Time',aaTime[-1],'Time/D')
    aTrLC[-1].Branch('GammaCounts',aaGamCountTi[-1],'GammaCounts/I')
    aTrLC[-1].Branch('WeightedGammaCounts',aaGamCountTiWeight[-1],'WeightedGammaCounts/D')
    aTrLC[-1].Branch('Livetime',aaLivetimeTi[-1],'Livetime/D')
    aTrLC[-1].Branch('GalL',aaGalLTi[-1],'GalL/D')
    aTrLC[-1].Branch('GalB',aaGalBTi[-1],'GalB/D')
    aTrLC[-1].Branch('AngularSeparationSun',aaAngSepSunTi[-1],'AngularSeparatoinSun/D')
    aTrLC[-1].Branch('AngularSeparationMoon',aaAngSepMoonTi[-1],'AngularSeparatoinMoon/D')
    aTrLC[-1].Branch('CosTheta',aaCthTi[-1],'CosTheta/D')
    # Evt tree
    aTrEvt.append(ROOT.TTree('trEvt_{0}'.format(aStrRegion[-1]), '{0} events'.format(aStrRegion[-1])))
    aaTimeInIt.append(np.zeros(1, dtype=float))
    aaEnergy.append(np.zeros(1, dtype=float))
    aaAngDist.append(np.zeros(1, dtype=float))
    aTrEvt[-1].Branch('TimeInInterval',aaTimeInIt[-1],'TimeInInterval/D')
    aTrEvt[-1].Branch('Energy',aaEnergy[-1],'Energy/D')
    aTrEvt[-1].Branch('AngularDistance',aaAngDist[-1],'AngularDistance/D)')

iEvt = 0
for iTI in range(nTI):
    #print "Time interval", iTI
    metTI = (aSTART[iTI] + aSTOP[iTI])/2.0
    mjdTI = ConvertMetToMjd(metTI)
    for iRegion in range(nRegion):
        #print "  ", aStrRegion[iRegion]
        aaTime[iRegion][0] = aSTART[iTI]
        #print aaTime[iRegion][0]
        aaGamCountTi[iRegion][0] = 0
        aaGamCountTiWeight[iRegion][0] = 0.0
        tTI = Time(mjdTI, format='mjd') + aDtRegion[iRegion]
#         yrTI = tTI.iso[0:4]
#         monTI = tTI.iso[5:7]
#         dTI = tTI.iso[8:10]
#         hTI = tTI.iso[11:13]
#         minTI = tTI.iso[14:16]
#         secTI = tTI.iso[17:]
        #utcTI = ts.utc(yrTI, monTI, dTI, hTI, minTI, secTI)
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
        coordsSCZ = SkyCoord(aRA_SCZ[iTI], aDEC_SCZ[iTI], unit="deg")
        angSCZ = coordsSCZ.separation(coordsJ)
        degSCZ = float(angSCZ.to_string(unit=u.deg, decimal=True))
        aaCthTi[iRegion][0] = math.cos(math.radians(degSCZ))
        aaGalLTi[iRegion][0] = degLJ
        aaGalBTi[iRegion][0] = degBJ
        aaAngSepSunTi[iRegion][0] = degS
        aaAngSepMoonTi[iRegion][0] = degM
        aaLivetimeTi[iRegion][0] = aLIVETIME[iTI]
        while iEvt<nEvt and aTIME[iEvt] < aSTART[iTI]:
            iEvt = iEvt + 1
        while iEvt<nEvt and aTIME[iEvt] < aSTOP[iTI]:
            coordsEvt = SkyCoord(aRA[iEvt], aDEC[iEvt], unit="deg")
            angEvt = coordsEvt.separation(coordsJ)
            degEvt = float(angEvt.to_string(unit=u.deg, decimal=True))
            #print aTIME[iEvt], aSTART[iTI]
            if degEvt < 1.0:
                aaGamCountTi[iRegion][0] = aaGamCountTi[iRegion][0]+1
                aaGamCountTiWeight[iRegion][0] = aaGamCountTiWeight[iRegion][0] + pow(aENERGY[iEvt]/1000., 2)
                aaTimeInIt[iRegion][0] = aTIME[iEvt] - aSTART[iTI]
                aaEnergy[iRegion][0] = aENERGY[iEvt]
                aaAngDist[iRegion][0] = degEvt
                aTrEvt[iRegion].Fill()
            iEvt = iEvt + 1
        aTrLC[iRegion].Fill()
        if aaGamCountTi[iRegion][0]>0:
            print "###################"
            print "#", aaGamCountTi[iRegion][0], aStrRegion[iRegion], "event!!! #"
            print "# MET", aaTime[iRegion][0], "#"
            print "# (", degRaJ, degDecJ, ") #"
            print "###################"

    if iTI%(nTI/200)==0:
        rate = int((iTI*100.)/nTI+0.5)
        if rate>0:
            nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
            meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
        else:
            meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
        sys.stdout.write(meter)
        sys.stdout.flush()
    
for kRegion in range(nRegion):
    aFileOut[kRegion].cd()
    aTrLC[kRegion].Write()
    aTrEvt[kRegion].Write()
print ""        
print "Finished!"
