#!/usr/bin/env python
"""Count number of events within a certain angular distance from a sky position.
"""
import sys
from array import array
import math
import numpy as np
import click
import ROOT
from ROOT import TTree
import xml.etree.ElementTree as ET
from pCalcAngDist_Catalogue import *
from pAnalysisConfig import *
from pFindGrbInfo import *
from pLsList import ls_list

@click.command()
@click.argument('grbname', type=str)
#@click.argument('distcut', type=float)
#@click.argument('rasrc', type=float)
#@click.argument('decsrc', type=float)
@click.argument('start', type=float)
@click.argument('stop', type=float)
@click.argument('datafiles', type=str)
@click.argument('suffix', type=str)
@click.option('--fixpsfenergy', '-e', type=float, default=0.0, help="Set energy in log scale if you will fix the PSF cut on energy.")
@click.option('--fixpsfinclin', '-i', type=float, default=0.0, help="Set cos(inclination angle) if you will fix the PSF cut on inclination.")
def main(grbname, start, stop, datafiles, suffix, fixpsfenergy, fixpsfinclin):

    # Catalogue Data
    pathList = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/PublicTableGRBs.xml" #"/nfs/farm/g/glast/u/mtakahas/data/catalogue/PublicTableGRBs.xml"
    fileList = ET.parse(pathList)
    rtXml = fileList.getroot()
    dict_grb = find_grb_info(grbname, rtXml)
    rasrc = dict_grb["RA"]
    decsrc = dict_grb["DEC"]
    ZENITH_CUT = 90.
    CALONLY_CLASS = 0 # 0: R100, 1:R30, 2:R10
    listFileIn = ls_list(datafiles) 
    
    #IRF
    listPathFilePerf = [['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                        ['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root']]
    htgPerf = CutPerformanceHtg(listPathFilePerf)

    vecSrc = vecterSkyPosition([rasrc, decsrc])
    fileOut = ROOT.TFile("CountsRawEvents_{0}.root".format(suffix), "UPDATE")
    fileOut.cd()
    htgOut = ROOT.TH2D("htgOut", "Raw event coutns within PSF68 from ({0}, {1})".format(rasrc, decsrc), 7, 4.35, 5.75, 10, 0.0, 1.0)
    print htgOut.GetTitle()

    for pathfilein in listFileIn:
        filein = ROOT.TFile(pathfilein, 'READ')
        print filein.GetName()
        trMer = filein.Get("MeritTuple")
        nevt = trMer.GetEntries()
        print trMer.GetName(), 'has', nevt, 'events.'
        for iEvt in range(nevt):
            trMer.GetEntry(iEvt)
            if trMer.Cal1RawEnergySum>=20000 and (trMer.TkrNumTracks==0 or (math.log10(max(trMer.CalTrackAngle,1E-4)) > (0.529795)*(trMer.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((trMer.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((trMer.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((trMer.EvtJointLogEnergy-3.000000)/0.916667,3))))*(trMer.EvtJointLogEnergy >= 3.000000 and trMer.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(trMer.EvtJointLogEnergy >  5.750000)) ) and trMer.FT1CalZenithTheta<ZENITH_CUT and trMer.EvtElapsedTime>=start and trMer.EvtElapsedTime<stop and trMer.FswGamState==0:
                vecEvt = vecterSkyPosition([trMer.FT1CalRa, trMer.FT1CalDec])
                radDist = get_ang_dist_vectors(vecSrc, vecEvt)
                degDist = math.degrees(radDist)

                if fixpsfenergy==0:
                    epsf = math.log10(trMer.WP8CalOnlyEnergy)
                else:
                    epsf = fixpsfenergy
                if fixpsfinclin==0:
                    cthpsf = trMer.Cal1MomZDir
                else:
                    cthpsf = fixpsfinclin
                #dictDistCut = { 'PSF95': (htgPerf.getPSF95_cth(1, CALONLY_CLASS, epsf, cthpsf) + dict_grb["ERROR_RADIUS"]),  'PSF68': (htgPerf.getPSF68_cth(1, CALONLYCLASS, epsf, cthpsf) + dict_grb["ERROR_RADIUS"])}
                distcut = htgPerf.getPSF68_cth(1, CALONLY_CLASS, epsf, cthpsf) + dict_grb["ERROR_RADIUS"]
                if degDist < distcut:
                    htgOut.Fill(math.log10(trMer.WP8CalOnlyEnergy), trMer.Cal1MomZDir)
        print htgOut.GetEntries(), 'events have been filled.'
    fileOut.cd()
    htgOut.Write()


if __name__ == '__main__':
    main()
