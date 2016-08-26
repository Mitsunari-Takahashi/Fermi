#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
from ROOT import TH2D
import numpy as np
from pAnalysisConfig import *
from pColor import *
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)


def main():
    fileRej = ROOT.TFile(par[5], 'READ')
    aHtgRej = []
    for icc in range(len(aaStrSelect[1])):
        aHtgRej.append(fileRej.Get('h2Rejection{0}'.format(icc)))

    # ----- Event class setup -----
    par = sys.argv
    cfg = ClassConfig('Both', [10, 3, 1], 1)
    aCutEGB = cfg.aCutEGB
    aaStrSelect = cfg.aaStrSelect
    aStrSelect = cfg.aStrSelect

    vecTgt = np.array([cos(radians(decSrc))*cos(radians(raSrc)), cos(radians(decSrc))*sin(radians(raSrc)), sin(radians(decSrc))])

    # Background estimation
    print "Estimation of unifomal backgrounds"
    nEvtPrecutPSF68 = [0., 0., 0.]
    nEvtPrecutPSF95 = [0., 0., 0.]
    nEvtPostcutPSF68 = [0., 0., 0.]
    nEvtPostcutPSF95 = [0., 0., 0.]
#chDat = ROOT.TChain('MeritTuple')
    for pathFileDat in listFileDat:
    #chDat.Add(fileDat)
        fileDat = ROOT.TFile(pathFileDat, 'READ')
        print fileDat.GetName()
        chDat = fileDat.Get('MeritTuple')
        nEvtBll = chDat.GetEntries()
        print chDat.GetName(), "has", nEvtBll, "events."

        for jEvt in range(nEvtBll):
            chDat.GetEntry(jEvt)
            if chDat.EvtElapsedTime>=metStart and chDat.EvtElapsedTime<metStop and chDat.FT1CalZenithTheta<90 and chDat.Cal1MomZDir>=0.2 and chDat.Cal1RawEnergySum>=20000 and chDat.EvtJointLogEnergy>=4.35 and chDat.EvtJointLogEnergy<5.75 and (chDat.TkrNumTracks==0 or (math.log10(max(chDat.CalTrackAngle,1E-4)) > (0.529795)*(chDat.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((chDat.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((chDat.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((chDat.EvtJointLogEnergy-3.000000)/0.916667,3))))*(chDat.EvtJointLogEnergy  >= 3.000000 and chDat.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(chDat.EvtJointLogEnergy >  5.750000))) and chDat.FswGamState == 0:
                vecEvtB = np.array([cos(radians(chDat.FT1CalDec))*cos(radians(chDat.FT1CalRa)), cos(radians(chDat.FT1CalDec))*sin(radians(chDat.FT1CalRa)), sin(radians(chDat.FT1CalDec))])
                radThetaB = acos(np.dot(vecTgt, vecEvtB))
                degDistB = degrees(radThetaB)
                aDictDistCutB = []
                for kcc in range(len(aaStrSelect[1])):
                    aDictDistCutB.append({ 'PSF95': (htgPerf.getPSF95_cth(1, kcc-1, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir) + err_rad), 'PSF68': (htgPerf.getPSF68_cth(1, kcc-1, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir) + err_rad) })
                    if degDistB<aDictDistCutB[kcc]['PSF95']:
                        nEvtPrecutPSF95[kcc] = nEvtPrecutPSF95[kcc] + 1
                        nEvtPostcutPSF95[kcc] = nEvtPostcutPSF95[kcc]+aHtgRej[kcc].GetBinContent(aHtgRej[kcc].GetXaxis().FindBin(chDat.EvtJointLogEnergy), aHtgRej[kcc].GetYaxis().FindBin(chDat.Cal1MomZDir))
                        if degDistB<aDictDistCutB[kcc]['PSF68']:
                            nEvtPrecutPSF68[kcc] = nEvtPrecutPSF68[kcc] + 1
                            nEvtPostcutPSF68[kcc] = nEvtPostcutPSF68[kcc]+aHtgRej[kcc].GetBinContent(aHtgRej[kcc].GetXaxis().FindBin(chDat.EvtJointLogEnergy), aHtgRej[kcc].GetYaxis().FindBin(chDat.Cal1MomZDir))
        for kcc in range(len(aaStrSelect[1])):
            print "Number of events within PSF95 after/before cut", aaStrSelect[1][kcc], ":", nEvtPostcutPSF95[kcc],"/",nEvtPrecutPSF95[kcc]
            print "Number of events within PSF68 after/before cut", aaStrSelect[1][kcc], ":", nEvtPostcutPSF68[kcc],"/",nEvtPrecutPSF68[kcc]
