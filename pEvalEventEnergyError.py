#!/usr/bin/env python
"""For evaluation of energy uncertainty of one event.
"""
import sys
import math
from math import log10, acos, atan2, pi, degrees, sin, cos, asin
#import numpy as np
import click
from ctypes import *
from array import array
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
from array import array
from pColor import *
ROOT.gROOT.SetBatch() 

@click.command()
@click.argument('e', type=float)
@click.argument('bep', type=float)
@click.option('--fileCut', '-c', default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/hBDT_CutValues_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.root')
@click.option('--fileBDT', '-b', default='/u/gl/mtakahas/work/data/MC/AG200909_62_2016Jun_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.root')
@click.option('--fileOut', '-o', default='hEnergyUncertainty.root')
@click.option('--htgCut', '-t', default='htgCutBDT0')
@click.option('--bdt', '-v', default='S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060')
def main(e, bep, filecut, filebdt, fileout, htgcut, bdt):
    eLin = pow(10, e-3)
    fileAG = ROOT.TFile('/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2016Jun.root')
    trAG = fileAG.Get('MeritTuple')
    nEvt = trAG.GetEntries()
    fileBDT = ROOT.TFile(filebdt)
    trBDT = fileBDT.Get('MeritTuple')
    trAG.AddFriend(trBDT, 'trTempBDT')
    bdtR = c_float() 
    print bdt
    trAG.SetBranchAddress('{0}'.format(bdt), bdtR)
    aPowIndex = [0, 1, 2, 3]
    dictHtg2D = {}
    dictHtgProjX = {}
    fileOutEU = ROOT.TFile(fileout, 'RECREATE')
    for pi in aPowIndex:
        dictHtg2D[pi] = ROOT.TH2D('htg2D_{0}'.format(int(pi*10)), 'Power-law index -{0};EvtJointLogEnergy-{1};McLogEnergy-{1}'.format(pi, e), 100, -1.0, 1.0, 100, -1.0, 1.0)
    fileCutBDT = ROOT.TFile(filecut)
    htgCut = fileCutBDT.Get(htgcut)
    for iEvt in range(nEvt):
        trAG.GetEntry(iEvt)
        if trAG.Cal1MomZCrossSide840>=0.0 and trAG.Cal1MomNumIterations>0 and trAG.Cal1TransRms>=10 and trAG.Cal1TransRms<70 and trAG.Acd2Cal1VetoSigmaHit>0 and trAG.Cal1MomZDir>=0.2 and -log10(1.0-(1.0+bdtR.value)/2.)>=htgCut.GetBinContent(htgCut.GetXaxis().FindBin(trAG.EvtJointLogEnergy),htgCut.GetYaxis().FindBin(trAG.Cal1MomZDir)) and trAG.Cal1RawEnergySum>=20000 and (trAG.TkrNumTracks==0 or (math.log10(max(trAG.CalTrackAngle,1E-4)) > (0.529795)*(trAG.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(math.pow((trAG.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(math.pow((trAG.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(math.pow((trAG.EvtJointLogEnergy-3.000000)/0.916667,3))))*(trAG.EvtJointLogEnergy >= 3.000000 and trAG.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(trAG.EvtJointLogEnergy >  5.750000))) and abs((trAG.WP8CalOnlyBEPCaseE_myBDT+1.)/2.-bep)<0.01 and trAG.FswGamState==0: #and (trAG.GltGemSummary&0x20)==0):
            for pi in aPowIndex:
                dictHtg2D[pi].Fill(trAG.McLogEnergy-e, trAG.EvtJointLogEnergy-e, pow(10, (-pi+1)*trAG.McEnergy/pow(10, e)))
    aPos = [0.68, 0.95]
    aPosLeft = array('d', [1.-aPos[0], 1-aPos[1]])
    aPosRight = array('d', [aPos[0], aPos[1]])
    aQuantLeft = array('d', [0., 0.])
    aQuantRight = array('d', [0., 0.])
    fileOutEU.cd()
    for pi in aPowIndex:
        print '=====  Index:', pi, '  ====='
        dictHtg2D[pi].Write()
        nBinCutLowY = dictHtg2D[pi].GetYaxis().FindBin(-0.005)
        nBinCutUpY = dictHtg2D[pi].GetYaxis().FindBin(0.005)
        print 'Bin:', nBinCutLowY, '-', nBinCutUpY
        dictHtgProjX[pi] = dictHtg2D[pi].ProjectionX('htgProjX_{0}'.format(int(pi*10)), nBinCutLowY, nBinCutUpY)
        dictHtgProjX[pi].Write()
        hLeft = dictHtgProjX[pi].Clone('hTempLeft')
        hRight = dictHtgProjX[pi].Clone('hTempRight')
        for iBin in range(dictHtgProjX[pi].GetNbinsX()):
            if dictHtgProjX[pi].GetXaxis().GetBinCenter(iBin+1)>0:
                hLeft.SetBinContent(iBin+1, 0)
            elif dictHtgProjX[pi].GetXaxis().GetBinCenter(iBin+1)<0:
                hRight.SetBinContent(iBin+1, 0)
        hLeft.GetQuantiles(len(aPos), aQuantLeft, aPosLeft)
        print aQuantLeft              
        hRight.GetQuantiles(len(aPos), aQuantRight, aPosRight)
        print aQuantRight              
        for (iPos, nPos) in enumerate(aPos):
            print '{0}% Error: {1}-{2}+{3} GeV'.format(int(nPos*100), eLin, eLin-pow(10, e+aQuantLeft[iPos]-3), pow(10, e+aQuantRight[iPos]-3)-eLin)
        

if __name__ == '__main__':
    main()
