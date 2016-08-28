#!/usr/bin/env python

import os
import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
from ROOT import TH2D
import healpy as hp
from healpy import pixelfunc as hppf
import numpy as np
import datetime
import click
from ctypes import *
import yaml
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

#from pCutBDT import cutBDT
from pFindHEALPix import *
from pAnalysisConfig import *
from pFindGrbInfo import *
from pLsList import ls_list

@click.command()
@click.option('--suffix', '-s',default="", help="Suffix for name of output product. GRB name is added automatically..")
@click.argument('evtfiles', type=str)
#@click.argument('datafiles', type=str)
@click.argument('start', type=float)
@click.argument('stop', type=float)
@click.option('--nhpside', '-n', default=16, help="NSIDE number for healpix. Use the same number as one for livetime calculation. ")
@click.option('--threshold', '-t', default=0, help="Flux threshold of catalogue sources which should be removed from Galactic OFF regions.")
def main(evtfiles, start, stop, suffix, nhpside, threshold):
    # Region setup
    pathCatalogue = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/gll_psch_v09.fit"
    aHpxGalOFF = find_galoff_healpxs(nhpside, threshold, pathCatalogue)

    # ----- Event class setup -----
    cfg = ClassConfig('Both', [10, 3, 1], 1)
    aCutEGB = cfg.aCutEGB
    aaStrSelect = cfg.aaStrSelect
    aStrSelect = cfg.aStrSelect

    #IRF
    listPathFilePerf = [['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                        ['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root']]
    htgPerf = CutPerformanceHtg(listPathFilePerf)

    print "===================="
    listFileIn = ls_list(evtfiles)
    print listFileIn

    aliasSelections = yaml.load(open("{0}/config/pass8_event_selections.yaml".format(os.environ.get("EVENTSELECT")),'r'))

    chIn = ROOT.TChain('EVENTS')
    for fileIn in listFileIn:
        chIn.Add(fileIn)
    if suffix!="":
        suffix = "_" + suffix

    print 'Analysis time domain: MET', start, '-', stop
    nEventChain = chIn.GetEntries()
    nEvent = chIn.GetEntries('TIME>={0} && TIME<{1}'.format(start, stop))
    print "Total number of events in the time domain:", nEvent

    # TTree
    fileOut = ROOT.TFile("Plot_GalOff{0}.root".format(suffix), 'UPDATE')
    fileOut.cd()
    trHighB = ROOT.TTree("EVENTS_GalOff", "Friend TTree for Galactic OFF")
    cbFlag = c_bool()
    trHighB.Branch('FLAG', cbFlag, 'FLAG/O')

    mgr = ROOT.TMultiGraph("mgr", "Gamma-like events in Galactic OFF regions")
    greOn = []
    mgrZenith = ROOT.TMultiGraph("mgrZenith", "Zenith angle of events in Galactic OFF regions")
    greZenith = []

    aHtgEvt = []
    NBIN_CTH = 40
    EDGE_CTH_LOW = 0.2
    EDGE_CTH_UP = 1.0
    NBIN_ZEN = 180
    EDGE_ZEN_LOW = 0
    EDGE_ZEN_UP = 180
    NBIN_ENE = 7
    EDGE_ENE_LOW =  4.35
    EDGE_ENE_UP =  5.75

    for pC in range(len(aaStrSelect)):
        greOn.append([])
        greZenith.append([])
        aHtgEvt.append([])
        for qC in range(len(aaStrSelect[pC])):
            greOn[-1].append(ROOT.TGraphErrors())
            greOn[-1][-1].SetName("greOn_{0}_{1}".format(pC, qC))
            greOn[-1][-1].SetTitle("{0}".format(aaStrSelect[pC][qC]))
            greOn[-1][-1].SetMarkerStyle(20)
            greZenith[-1].append(ROOT.TGraphErrors())
            greZenith[-1][-1].SetName("greZenith_{0}_{1}".format(pC, qC))
            greZenith[-1][-1].SetTitle("{0}".format(aaStrSelect[pC][qC]))
            greZenith[-1][-1].SetMarkerStyle(20)
            if pC==0:
                greOn[-1][-1].SetMarkerColor(13-12*qC)
                greZenith[-1][-1].SetMarkerColor(13-12*qC)
            elif pC==1:
                greOn[-1][-1].SetMarkerColor(kRed+3*(qC-2))
                greZenith[-1][-1].SetMarkerColor(kRed+3*(qC-2))
            mgr.Add(greOn[-1][-1])
            mgrZenith.Add(greZenith[-1][-1])  
            aHtgEvt[-1].append(ROOT.TH3D("htgEvt_GalOff{0}".format(aaStrSelect[pC][qC]), "{0} events in Galactic OFF regions (MET {1} - {2});Cos(Inclination angle);Zenith angle [deg];log_{{10}}Energy [MeV]".format(aaStrSelect[pC][qC], start, stop), NBIN_CTH, EDGE_CTH_LOW, EDGE_CTH_UP, NBIN_ZEN, EDGE_ZEN_LOW, EDGE_ZEN_UP, NBIN_ENE, EDGE_ENE_LOW, EDGE_ENE_UP))
    timeStart = datetime.datetime.now()
    for iEvent in range(nEventChain):
        cbFlag.value = 0
        chIn.GetEntry(iEvent)
        #vecEvt = np.array([cos(radians(chIn.dec))*cos(radians(chIn.ra)), cos(radians(chIn.dec))*sin(radians(chIn.ra)), sin(radians(chIn.dec))])
        npix = hppf.ang2pix(nhpside, math.pi/2.-math.radians(chIn.dec), math.radians(chIn.ra))
        if (npix in aHpxGalOFF) and chIn.t>=start and chIn.t<stop:
            cbFlag.value = 1
            if chIn.c == 1:
                if chIn.s == 4:
                    if chIn.z<90:
                        greOn[0][0].SetPoint(greOn[0][0].GetN(), chIn.t, pow(10, chIn.e-3))
                    greZenith[0][0].SetPoint(greZenith[0][0].GetN(), chIn.t, chIn.z)
                elif chIn.s == 128:
                    if chIn.z<90:
                        greOn[0][1].SetPoint(greOn[0][1].GetN(), chIn.t, pow(10, chIn.e-3))
                    greZenith[0][1].SetPoint(greZenith[0][1].GetN(), chIn.t, chIn.z)
                if chIn.s >= 4:
                    aHtgEvt[0][0].Fill(chIn.cth, chIn.z, chIn.e)
                if chIn.s >= 128:
                    aHtgEvt[0][1].Fill(chIn.cth, chIn.z, chIn.e)
            elif chIn.c == 2:
                if chIn.s == 4096:
                    if chIn.z<90:
                        greOn[1][0].SetPoint(greOn[1][0].GetN(), chIn.t, pow(10, chIn.e-3))
                    greZenith[1][0].SetPoint(greZenith[1][0].GetN(), chIn.t, chIn.z)
                elif chIn.s == 8192:
                    if chIn.z<90:
                        greOn[1][1].SetPoint(greOn[1][1].GetN(), chIn.t, pow(10, chIn.e-3))
                    greZenith[1][1].SetPoint(greZenith[1][1].GetN(), chIn.t, chIn.z)
                elif chIn.s == 16384:
                    if chIn.z<90:
                        greOn[1][2].SetPoint(greOn[1][2].GetN(), chIn.t, pow(10, chIn.e-3))
                    greZenith[1][2].SetPoint(greZenith[1][2].GetN(), chIn.t, chIn.z)
                if chIn.s >= 4096:
                    aHtgEvt[1][0].Fill(chIn.cth, chIn.z, chIn.e)
                if chIn.s >= 8192:
                    aHtgEvt[1][1].Fill(chIn.cth, chIn.z, chIn.e)
                if chIn.s >= 16384:
                    aHtgEvt[1][2].Fill(chIn.cth, chIn.z, chIn.e)

        trHighB.Fill()
        if iEvent%(nEventChain/20)==0:
            rate = int((iEvent*100.)/nEventChain+0.5)
            if rate>0:
                nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
            else:
                meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
            sys.stdout.write(meter)
            sys.stdout.flush()
    trHighB.AddFriend(chIn)
    cEvent = ROOT.TCanvas("cEvent", "Gamma-like events in Galactic OFF regions")
    cEvent.cd()
    mgr.Draw("AP")
    mgr.GetXaxis().SetTitle("Time [s]")
    mgr.GetYaxis().SetTitle("Energy [GeV]")
    leg = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            leg.AddEntry(greOn[pD][qD], greOn[pD][qD].GetTitle(), "p")
    leg.Draw("same")
    cZenith = ROOT.TCanvas("cZenith", "Zenith angle of Galactic OFF events")
    cZenith.cd()
    mgrZenith.Draw("AP")
    mgrZenith.GetXaxis().SetTitle("Time [s]")
    mgrZenith.GetYaxis().SetTitle("Zenith angle [deg]")
    legZenith = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            legZenith.AddEntry(greZenith[pD][qD], greZenith[pD][qD].GetTitle(), "p")

    legZenith.Draw("same")
    print ""
    fileOut.cd()
    cEvent.Write()
    cZenith.Write()
    trHighB.Write()
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            aHtgEvt[pD][qD].Sumw2()
            aHtgEvt[pD][qD].Write()
    print "Finished!"
    

if __name__ == '__main__':
    main()
