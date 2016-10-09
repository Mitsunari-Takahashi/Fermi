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
from pCalcAngDist_Catalogue import *
from pLsList import ls_list

@click.command()
@click.argument('datafile', type=str)
@click.argument('suffix', type=str)
def main(datafile, suffix):

    #listFileIn = ls_list(datafile) 
    fileOut = ROOT.TFile("CompareCalOnlyFraction_{0}.root".format(suffix), "UPDATE")
    fileOut.cd()

    fileDat = ROOT.TFile(datafile)
    ch = fileDat.Get("MeritTuple")
    path_bdtfile = ROOT.TString(datafile)
    path_bdtfile.ReplaceAll(".root", "_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.root")
    fileBdt = ROOT.TFile(path_bdtfile.Data())
    ch.AddFriend("BDT = MeritTuple", fileBdt)

    ZENITH_CUT = 90.
    NBIN_ENE = 7
    LOWEDGE_ENE = 4.35
    UPEDGE_ENE = 5.75
    NBIN_CTH = 80
    LOWEDGE_CTH = 0.2
    UPEDGE_CTH = 1.0
    NBIN_CTE = 200
    LOWEDGE_CTE = 0
    UPEDGE_CTE = 200
    NBIN_BDT = 50
    LOWEDGE_BDT = 0
    UPEDGE_BDT = 5

    hRaw = ROOT.TH2F('hRaw', 'Raw number of events;log_{10}WP8CalOnlyEnergy;Cal1MomZDir', NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hCalOnly = ROOT.TH2F('hCalOnly', 'Number of CalOnly events;log_{10}WP8CalOnlyEnergy;Cal1MomZDir', NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hRaw_CalTwrEdge = ROOT.TH2F('hRaw_CalTwrEdge', 'Raw number of events;CalTwrEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hCalOnly_CalTwrEdge = ROOT.TH2F('hCalOnly_CalTwrEdge', 'Number of CalOnly events;CalTwrEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
#    hRaw_CalLATEdge = ROOT.TH2F('hRaw_CalLATEdge', 'Raw number of events;CalLATEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
#    hCalOnly_CalLATEdge = ROOT.TH2F('hCalOnly_CalLATEdge', 'Number of CalOnly events;CalLATEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hRaw_BDT = ROOT.TH2F('hRaw_BDT', 'Raw number of events;-log(1-(1+BDT)/2);Cal1MomZDir', NBIN_BDT, LOWEDGE_BDT, UPEDGE_BDT, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hCalOnly_BDT = ROOT.TH2F('hCalOnly_BDT', 'Number of CalOnly events;-log(1-(1+BDT)/2);Cal1MomZDir', NBIN_BDT, LOWEDGE_BDT, UPEDGE_BDT, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)

    NEVT = ch.GetEntries()
    print NEVT, 'events'
    for iEvt in range(NEVT):
        ch.GetEntry(iEvt)
        if ch.FT1CalZenithTheta<ZENITH_CUT and ch.Cal1RawEnergySum>=20000 and ch.FswGamState==0:
            hRaw.Fill(math.log10(max(1, ch.WP8CalOnlyEnergy)), ch.Cal1MomZDir)
            hRaw_CalTwrEdge.Fill(ch.CalTwrEdge, ch.Cal1MomZDir)
            hRaw_BDT.Fill(-math.log10(1.0-(1.0+ch.S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060)/2.0), ch.Cal1MomZDir)
            #hRaw_CalLATEdge.Fill(ch.CalLATEdge, ch.Cal1MomZDir)
            if  (ch.TkrNumTracks==0 or (math.log10(max(ch.CalTrackAngle,1E-4)) > (0.529795)*(ch.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((ch.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((ch.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((ch.EvtJointLogEnergy-3.000000)/0.916667,3))))*(ch.EvtJointLogEnergy >= 3.000000 and ch.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(ch.EvtJointLogEnergy >  5.750000)) ):
                hCalOnly.Fill(math.log10(max(1, ch.WP8CalOnlyEnergy)), ch.Cal1MomZDir)
                hCalOnly_CalTwrEdge.Fill(ch.CalTwrEdge, ch.Cal1MomZDir)
                #hCalOnly_CalLATEdge.Fill(ch.CalLATEdge, ch.Cal1MomZDir)
                hCalOnly_BDT.Fill(-math.log10(1.0-(1.0+ch.S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060)/2.0), ch.Cal1MomZDir)
    fileOut.cd()

    hRaw.Write()
    hCalOnly.Write()
    hFraction = hCalOnly.Clone("hFraction")
    hFraction.SetTitle("Fraction of CalOnly events")
    hFraction.Divide(hRaw)
    hFraction.Write()

    hRaw_CalTwrEdge.Write()
    hCalOnly_CalTwrEdge.Write()
    hFraction_CalTwrEdge = hCalOnly_CalTwrEdge.Clone("hFraction_CalTwrEdge")
    hFraction_CalTwrEdge.SetTitle("Fraction of CalOnly events")
    hFraction_CalTwrEdge.Divide(hRaw_CalTwrEdge)
    hFraction_CalTwrEdge.Write()

    hRaw_BDT.Write()
    hCalOnly_BDT.Write()
    hFraction_BDT = hCalOnly_BDT.Clone("hFraction_BDT")
    hFraction_BDT.SetTitle("Fraction of CalOnly events")
    hFraction_BDT.Divide(hRaw_BDT)
    hFraction_BDT.Write()


if __name__ == '__main__':
    main()
