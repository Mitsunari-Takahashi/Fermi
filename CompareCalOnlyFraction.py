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
@click.argument('datafiles', type=str)
@click.argument('suffix', type=str)
def main(datafiles, suffix):

    #listFileIn = ls_list(datafiles) 
    fileOut = ROOT.TFile("CompareCalOnlyFraction_{0}.root".format(suffix), "UPDATE")
    fileOut.cd()

    ch = ROOT.TChain("MeritTuple")
    ch.Add(datafiles)

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
    hRaw = ROOT.TH2F('hRaw', 'Raw number of events;log_{{10}}WP8CalOnlyEnergy;Cal1MomZDir', NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hCalOnly = ROOT.TH2F('hCalOnly', 'Number of CalOnly events;log_{{10}}WP8CalOnlyEnergy;Cal1MomZDir', NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hRaw_CalTwrEdge = ROOT.TH2F('hRaw_CalTwrEdge', 'Raw number of events;CalTwrEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
    hCalOnly_CalTwrEdge = ROOT.TH2F('hCalOnly_CalTwrEdge', 'Number of CalOnly events;CalTwrEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
#    hRaw_CalLATEdge = ROOT.TH2F('hRaw_CalLATEdge', 'Raw number of events;CalLATEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)
#    hCalOnly_CalLATEdge = ROOT.TH2F('hCalOnly_CalLATEdge', 'Number of CalOnly events;CalLATEdge;Cal1MomZDir', NBIN_CTE, LOWEDGE_CTE, UPEDGE_CTE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH)

    NEVT = ch.GetEntries()
    print NEVT, 'events'
    for iEvt in range(NEVT):
        ch.GetEntry(iEvt)
        if ch.FT1CalZenithTheta<ZENITH_CUT and ch.Cal1RawEnergySum>=20000 and ch.FswGamState==0:
            hRaw.Fill(math.log10(ch.WP8CalOnlyEnergy), ch.Cal1MomZDir)
            hRaw_CalTwrEdge.Fill(ch.CalTwrEdge, ch.Cal1MomZDir)
            #hRaw_CalLATEdge.Fill(ch.CalLATEdge, ch.Cal1MomZDir)
            if  (ch.TkrNumTracks==0 or (math.log10(max(ch.CalTrackAngle,1E-4)) > (0.529795)*(ch.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((ch.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((ch.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((ch.EvtJointLogEnergy-3.000000)/0.916667,3))))*(ch.EvtJointLogEnergy >= 3.000000 and ch.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(ch.EvtJointLogEnergy >  5.750000)) ):
                hCalOnly.Fill(math.log10(max(1, ch.WP8CalOnlyEnergy)), ch.Cal1MomZDir)
                hCalOnly_CalTwrEdge.Fill(ch.CalTwrEdge, ch.Cal1MomZDir)
                #hCalOnly_CalLATEdge.Fill(ch.CalLATEdge, ch.Cal1MomZDir)
    fileOut.cd()
    hRaw.Write()
    hCalOnly.Write()
    hFraction = hCalOnly.Clone("hFraction")
    hFraction.SetTitle("Fraction of CalOnly events")
    hFraction.Write()
    hRaw_CalTwrEdge.Write()
    hCalOnly_CalTwrEdge.Write()
    hFraction_CalTwrEdge = hCalOnly_CalTwrEdge.Clone("hFraction_CalTwrEdge")
    hFraction_CalTwrEdge.SetTitle("Fraction of CalOnly events")
    hFraction_CalTwrEdge.Write()

if __name__ == '__main__':
    main()
