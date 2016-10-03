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
@click.argument('distcut', type=float)
@click.argument('rasrc', type=float)
@click.argument('decsrc', type=float)
@click.argument('start', type=float)
@click.argument('stop', type=float)
@click.argument('datafiles', type=str)
@click.argument('suffix', type=str)
def main(distcut, rasrc, decsrc, start, stop, datafiles, suffix):

    listFileIn = ls_list(datafiles) 

    vecSrc = vecterSkyPosition([rasrc, decsrc])
    fileOut = ROOT.TFile("CountsRawEvents_{0}.root".format(suffix), "UPDATE")
    fileOut.cd()
    htgOut = ROOT.TH1D("htgOut", "Raw event coutns within {0} deg from ({1}, {2})".format(distcut, rasrc, decsrc), 7, 4.35, 5.75)
    print htgOut.GetTitle()

    for pathfilein in listFileIn:
        filein = ROOT.TFile(pathfilein, 'READ')
        print filein.GetName()
        trMer = filein.Get("MeritTuple")
        nevt = trMer.GetEntries()
        print trMer.GetName(), 'has', nevt, 'events.'
        for iEvt in range(nevt):
            trMer.GetEntry(iEvt)
            if trMer.Cal1RawEnergySum>=20000 and (trMer.TkrNumTracks==0 or (math.log10(max(trMer.CalTrackAngle,1E-4)) > (0.529795)*(trMer.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((trMer.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((trMer.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((trMer.EvtJointLogEnergy-3.000000)/0.916667,3))))*(trMer.EvtJointLogEnergy >= 3.000000 and trMer.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(trMer.EvtJointLogEnergy >  5.750000)) ) and trMer.EvtElapsedTime>=start and trMer.EvtElapsedTime<stop and trMer.FswGamState==0:
                vecEvt = vecterSkyPosition([trMer.FT1CalRa, trMer.FT1CalDec])
                radDist = get_ang_dist_vectors(vecSrc, vecEvt)
                degDist = math.degrees(radDist)
                if degDist < distcut:
                    htgOut.Fill(trMer.EvtJointLogEnergy)
        print htgOut.GetEntries(), 'events have been filled.'
    fileOut.cd()
    htgOut.Write()
    print "=========="
    for ixbin in range(1, htgOut.GetNbinsX()+1):
        print htgOut.GetXaxis().GetBinLowEdge(ixbin), '-', htgOut.GetXaxis().GetBinUpEdge(ixbin)
        print htgOut.GetBinContent(ixbin), 'events'


if __name__ == '__main__':
    main()
