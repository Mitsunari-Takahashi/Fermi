#!/usr/bin/env python

import ROOT
import sys
import click
import subprocess
from pLsList import ls_list
ROOT.gROOT.SetBatch()


@click.command()
@click.argument('name', type=str)
@click.option('definition', type=str)
@click.option('--ag', '-g', type=str, help="Path of AG event file")
@click.option('--bkg', '-b', type=str, help="Path of BKG event file")
@click.option('--logy', '-l', is_flag=True, help="Plot y-axis in logarithmic scale")

def main(name, definition, ag, bkg):
    fileAG = ROOT.TFile(ag, "update")
    print fileAG.GetName(), "is opened."
    trAG = fileAG.Get("MeritTuple")
    print trAG.GetName(), "is found."
    fileBKG = ROOT.TFile(bkg, "update")
    print fileBKG.GetName(), "is opened."
    trBKG = fileBKG.Get("MeritTuple")
    print trAG.GetName(), "is found."

    strCutBase = "(Cal1RawEnergySum>=20000 && FswGamState==0)"
    strCutCalOnly = "(TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy  >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) )"
    strCutHad = "( (McSourceId >= 1000 && McSourceId < 2000) || (McSourceId>=4000 && McSourceId<7000) )"
    strCutLep = "(McSourceId >= 2000 && McSourceId < 4000)"
    minAG = trAG.GetMinimum(strBase+strCalOnly)
    maxAG = trAG.GetMaximum(strBase+strCalOnly)
    minBKG = trBKG.GetMinimum(strBase+strCalOnly)
    maxBKG = trBKG.GetMaximum(strBase+strCalOnly)
    minCommon = min(minAG, minBKG)
    maxCommon = max(maxAG, maxBKG)
    NBIN_ENE = 7
    LOWEDGE_ENE = 4.35
    UPEDGE_ENE = 5.75
    NBIN_CTH = 4
    LOWEDGE_CTH = 0.2
    UPEDGE_CTH = 1.0
    NBIN_VAR = 100

    TH3F *h3Gam = ROOT.TH3F("h3Gam_"+name, "Gammas", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, minCommon, maxCommon)
    h3Gam.SetLineColor(kBlue)
    TH3F *h3Had = ROOT.TH3F("h3Had_"+name, "Hadrons", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, minCommon, maxCommon)
    h3Had.SetLineColor(kRed)
    TH3F *h3Lep = ROOT.TH3F("h3Lep_"+name, "Leptons", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, minCommon, maxCommon)
    h3Lep.SetLineColor(kMagenta)

    trAG.Draw(name+":-McZDir:McLogEnergy>>"+h3Gam.GetName(), strCutBase+strCalOnly, "goff")
    trBKG.Draw(name+":-McZDir:McLogEnergy>>"+h3Had.GetName(), strCutBase+strCutCalOnly+strCutHad, "goff")
    trBKG.Draw(name+":-McZDir:McLogEnergy>>"+h3Lep.GetName(), strCutBase+strCutCalOnly+strCutLep, "goff")

    cPlot = ROOT.TCanvas(name, definition, NBIN_ENE*150, NBIN_CTH*150)
    cPlot.Divide(NBIN_ENE, NBIN_CTH)

    arHsPlot = []

    for iE in range(NBIN_ENE):
        aHsPlot.append([])
        for iC in range(NBIN_CTH):
            aHsPlot[-1].append(ROOT.THStack("hs{0}_{1}_{2}".format(name, iE, iC), "{0}<=McLogEnergy<{1} and {2}<=-McZDir<{3}".format(h3Gam.GetXaxis().GetBinLowEdge(iE+1), h3Gam.GetXaxis().GetBinUpEdge(iE+1), h3Gam.GetYaxis().GetBinLowEdge(iC+1), h3Gam.GetYaxis().GetBinUpEdge(iC+1))))
            aHsPlot[-1].Add(h3Gam.ProjectionZ("{0}_Z".format(h3Gam.GetName()), iE+1, iE+1, iC+1, iC+1))
            aHsPlot[-1].Add(h3Had.ProjectionZ("{0}_Z".format(h3Had.GetName()), iE+1, iE+1, iC+1, iC+1))
            aHsPlot[-1].Add(h3Lep.ProjectionZ("{0}_Z".format(h3Lep.GetName()), iE+1, iE+1, iC+1, iC+1))
            cPlot.cd(iE*iC+1)
            aHsPlot.Draw("nostack")
            if logy==True:
                cPlot.cd(iE*iC+1).SetLogy()

    fileOut = ROOT.TFile("PlotVariables.root", "UPDATE")
    fileOut.cd()
    cPlot.Write()

if __name__ == '__main__':
    main()
