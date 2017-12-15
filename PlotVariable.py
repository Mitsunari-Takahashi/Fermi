#!/usr/bin/env python

import ROOT
import sys
import math
import click
import subprocess
from pLsList import ls_list
ROOT.gROOT.SetBatch()


@click.command()
@click.argument('name', type=str)
@click.argument('definition', type=str)
@click.option('--ag', '-g', type=str, default='/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2016Jun.root', help="Path of AG event file")
@click.option('--bkg', '-b', type=str, default='/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_Combined.root', help="Path of BKG event file")
@click.option('--logy', '-l', is_flag=True, help="Plot y-axis in logarithmic scale")
@click.option('--nbine', '-e', type=int, default=7, help="Number of energy bin")
@click.option('--nbincth', '-c', type=int, default=4, help="Number of cos(theta) bin")

def main(name, definition, ag, bkg, nbine, nbincth, logy):
    NBIN_ENE = nbine
    LOWEDGE_ENE = 4.35
    UPEDGE_ENE = 5.75
    NBIN_CTH = nbincth
    LOWEDGE_CTH = 0.2
    UPEDGE_CTH = 1.0
    NBIN_VAR = 100
    ncx = 1
    ncy = 1
    if NBIN_ENE==1 or NBIN_CTH==1:
        ncy = int(math.sqrt(NBIN_ENE*NBIN_CTH))
        ncx = int(math.ceil(NBIN_CTH/ncy))
    else:
        ncy = NBIN_CTH
        ncx = NBIN_ENE
    cPlot = ROOT.TCanvas(name, definition, int(ncx*160), int(ncy*160))
    cPlot.Divide(ncx, ncy)
#    c = ROOT.TCanvas(name, definition, int(ncx*160), int(ncy*160))
#    cPlot.Divide(ncx, ncy)

    fileAG = ROOT.TFile(ag, "read")
    print fileAG.GetName(), "is opened."
    trAG = fileAG.Get("MeritTuple")
    print trAG.GetName(), "is found."
    fileBKG = ROOT.TFile(bkg, "read")
    print fileBKG.GetName(), "is opened."
    trBKG = fileBKG.Get("MeritTuple")
    print trAG.GetName(), "is found."

    strCutBase = "( Cal1RawEnergySum>=20000 && FswGamState==0 )" # && Cal1MomZCrossSide840>=0.0 && Cal1MomZDir>=0.2 && Cal1MomNumIterations>0
    strCutCalOnly = "( TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy  >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000))  )"
    strCutHad = "( (McSourceId >= 1000 && McSourceId < 2000) || (McSourceId>=4000 && McSourceId<7000) )"
    strCutLep = "( McSourceId >= 2000 && McSourceId < 4000 )"
    minAG = trAG.GetMinimum(strCutBase+strCutCalOnly)
    maxAG = trAG.GetMaximum(strCutBase+strCutCalOnly)
    minBKG = trBKG.GetMinimum(strCutBase+strCutCalOnly)
    maxBKG = trBKG.GetMaximum(strCutBase+strCutCalOnly)
    minCommon = min(minAG, minBKG)
    maxCommon = max(maxAG, maxBKG)

    fileOut = ROOT.TFile("PlotVariables.root", "UPDATE")
    fileOut.cd()

    h3Gam = ROOT.TH3F("h3Gam_"+name, "Gammas", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, minCommon, maxCommon)
    h3Gam.SetLineColor(ROOT.kBlue)
    h3Gam.SetLineWidth(2)
    h3Gam.SetFillStyle(3004)
    h3Gam.SetFillColor(ROOT.kBlue)
    h3Gam.Write()
    h3Had = ROOT.TH3F("h3Had_"+name, "Hadrons", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, minCommon, maxCommon)
    h3Had.SetLineColor(ROOT.kRed)
    h3Had.SetLineWidth(2)
    h3Had.SetFillStyle(3005)
    h3Had.SetFillColor(ROOT.kRed)
    h3Had.Write()
    h3Lep = ROOT.TH3F("h3Lep_"+name, "Leptons", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, minCommon, maxCommon)
    h3Lep.SetLineColor(ROOT.kMagenta)
    h3Lep.SetLineWidth(2)
    h3Lep.SetFillStyle(3006)
    h3Lep.SetFillColor(ROOT.kMagenta)
    h3Lep.Write()

    trAG.Draw(definition+":-McZDir:McLogEnergy>>"+h3Gam.GetName(), strCutBase+" && "+strCutCalOnly, "goff")
    print h3Gam.GetEntries(), "events filled into", h3Gam.GetName()
    trBKG.Draw(definition+":-McZDir:McLogEnergy>>"+h3Had.GetName(), strCutBase+" && "+strCutCalOnly+" && "+strCutHad, "goff")
    print h3Had.GetEntries(), "events filled into", h3Had.GetName()
    trBKG.Draw(definition+":-McZDir:McLogEnergy>>"+h3Lep.GetName(), strCutBase+" && "+strCutCalOnly+" && "+strCutLep, "goff")
    print h3Lep.GetEntries(), "events filled into", h3Lep.GetName()
    
    aHsPlot = []
    aHtgGam = []
    aHtgHad = []
    aHtgLep = []
    ipad = 0
    entoropy_had = 1
    entoropy_lep = 1
    for iE in range(NBIN_ENE):
        print "Energy bin No.", iE+1, " McLogEnergy:", h3Gam.GetXaxis().GetBinLowEdge(iE+1), " - ", h3Gam.GetXaxis().GetBinUpEdge(iE+1)
        aHsPlot.append([])
        aHtgGam.append([])
        aHtgHad.append([])
        aHtgLep.append([])
        for iC in range(NBIN_CTH):
            print "  Inclination bin No.", iC+1, "-McZDir:", h3Gam.GetYaxis().GetBinLowEdge(iC+1), " - ", h3Gam.GetYaxis().GetBinUpEdge(iC+1)
            aHsPlot[-1].append(ROOT.THStack("hs{0}_{1}_{2}".format(name, iE+1, iC+1), "{0:1.2f}<=McLogEnergy<{1:1.2f} and {2:1.2f}<=-McZDir<{3:1.2f}".format(h3Gam.GetXaxis().GetBinLowEdge(iE+1), h3Gam.GetXaxis().GetBinUpEdge(iE+1), h3Gam.GetYaxis().GetBinLowEdge(iC+1), h3Gam.GetYaxis().GetBinUpEdge(iC+1))))
            aHtgGam[-1].append(h3Gam.ProjectionZ("{0}_projZ_{1}_{2}".format(h3Gam.GetName(), iE+1, iC+1), iE+1, iE+1, iC+1, iC+1))
            print '    Integral of {0}: {1}'.format(aHtgGam[-1][-1].GetName(), aHtgGam[-1][-1].Integral())
            if aHtgGam[-1][-1].Integral()>0:
                aHtgGam[-1][-1].Scale(1./aHtgGam[-1][-1].Integral(1, aHtgGam[-1][-1].GetNbinsX()))
                aHsPlot[-1][-1].Add(aHtgGam[-1][-1])
                print '    Integral of {0}: {1}'.format(aHtgGam[-1][-1].GetName(), aHtgGam[-1][-1].Integral())
            else:
                print '    Integral of {0} is zero!'.format(aHtgGam[-1][-1].GetName())
            aHtgHad[-1].append(h3Had.ProjectionZ("{0}_projZ_{1}_{2}".format(h3Had.GetName(), iE+1, iC+1), iE+1, iE+1, iC+1, iC+1))
            print '    Integral of {0}: {1}'.format(aHtgHad[-1][-1].GetName(), aHtgHad[-1][-1].Integral())
            if aHtgHad[-1][-1].Integral()>0:
                aHtgHad[-1][-1].Scale(1./aHtgHad[-1][-1].Integral(1, aHtgGam[-1][-1].GetNbinsX()))
                aHsPlot[-1][-1].Add(aHtgHad[-1][-1])
                print '    Integral of {0}: {1}'.format(aHtgHad[-1][-1].GetName(), aHtgHad[-1][-1].Integral())
            else:
                print '    Integral of {0} is zero!'.format(aHtgHad[-1][-1].GetName())
            aHtgLep[-1].append(h3Lep.ProjectionZ("{0}_projZ_{1}_{2}".format(h3Lep.GetName(), iE+1, iC+1), iE+1, iE+1, iC+1, iC+1))
            print '    Integral of {0}: {1}'.format(aHtgLep[-1][-1].GetName(), aHtgLep[-1][-1].Integral())
            if aHtgLep[-1][-1].Integral()>0:
                aHtgLep[-1][-1].Scale(1./aHtgLep[-1][-1].Integral(1, aHtgGam[-1][-1].GetNbinsX()))
                aHsPlot[-1][-1].Add(aHtgLep[-1][-1])
                print '    Integral of {0}: {1}'.format(aHtgLep[-1][-1].GetName(), aHtgLep[-1][-1].Integral())
            else:
                print 'Integral of {0} is zero!'.format(aHtgLep[-1][-1].GetName())

            # Calculate Minimum Entropy
            # entropy_had = 1000
            # entropy_lep = 1000
            # var_min_ent_had = -1
            # var_min_ent_lep = -1
            # for ibvar in range(1, NBIN_VAR+1):
            #     mevt_gam = aHtgGam[-1][-1].Integral(ibvar, NBIN_VAR)
            #     mevt_had = aHtgHad[-1][-1].Integral(ibvar, NBIN_VAR)
            #     mevt_lep = aHtgLep[-1][-1].Integral(ibvar, NBIN_VAR)
            #     if mevt_had>0:
            #         frac_had = 1 / (mevt_gam/mevt_had+1)
            #         entropy_had_temp = - frac_had*math.log(frac_had) - (1.-frac_had)*math.log(1.-frac_had)
            #         if entropy_had_temp<entropy_had:
            #             entropy_had = entropy_had_temp
            #             var_min_ent_had = aHtgHad[-1][-1].GetZaxis().GetBinLowEdge(ibvar)
            #     if mevt_lep>0:
            #         frac_lep = 1 / (mevt_gam/mevt_lep+1)
            #         entropy_lep_temp = - frac_lep*math.log(frac_lep) - (1.-frac_lep)*math.log(1.-frac_lep)
            #         if entropy_lep_temp<entropy_lep:
            #             entropy_lep = entropy_lep_temp
            #             var_min_ent_lep = aHtgLep[-1][-1].GetZaxis().GetBinLowEdge(ibvar)

            ipad = iE+7*iC+1#ipad + 1
            cPlot.cd(ipad)
            aHsPlot[-1][-1].Draw("nostack")
            if logy==True:
                ROOT.gPad.SetLogy()
            # latex = ROOT.TLatex()
            # latex.SetTextSize(0.025)
            # latex.DrawLatex(.3,.8,"H_{{had}} = {0:1.2f} ({1}={2:1.1f})".format(entropy_had, definition, var_min_ent_had))
            # latex.DrawLatex(.3,.7,"H_{{lep}} = {0:1.2f} ({1}={2:1.1f})".format(entropy_lep, definition, var_min_ent_lep))
    cPlot.Write()

if __name__ == '__main__':
    main()
