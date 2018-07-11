#!/usr/bin/env python

import ROOT
from ROOT import gROOT, TDirectoryFile, TH3, TH3D, TH2D, TH1D
import sys
import math
from math import log10
import click
import subprocess
from pLsList import ls_list
from ctypes import *
gROOT.SetBatch()


@click.command()
@click.argument('name', type=str)
@click.argument('definition', type=str)
@click.option('--ag', '-g', type=str, default='/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2016Jun.root', help="Path of AG event file")
@click.option('--bkg', '-b', type=str, default='/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_Combined.root', help="Path of BKG event file")
@click.option('--logy', is_flag=True, help="Plot y-axis in logarithmic scale")
@click.option('--nbine', '-e', type=int, default=7, help="Number of energy bin")
@click.option('--nbincth', '-c', type=int, default=2, help="Number of cos(theta) bin")
@click.option('--vmin', type=float, default=0, help="Lower edge of the histogram")
@click.option('--vmax', type=float, default=10, help="Lower edge of the histogram")
@click.option('--outfile', '-o', type=str, default='./MVA_Variables_CalOnly', help="Name base of output file.")
@click.option('--force', '-f', is_flag=True, help="Force to resample values")
def main(name, definition, ag, bkg, nbine, nbincth, logy, vmin, vmax, outfile, force):
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

    cPlot = ROOT.TCanvas(name, definition, int(ncx*200), int(ncy*200))
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

    fileOut = ROOT.TFile("{0}.root".format(outfile), "UPDATE")
    print 'Output file', fileOut.GetName()
    dirOut = fileOut.GetDirectory(name)
    if dirOut == None:
        dirOut = fileOut.mkdir(name)
    print 'Output directory', dirOut.GetName()
    fileOut.cd(dirOut.GetName())

    h3Gam = dirOut.Get("h3Gam_"+name)
    if force is True or h3Gam == None:
        del h3Gam
        h3Gam = TH3D("h3Gam_"+name, "Gammas", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, vmin, vmax)
        print h3Gam.GetName(), h3Gam.GetXaxis().GetNbins(), h3Gam.GetXaxis().GetBinLowEdge(1), h3Gam.GetXaxis().GetBinUpEdge(h3Gam.GetXaxis().GetNbins()), h3Gam.GetYaxis().GetNbins(), h3Gam.GetYaxis().GetBinLowEdge(1), h3Gam.GetYaxis().GetBinUpEdge(h3Gam.GetYaxis().GetNbins()), h3Gam.GetZaxis().GetNbins(), h3Gam.GetZaxis().GetBinLowEdge(1), h3Gam.GetZaxis().GetBinUpEdge(h3Gam.GetZaxis().GetNbins())
        trAG.Draw("{0}:-McZDir:McLogEnergy>>{1}".format(definition, h3Gam.GetName()), " && ".join([strCutCalOnly,strCutBase]), "goff")
        h3Gam.SetLineColor(ROOT.kBlue)
        h3Gam.SetLineWidth(2)
        h3Gam.SetFillStyle(3004)
        h3Gam.SetFillColor(ROOT.kBlue)
        h3Gam.Write()
    print h3Gam.GetEntries(), "events filled into", h3Gam.GetName(), h3Gam.GetXaxis().GetNbins(), h3Gam.GetXaxis().GetBinLowEdge(1), h3Gam.GetXaxis().GetBinUpEdge(h3Gam.GetXaxis().GetNbins()), h3Gam.GetYaxis().GetNbins(), h3Gam.GetYaxis().GetBinLowEdge(1), h3Gam.GetYaxis().GetBinUpEdge(h3Gam.GetYaxis().GetNbins()), h3Gam.GetZaxis().GetNbins(), h3Gam.GetZaxis().GetBinLowEdge(1), h3Gam.GetZaxis().GetBinUpEdge(h3Gam.GetZaxis().GetNbins())

    h3Had = dirOut.Get("h3Had_"+name)
    if force is True or h3Had == None:
        del h3Had
        h3Had = TH3D("h3Had_"+name, "Hadrons", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, vmin, vmax)
        print h3Had.GetName(), h3Had.GetXaxis().GetNbins(), h3Had.GetXaxis().GetBinLowEdge(1), h3Had.GetXaxis().GetBinUpEdge(h3Had.GetXaxis().GetNbins()), h3Had.GetYaxis().GetNbins(), h3Had.GetYaxis().GetBinLowEdge(1), h3Had.GetYaxis().GetBinUpEdge(h3Had.GetYaxis().GetNbins()), h3Had.GetZaxis().GetNbins(), h3Had.GetZaxis().GetBinLowEdge(1), h3Had.GetZaxis().GetBinUpEdge(h3Had.GetZaxis().GetNbins())
        h3Had.SetLineColor(ROOT.kRed)
        h3Had.SetLineWidth(2)
        h3Had.SetFillStyle(3005)
        h3Had.SetFillColor(ROOT.kRed)
        trBKG.Draw(definition+":-McZDir:McLogEnergy>>"+h3Had.GetName(), " && ".join([strCutBase, strCutCalOnly, strCutHad]), "goff")
        h3Had.Write()
    print h3Had.GetEntries(), "events filled into", h3Had.GetName(), h3Had.GetXaxis().GetNbins(), h3Had.GetXaxis().GetBinLowEdge(1), h3Had.GetXaxis().GetBinUpEdge(h3Had.GetXaxis().GetNbins()), h3Had.GetYaxis().GetNbins(), h3Had.GetYaxis().GetBinLowEdge(1), h3Had.GetYaxis().GetBinUpEdge(h3Had.GetYaxis().GetNbins()), h3Had.GetZaxis().GetNbins(), h3Had.GetZaxis().GetBinLowEdge(1), h3Had.GetZaxis().GetBinUpEdge(h3Had.GetZaxis().GetNbins())

    h3Lep = dirOut.Get("h3Lep_"+name)
    if force is True or h3Lep == None:
        del h3Lep
        h3Lep = TH3D("h3Lep_"+name, "Leptons", NBIN_ENE, LOWEDGE_ENE, UPEDGE_ENE, NBIN_CTH, LOWEDGE_CTH, UPEDGE_CTH, NBIN_VAR, vmin, vmax)
        h3Lep.SetLineColor(ROOT.kMagenta)
        h3Lep.SetLineWidth(2)
        h3Lep.SetFillStyle(3006)
        h3Lep.SetFillColor(ROOT.kMagenta)
        trBKG.Draw(definition+":-McZDir:McLogEnergy>>"+h3Lep.GetName(), " && ".join([strCutBase, strCutCalOnly, strCutLep]), "goff")
        h3Lep.Write()
    print h3Lep.GetEntries(), "events filled into", h3Lep.GetName(), h3Lep.GetXaxis().GetNbins(), h3Lep.GetXaxis().GetBinLowEdge(1), h3Lep.GetXaxis().GetBinUpEdge(h3Lep.GetXaxis().GetNbins()), h3Lep.GetYaxis().GetNbins(), h3Lep.GetYaxis().GetBinLowEdge(1), h3Lep.GetYaxis().GetBinUpEdge(h3Lep.GetYaxis().GetNbins()), h3Lep.GetZaxis().GetNbins(), h3Lep.GetZaxis().GetBinLowEdge(1), h3Lep.GetZaxis().GetBinUpEdge(h3Lep.GetZaxis().GetNbins())
    
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
            aHtgGam[-1][-1].SetLineColor(ROOT.kBlue)
            aHtgGam[-1][-1].SetLineWidth(2)
            aHtgGam[-1][-1].SetFillStyle(3004)
            aHtgGam[-1][-1].SetFillColor(ROOT.kBlue)
            if aHtgGam[-1][-1].Integral()>0:
                aHtgGam[-1][-1].Scale(1./aHtgGam[-1][-1].Integral(1, aHtgGam[-1][-1].GetNbinsX())/aHtgGam[-1][-1].GetXaxis().GetBinWidth(1))
                aHtgGam[-1][-1].Write()
                aHsPlot[-1][-1].Add(aHtgGam[-1][-1])
                print '    Integral of {0}: {1}'.format(aHtgGam[-1][-1].GetName(), aHtgGam[-1][-1].Integral())
            else:
                print '    Integral of {0} is zero!'.format(aHtgGam[-1][-1].GetName())
            aHtgHad[-1].append(h3Had.ProjectionZ("{0}_projZ_{1}_{2}".format(h3Had.GetName(), iE+1, iC+1), iE+1, iE+1, iC+1, iC+1))
            print '    Integral of {0}: {1}'.format(aHtgHad[-1][-1].GetName(), aHtgHad[-1][-1].Integral())
            aHtgHad[-1][-1].SetLineColor(ROOT.kRed)
            aHtgHad[-1][-1].SetLineWidth(2)
            aHtgHad[-1][-1].SetFillStyle(3005)
            aHtgHad[-1][-1].SetFillColor(ROOT.kRed)
            if aHtgHad[-1][-1].Integral()>0:
                aHtgHad[-1][-1].Scale(1./aHtgHad[-1][-1].Integral(1, aHtgGam[-1][-1].GetNbinsX())/aHtgHad[-1][-1].GetXaxis().GetBinWidth(1))
                aHtgHad[-1][-1].Write()
                aHsPlot[-1][-1].Add(aHtgHad[-1][-1])
                print '    Integral of {0}: {1}'.format(aHtgHad[-1][-1].GetName(), aHtgHad[-1][-1].Integral())
            else:
                print '    Integral of {0} is zero!'.format(aHtgHad[-1][-1].GetName())
            aHtgLep[-1].append(h3Lep.ProjectionZ("{0}_projZ_{1}_{2}".format(h3Lep.GetName(), iE+1, iC+1), iE+1, iE+1, iC+1, iC+1))
            print '    Integral of {0}: {1}'.format(aHtgLep[-1][-1].GetName(), aHtgLep[-1][-1].Integral())
            aHtgLep[-1][-1].SetLineColor(ROOT.kMagenta)
            aHtgLep[-1][-1].SetLineWidth(2)
            aHtgLep[-1][-1].SetFillStyle(3006)
            aHtgLep[-1][-1].SetFillColor(ROOT.kMagenta)
            if aHtgLep[-1][-1].Integral()>0:
                aHtgLep[-1][-1].Scale(1./aHtgLep[-1][-1].Integral(1, aHtgGam[-1][-1].GetNbinsX())/aHtgLep[-1][-1].GetXaxis().GetBinWidth(1))
                aHtgLep[-1][-1].Write()
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
            aHsPlot[-1][-1].Write()
            if logy==True:
                ROOT.gPad.SetLogy()
            # latex = ROOT.TLatex()
            # latex.SetTextSize(0.025)
            # latex.DrawLatex(.3,.8,"H_{{had}} = {0:1.2f} ({1}={2:1.1f})".format(entropy_had, definition, var_min_ent_had))
            # latex.DrawLatex(.3,.7,"H_{{lep}} = {0:1.2f} ({1}={2:1.1f})".format(entropy_lep, definition, var_min_ent_lep))
    cPlot.Write()

if __name__ == '__main__':
    main()
