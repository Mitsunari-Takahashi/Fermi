#!/usr/bin/env python

import sys
import ROOT
import click
import pPoissonianProb
ROOT.gROOT.SetBatch()
#par = sys.argv


@click.command()
@click.argument('evton', type=str)
@click.argument('evtoff', type=str)
@click.argument('expon', type=str)
@click.argument('expoff', type=str)
@click.argument('evtenergy', type=float)
@click.argument('evtinclin', type=float)
@click.argument('evtzenith', type=float)
@click.option('--suffix', '-s', default='')
@click.option('--htgevton', default='htgEvt_GRB160509374_CalOnly_R100_PSF68')
@click.option('--htgevtoff', default='htgEvt_GRB160509374_CalOnly_R100_PSF68')
@click.option('--htgexpon', default='htgExp_GRB160509374_0_scaled_ON5000sec_CalOnlyR100')
@click.option('--htgexpoff', default='htgExp_GRB160509374_0_scaled_summed_PointingOFF')
@click.option('--mode', '-m', type=click.Choice(['PointingOFF', 'GalacticOFF']))
def main(evton, evtoff, expon, expoff, evtenergy, evtinclin, evtzenith, htgevton, htgevtoff, htgexpon, htgexpoff, mode, suffix):

    N_REBIN_INCLIN = 5
#    N_REBIN_ENERGY = 4

    fileEvtOn = ROOT.TFile(evton)
    htg3EvtOn = fileEvtOn.Get(htgevton)
    htg3EvtOn.Rebin3D(N_REBIN_INCLIN, 1, 4)
    print htg3EvtOn.GetName(), htg3EvtOn.GetNbinsX(), htg3EvtOn.GetNbinsY(), htg3EvtOn.GetNbinsZ()

    fileEvtOff = ROOT.TFile(evtoff)
    htg3EvtOff = fileEvtOff.Get(htgevtoff)
    htg3EvtOff.Rebin3D(N_REBIN_INCLIN, 1, 1)
    print htg3EvtOff.GetName(), htg3EvtOff.GetNbinsX(), htg3EvtOff.GetNbinsY(), htg3EvtOff.GetNbinsZ()

    fileExpOn = ROOT.TFile(expon)
    htg3ExpOn = fileExpOn.Get(htgexpon)
    htg3ExpOn.Rebin3D(N_REBIN_INCLIN, 1, 1)
    print htg3ExpOn.GetName(), htg3ExpOn.GetNbinsX(), htg3ExpOn.GetNbinsY(), htg3ExpOn.GetNbinsZ()

    fileExpOff = ROOT.TFile(expoff)
    htg3ExpOff = fileExpOff.Get(htgexpoff)
    htg3ExpOff.Rebin3D(N_REBIN_INCLIN, 1, 1)
    print htg3ExpOff.GetName(), htg3ExpOff.GetNbinsX(), htg3ExpOff.GetNbinsY(), htg3ExpOff.GetNbinsZ()        

    nameFileOut = evtoff.replace("Plot", "ExpectedBKG")
    if suffix!="":
        suffix = "_" + suffix
        nameFileOut = nameFileOut.replace(".root", "{0}.root".format(suffix))
    print nameFileOut
    fileOut = ROOT.TFile(nameFileOut, "UPDATE")
    fileOut.cd()
    htg3ExpOff.Write()

    # Integration over inclination and zenith
    DEG_ZEN_TOL = 10.
    if mode=="PointingOFF":
        nInclinLow = 1
        print htg3EvtOff.GetNbinsX()
        nInclinUp = htg3EvtOff.GetNbinsX()
        nZenithLow = htg3EvtOff.GetYaxis().FindBin(evtzenith-DEG_ZEN_TOL)
        nZenithUp = htg3EvtOff.GetYaxis().FindBin(evtzenith+DEG_ZEN_TOL)
        str_title = " within |zenith-{0}|<={1}".format(evtzenith, DEG_ZEN_TOL)
    elif mode=="GalacticOFF":
        nInclinLow = htg3EvtOff.GetXaxis().FindBin(evtinclin)
        nInclinUp = htg3EvtOff.GetXaxis().FindBin(evtinclin)
        nZenithLow = htg3EvtOff.GetYaxis().FindBin(evtzenith-DEG_ZEN_TOL)
        nZenithUp = htg3EvtOff.GetYaxis().FindBin(evtzenith+DEG_ZEN_TOL)
        str_title = " within {0}<=cos#theta<{1}, within |zenith-{2}|<={3}".format(htg3EvtOn.GetXaxis().GetBinLowEdge(nInclinLow), htg3EvtOn.GetXaxis().GetBinUpEdge(nInclinUp), evtzenith, DEG_ZEN_TOL)
    print "Inclination:", nInclinLow, "-", nInclinUp
    print "Zenith:", nZenithLow, "-", nZenithUp
    # htg1EvtOn = htg3EvtOn.ProjectionZ("{0}_projE_{1}".format(htg3EvtOn.GetName(), mode), nInclinLow, nInclinUp, nZenithLow, nZenithUp)
    # htg1EvtOn.SetTitle("{0} {1}".format(htg3EvtOn.GetTitle(), str_title))
    # htg1EvtOn.Write()
    htg1EvtOff = htg3EvtOff.ProjectionZ("{0}_projE_{1}".format(htg3EvtOff.GetName(), mode), nInclinLow, nInclinUp, nZenithLow, nZenithUp)
    htg1EvtOff.SetTitle("{0} {1}".format(htg3EvtOff.GetTitle(), str_title))
    htg1EvtOff.Write()
    htg1ExpOn = htg3ExpOn.ProjectionZ("{0}_projE_{1}".format(htg3ExpOn.GetName(), mode), nInclinLow, nInclinUp, nZenithLow, nZenithUp)
    htg1ExpOn.SetTitle("{0} {1}".format(htg3ExpOn.GetTitle(), str_title))
    htg1ExpOn.Write()
    htg1ExpOff = htg3ExpOff.ProjectionZ("{0}_projE_{1}".format(htg3ExpOff.GetName(), mode), nInclinLow, nInclinUp, nZenithLow, nZenithUp)
    htg1ExpOff.SetTitle("{0} {1}".format(htg3ExpOff.GetTitle(), str_title))
    htg1ExpOff.Write()

    htg1ExpRatio = htg1ExpOn.Clone("htgExpRatio_{0}".format(mode))
    htg1ExpRatio.SetTitle("Exposure ratio of ON/OFF")
    htg1ExpRatio.Divide(htg1ExpOff)
    htg1ExpRatio.Write()

    nameHtg1ExBKG = htgevton.replace("Evt", "ExBKG")
    titleHtg1ExBKG = htg3EvtOn.GetTitle()
    titleHtg1ExBKG.replace("events", "expected background events")
    htg1ExBKG = htg1EvtOff.Clone(nameHtg1ExBKG)
    htg1ExBKG.SetTitle(titleHtg1ExBKG)
#    htg3ExBKG.Rebin3D(1, 1, 4)
#    print htg3ExBKG.GetNbinsX(), htg3ExBKG.GetNbinsY(), htg3ExBKG.GetNbinsZ()
#    print htg3ExpRatio.GetNbinsX(), htg3ExpRatio.GetNbinsY(), htg3ExpRatio.GetNbinsZ()
    htg1ExBKG.Multiply(htg1ExpRatio)
    htg1ExBKG.Write()

    # print htg1ExBKG
    # htg1ExBKG_cum = htg1ExBKG.GetCumulative(ROOT.kFALSE)
    htg1ExBKG_cum = htg1ExBKG.Clone("{0}_cum".format(htg1ExBKG.GetName()))
    htg1ExBKG_cum.SetTitle("{0} (Cumulative)".format(htg1ExBKG.GetTitle()))
    for jx in range(htg1ExBKG_cum.GetNbinsX()-1, 0, -1):
        htg1ExBKG_cum.SetBinContent(jx, htg1ExBKG_cum.GetBinContent(jx)+htg1ExBKG_cum.GetBinContent(jx+1))
    htg1ExBKG_cum.Write()

    NTHRESHOLD = 1
    htgProb = htg1ExBKG_cum.Clone("htgProb_{0}".format(mode))
    htgProb.SetTitle("Poissonian coincidence probability of {0}".format(htg1ExBKG_cum.GetTitle()))
    for ibin in range(htg1ExBKG_cum.GetNbinsX()):
        print "log10(Energy) >", htg1ExBKG_cum.GetXaxis().GetBinLowEdge(ibin+1)
        if htg1ExBKG_cum.GetBinContent(ibin+1)>0:
            prob = pPoissonianProb.calcProb([htg1ExBKG_cum.GetBinContent(ibin+1), NTHRESHOLD])
        else:
            prob = 0
        print prob*100., "%"
        htgProb.SetBinContent(ibin+1, prob)
    htgProb.Write("", ROOT.TObject.kOverwrite)
    

if __name__ == '__main__':
    main()  
    

