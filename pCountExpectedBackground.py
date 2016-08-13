#!/usr/bin/env python

import sys
import ROOT
import pPoissonianProb
ROOT.gROOT.SetBatch()
#par = sys.argv

fileOFF = ROOT.TFile("GRB160509374_hOFF.root", "update")
hExpRatio = fileOFF.Get("htgExpZ90_GRB160509374_0")
dictHtgOFF = {"PSF95":fileOFF.Get("hOFF_CalOnlyR100_PSF95"), "PSF68":fileOFF.Get("hOFF_CalOnlyR100_PSF68")}
#fPoi = ROOT.TF1("fPoi", "TMath::Poisson(x, [0])", 0, 100)
aThreshold = [1, 2]

for kHtgOFF, vHtgOFF in dictHtgOFF.iteritems():
    print '#####  ', kHtgOFF, '  #####'
    vHtgOFF.Multiply(hExpRatio)
    hExBKG = vHtgOFF.GetCumulative(0)
    hExBKG.SetName("hExpectedBackground_CalOnlyR100_{0}".format(kHtgOFF))
    hExBKG.SetTitle("Expected background counts of CalOnlyR100 within {0} (Cumulative)".format(kHtgOFF))
    fileOFF.cd()
    hExBKG.Write()
    aHtgProb = {}
    aProb = {}
    for ith in aThreshold:
        print '-----  Threshold:', ith, 'events  -----'
        aHtgProb[ith] = hExBKG.Clone("hProbPoisson_CalOnlyR100_{0}_th{1}".format(kHtgOFF, ith))
        aHtgProb[ith].SetTitle("Background coincidence probability (Poissonian) of {0} CalOnlyR100 events within {1}".format(ith, kHtgOFF))
        for ibin in range(hExBKG.GetNbinsX()):
            print hExBKG.GetXaxis().GetBinLowEdge(ibin+1), "-", hExBKG.GetXaxis().GetBinLowEdge(ibin+2)
            if hExBKG.GetBinContent(ibin+1)>0:
                aProb[ith] = pPoissonianProb.calcProb([hExBKG.GetBinContent(ibin+1), ith])
            else:
                aProb[ith] = 0
            print aProb[ith]*100., "%"
            aHtgProb[ith].SetBinContent(ibin+1, aProb[ith])
        aHtgProb[ith].Write("", ROOT.TObject.kOverwrite)
    
    
    

