#!/usr/bin/env python

import sys
import os.path
import ROOT
from ROOT import TFile, TTree, TChain, TH1, TH2, TH3, TH1F, TH2F, TH3F, TGraph, TGraphErrors, TCanvas
import numpy as np
import commands
import click
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, log, log10


ROOT.gROOT.SetBatch()


def cut_bdt_one(dir_one):
    dir_one.cd()
    str_dir = dir_one.GetName()
    print str_dir
    graph_evtdens = dir_one.Get('graph_evtdens_{}'.format(str_dir))
    graph_cut = dir_one.Get('graph_cut_{}'.format(str_dir))
    hist_acc = dir_one.Get('sig_acc_cth_cut_{}'.format(str_dir))
    acc_ideal = hist_acc.GetMaximum()
    print 'Ideal acceptance: {}'.format(acc_ideal)

    hist_dist = TH1F('hist_dist', 'Distribution of the resudual event density for cut', 100, 0, 1.1*graph_evtdens.GetMaximum())

    d1, d2 = ROOT.Double(0), ROOT.Double(0)
    for ip in range(graph_evtdens.GetN()):
        graph_evtdens.GetPoint( ip, d1, d2 )
        hist_dist.Fill(d2)
    hist_dist.Fit("gaus")
    hist_dist.Write()
    mean = hist_dist.GetFunction("gaus").GetParameter(1)
    sigma = hist_dist.GetFunction("gaus").GetParameter(2)
    print 'Residual event density plateau: {m} +/- {s}'.format(m=mean, s=sigma)
    thresholds = [mean + 1.*sigma, mean + 2.*sigma]
    acc_allowed = {}
    cut_allowed = {}

    evt1, evt2 = ROOT.Double(0), ROOT.Double(0)
    cut1, cut2 = ROOT.Double(0), ROOT.Double(0)

    jp=1000
    graph_evtdens.GetPoint( jp, evt1, evt2 )
    evterr = graph_evtdens.GetErrorY( jp )
    graph_cut.GetPoint( jp, cut1, cut2 )
    print '{0}: {1} +/- {2}'.format(jp, evt2, evterr)

    for threshold in thresholds:
        print 'Threshold: {}'.format(threshold)
        acc_allowed[threshold] = []
        cut_allowed[threshold] = []
        for jp in range(graph_evtdens.GetN()):
            graph_evtdens.GetPoint( jp, evt1, evt2 )
            evterr = graph_evtdens.GetErrorY( jp )
            graph_cut.GetPoint( jp, cut1, cut2 )
            #print '{0}: {1} +/- {2}'.format(jp, evt2, evterr)
            if evt2-evterr<=threshold and evt2+evterr>=threshold:
                acc_allowed[threshold].append(evt1)
                cut_allowed[threshold].append(cut2)
        print 'Cut: {0}-{1}'.format(min(cut_allowed[threshold]), max(cut_allowed[threshold]))
        print 'Acceptance: {0}-{1}'.format(min(acc_allowed[threshold]), max(acc_allowed[threshold]))
                


def cut_bdt_flight_data(hist_covered, dir_top):
    for ix in range(1, 1+hist_covered.GetXaxis().GetNbins()): # Energy
        for iy in range(1, 1+hist_covered.GetYaxis().GetNbins()): # Inclination
            if hist_covered.GetBinContent(ix, iy)>0:
                dir_one = dir_top.GetDirectory('E{e}_Z{z}'.format(e=ix, z=iy))
                cut_bdt_one(dir_one=dir_one)



@click.command()
@click.argument('roc3d', type=str)
#@click.option('--outpath', '-o', default=None)
def main(roc3d):
    # Input
    file_roc = TFile(roc3d, "UPDATE")
    hist_covered = file_roc.Get("hist_covered")

    # Output
    #file_out = TFile(outpath, "RECREATE")
    #file_out.cd()

    cut_bdt_flight_data(hist_covered, file_roc)


if __name__ == '__main__':
    main()
