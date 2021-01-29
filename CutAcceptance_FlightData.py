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


def cut_acceptance_one(dir_one):
    dir_one.cd()
    str_dir = dir_one.GetName()
    print str_dir
    hist_acc = dir_one.Get('sig_acc_cth_cut_{}'.format(str_dir))
    acc_ideal = hist_acc.GetMaximum()
    print 'Ideal acceptance: {}'.format(acc_ideal)

    dict_cut = {}
    for th in [0.8, 0.9]:
        dict_cut[th] = hist_acc.GetBinLowEdge(hist_acc.FindLastBinAbove(th*acc_ideal))
    return ((dict_cut[0.8]+dict_cut[0.9])/2., (dict_cut[0.8]-dict_cut[0.9])/2.)


def cut_acceptance_flight_data(hist_covered, dir_top):#, erange=(1,28), cthrange=(1,40)):
    hist_acc_cut = hist_covered.Clone("hist_acc_cut") #TH2F("hist_acc_cut", "Acceptance cut", )
    for ix in range(1, 1+hist_covered.GetXaxis().GetNbins()): # Energy
        for iy in range(1, 1+hist_covered.GetYaxis().GetNbins()): # Inclination
            if hist_covered.GetBinContent(ix, iy)>0:
                dir_one = dir_top.GetDirectory('E{e}_Z{z}'.format(e=ix, z=iy))
                cut, cuterr = cut_acceptance_one(dir_one=dir_one)
                hist_acc_cut.SetBinContent(ix, iy, cut)
                if cut<=0:
                    hist_acc_cut.SetBinError(ix, iy, sys.maxint)
                elif cuterr>0:
                    hist_acc_cut.SetBinError(ix, iy, cuterr)
                else:
                    hist_acc_cut.SetBinError(ix, iy, cut)
            else:
                hist_acc_cut.SetBinContent(ix, iy, 0)
                hist_acc_cut.SetBinError(ix, iy, sys.maxint)
    dir_top.cd()
    hist_acc_cut.Write()


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

    #cut_acceptance_flight_data(file_roc)
    cut_acceptance_flight_data(hist_covered, file_roc)


if __name__ == '__main__':
    main()
