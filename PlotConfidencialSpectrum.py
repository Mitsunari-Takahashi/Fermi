#!/usr/bin/env python

import sys
#import os
#import numpy
#import yaml
#import datetime
#from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TMath
ROOT.gROOT.SetBatch()
import pColor


def PlotConfidencialSpectrum(gr2_conf, path_fileout, emin=1E4, emax=1.5E5, sigma_thresholds=[2, 1]):
    PROB_THRETHOLDS = [TMath.Prob(x**2,1) for x in sigma_thresholds]:
    fileout = ROOT.TFile(path_fileout, 'UPDATE')
    fileout.cd()
    cvs = ROOT.TCanvas('cvs', 'Confidential spectrum', 800, 500)
    cvs.SetLogx()
    cvs.SetGridx()
    cvs.SetGridy()
    cvs.cd()
    fnc_spe = ROOT.TF1('dNdE', '[0]/10000.*TMath::Power(x, [1])', 100., 1E6)
    htg2 = gr2_conf.Project("yx")
    for (ithre,pthre) in enumerate(PROB_THRETHOLDS):
        for ipx in range(1, htg2.GetXaxis().GetNbins()+1):
            for ipy in range(1, htg2.GetYaxis().GetNbins()+1):
                fnc_spe.SetParameter(0, TMath.Power(10, htg2.GetXaxis().GetBinCenter(ipx)))
                fnc_spe.SetParameter(1, htg2.GetYaxis().GetBinCenter(ipy))
                if htg2.GetBinContent(ipx, ipy)>thre:
                    fnc_spe.SetLineColor(ROOT.kBlue-10+6*ithre)
                    fnc_spe.DrawCopy('LSAME')
    cvs.Write()


@click.command()
@click.argument('filein', type=str)
@click.option('--namegr', type=str, default='htg2_unlikerfrac')
@click.option('--fileout', type=str, default=None)
def main(filein, namegr, fileout):
    file_conf = ROOT.TFile(filein)
    gr2_conf = file_conf.Get(namegr)
    gr2_conf.SetDirectory(0)
    file_conf.Close()
    if fileout==None:
        PlotConfidencialSpectrum(gr2_conf, file_conf)
    else:
        filenew = ROOT.TFile(fileout, 'UPDATE')
        PlotConfidencialSpectrum(gr2_conf, filenew)

if __name__ == '__main__':
    main()
