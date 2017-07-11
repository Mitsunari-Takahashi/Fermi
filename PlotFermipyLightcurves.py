#!/usr/bin/env python

import sys
import os
import numpy as np
#from fermipy.gtanalysis import GTAnalysis
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi, log10, sqrt, ceil, isnan
import click
import copy
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
ROOT.gROOT.SetBatch()
from pColor import *
from pMETandMJD import *


MEV2ERG = 1.6021766208e-6


def plot_fermipy_lightcurves(name, pathcsv, pathoutdir, thts=25, suffix=='', formats=['png', 'tex', 'pdf']):
    """Plot the results from AnalyzeGRB_fermipy.py in root format.
"""
    if suffix!='':
        suffix = '_' + suffix
    trlc = ROOT.TTree('trlc', 'Light curve Tree')
    trlc.ReadFile(pathcsv)
    fileout = ROOT.TFile('{0}/{1}{2}.root'.format(pathoutdir, name, suffix), 'RECREATE')
    fileout.cd()
    trlc.Write()
    csummary = ROOT.TCanvas('csummary', name, 800, 800)
    csummary.Divide(1,3)
    print 'Flux light curve'
    trlc.Draw("(stop+start)/2.:(stop-start)/2.:flux:flux_err", "ts>={0}".format(thts), "goff")
    grflux = ROOT.TGraphErrors(tr.lc.GetEntries("ts>={0}".format(thts)), trlc.GetV1(), trlc.GetV3(), trlc.GetV2(), trlc.GetV4())
    grflux.SetName('grflux')
    grflux.SetTitle('Flux light curve')
    grflux.SetLineColor(akColor(0))
    grflux.SetLineStyle(1)
    grflux.SetLineWidth(2)
    grflux.SetMarkerColor(akColor(0))
    grflux.SetMarkerStyle(21)
    grflux.SetMarkerSize(0.7)
    grflux.Write()
    print 'Flux_Ul95 light curve'
    trlc.Draw("(stop+start)/2.:(stop-start)/2.:flux_ul95", "ts<{0}".format(thts), "goff")
    grflux_ul95 = ROOT.TGraphErrors(tr.lc.GetEntries("ts<{0}".format(thts)), trlc.GetV1(), trlc.GetV3(), trlc.GetV2(), 0)
    grflux_ul95.SetName('grflux_ul95')
    grflux_ul95.SetTitle('Flux_ul95 light curve')
    grflux_ul95.SetLineColor(akColor(1))
    grflux_ul95.SetLineStyle(1)
    grflux_ul95.SetLineWidth(2)
    grflux_ul95.SetMarkerColor(akColor(1))
    grflux_ul95.SetMarkerStyle(23)
    grflux_ul95.SetMarkerSize(0.7)
    grflux_ul95.Write()
    mgrfluxul = ROOT.TMultiGraph('mgrfluxul', 'Flux light curve')
    mgrfluxul.Add(grflux)
    mgrfluxul.Add(grflux_ul95)
    mgrfluxul.Write()
    csummary.cd(1)
    mgrfluxul.Draw("APL")
    print 'TS light curve'
    trlc.Draw("(stop+start)/2.:(stop-start)/2.:ts", "".format(thts), "goff")
    grts = ROOT.TGraphErrors(tr.lc.GetEntries(), trlc.GetV1(), trlc.GetV3(), trlc.GetV2(), 0)
    grts.SetName('grts')
    grts.SetTitle('TS light curve')
    grts.SetLineColor(akColor(0))
    grts.SetLineStyle(1)
    grts.SetLineWidth(2)
    grts.SetMarkerColor(akColor(0))
    grts.SetMarkerStyle(25)
    grts.SetMarkerSize(0.7)
    grts.Write()
    csummary.cd(2)
    grts.Draw("APL")
    print 'Index light curve'
    trlc.Draw("(stop+start)/2.:(stop-start)/2.:index:index_err", "ts>={0}".format(thts), "goff")
    grindex = ROOT.TGraphErrors(tr.lc.GetEntries("ts>={0}".format(thts)), trlc.GetV1(), trlc.GetV3(), trlc.GetV2(), trlc.GetV4())
    grindex.SetName('grindex')
    grindex.SetTitle('Index light curve')
    grindex.SetLineColor(akColor(0))
    grindex.SetLineStyle(1)
    grindex.SetLineWidth(2)
    grindex.SetMarkerColor(akColor(0))
    grindex.SetMarkerStyle(20)
    grindex.SetMarkerSize(0.7)
    grindex.Write()
    csummary.cd(3)
    grindex.Draw("APL")
    csummary.Write()
    for fmt in formats:
        csummary.SaveAs('{0}{1}.{2}'.format(name, suffix, fmt))


@click.command()
@click.argument('srcname', type=str)
@click.argument('inputcsv', type=str)
@click.option('outdir', type=str, default='.')
@click.option('--suffix', type=str, default='')
@click.option('--tsthreshold', type=float, default=25)
@click.option('--formats', multiple=True, default=['png', 'tex', 'pdf'])
def main(srcname, outdir):
    plot_fermipy_lightcurves(name=srcname, pathcsv=inputcsv, pathoutdir=outdir, thts=tsthreshold, suffix=='', formats=formats)


if __name__ == '__main__':
    main()
