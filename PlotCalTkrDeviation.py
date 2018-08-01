#!/usr/bin/env python

import sys
import os
import os.path
import numpy as np
#import yaml
#import datetime
#from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
ROOT.gROOT.SetBatch()
from pColor import *
#from pMETandMJD import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL

##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


def find_quantiles(hist_angdev, list_p=[0.68, 0.95, 0.99]):
    p = np.array(list_p)
    q = np.zeros_like(p)
    list_qhist = [ ROOT.TH2D('qhist_p{0:2.0f}'.format(p*100), 'Angular deviation ({0:2.0f}) [deg]'.format(p), hist_angdev.GetZaxis().GetNbins(), hist_angdev.GetZaxis().GetXbins().GetArray(), hist_angdev.GetYaxis().GetNbins(), hist_angdev.GetYaxis().GetXbins().GetArray()) for p in list_p ]
    for jY in range(1, hist_angdev.GetYaxis().GetNbins()+1):
        for kZ in range(1, hist_angdev.GetZaxis().GetNbins()+1):
            hist_pAngDev = hist_angdev.ProjectionX("{0}_pAngDev_Cth{1}_Ene{2}".format(hist_angdev.GetName(), jY, kZ), jY, jY, kZ, kZ)
            hist_pAngDev.SetTitle('{0:1.2f}<=Cos#theta<{1:1.2f} and {2:1.2f}<=logEnergy<{3:1.2f}'.format(hist_pAngDev.GetYaxis().GetBinLowEdge(jY), hist_pAngDev.GetYaxis().GetBinUpEdge(jY), hist_pAngDev.GetZaxis().GetBinLowEdge(kZ), hist_pAngDev.GetZaxis().GetBinUpEdge(kZ)))
            hist_pAngDev.GetQuantiles(len(p), q, p)
            for qhist in list_qhist:
                qhist.SetBinContent(kZ, jY, q)
    return list_qhist
            

@click.command()
@click.argument('infile', type=str)
#@click.argument('dst', nargs=-1)
@click.option('--inhist', '-i', type=str, default='h')
@click.option('--outdir', '-o', type=str, default='.')
@click.option('--suffix', '-s', type=str, default='')
#@click.option('--values', type=(str, int))
#@click.option('--values', multiple=True)
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
#@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(infile, inhist, outdir, suffix):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    infile = TFile(infile, "READ")
    hist_angdev = infile.Get(inhist)
    list_qhist = find_quantiles(hist_angdev)
    outname = os.path.basename(infile).replace('.root', '_quantiles.root')
    outpath = '/'.join([outdir, outname])
    outfile  = ROOT.TFile(outpath, "RECREATE")
    outfile.cd()
    for qh in list_qhist:
        qh.Write()
    

if __name__ == '__main__':
    main()
