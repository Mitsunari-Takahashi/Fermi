#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
import numpy as np
import xml.etree.ElementTree as ET
import datetime
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
import healpy as hp
from healpy import pixelfunc as hppf
import commands
import click
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pMETandMJD import *
from pColor import *

#ROOT.gStyle.SetPadGridX(True)
#ROOT.gStyle.SetPadGridY(True)
#ROOT.gStyle.SetPadTickX(True)
#ROOT.gStyle.SetPadTickY(True)

#from pCutBDT import cutBDT
#from pAnalysisConfig import *

@click.command()
@click.argument('pathfilelt', type=str)
@click.argument('pathfileperf', type=str)
@click.argument('namehtglt2', type=str)
@click.option('--suffix', '-s', default='')
def main(pathfilelt, pathfileperf, namehtglt2, suffix):
    file_lt = ROOT.TFile(pathfilelt, 'READ')
    htg_lt2 = file_lt.Get(namehtglt2)
    print htg_lt2.GetName(), "is found."

    pathfileexp = pathfilelt.replace('Livetime', 'Exposure')
    file_exp = ROOT.TFile(pathfileexp, 'RECREATE')

    # Performance file
    file_perf = ROOT.TFile(pathfileperf, 'READ')
    htg_acc = file_perf.Get("acc_cth_hist")
    print htg_acc.GetName(), "is found."

    print "===================="
    file_exp.cd()
    name_htg_exp2 = namehtglt2.replace('Lt', 'Exp')
    title_htg_lt2 = htg_lt.GetTitle()
    title_htg_exp2 = title_htg_lt2.replace('Livetime', 'Exposure')
    NBIN_LT2_X = htg_lt2.GetXaxis().GetNbins()
    NBIN_LT2_Y = htg_lt2.GetYaxis().GetNbins()
    NBIN_ENERGY = 7
    EDGE_ENERGY_LOW = 4.35
    EDGE_ENERGY_UP = 5.75
    if suffix!="":
        name_htg_exp2 = name_htg_exp2 + '_' + suffix
        title_htg_exp2 = title_htg_exp2 + ' (' + suffix + ')'
    htg_exp2 = ROOT.TH3D(name_htg_exp2, title_htg_exp2, NBIN_LT2_X, htg_lt2.GetXaxis().GetBinLowEdge(1), htg_lt2.GetXaxis().GetBinUpEdge(NBIN_LT2_X), NBIN_LT2_Y, htg_lt2.GetYaxis().GetBinLowEdge(1), htg_lt2.GetYaxis().GetBinUpEdge(NBIN_LT2_Y), NBIN_ENERGY, EDGE_ENERGY_LOW, EDGE_ENERGY_UP)

    for iz in range(NBIN_ENERGY):
        for ix in range(NBIN_LT2_X):
            acc = htg_acc.GetBinContent(htg_acc.GetXaxis().FindBin(htg_exp2.GetZaxis().GetBinCenter(iz+1)), htg_acc.GetYaxis().FindBin(htg_exp2.GetXaxis().GetBinCenter(ix+1)))
            for iy in range(NBIN_LT2_Y):
                htg_exp2.SetBinContent( ix+1, iy+1, iz+1,  htg_lt2.GetBinContent(ix+1, iy+1)*acc )
    htg_exp2.Write()


if __name__ == '__main__':
    main()
