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
@click.option('--extend3d', '-e', is_flag=True)
def main(pathfilelt, pathfileperf, namehtglt2, suffix, extend3d):
    file_lt = ROOT.TFile(pathfilelt, 'READ')
    htg_lt3 = file_lt.Get(namehtglt2)
    print htg_lt3.GetName(), "is found."

    if extend3d==True:
        htg_ext3d = ROOT.TH3D("{0}_ext3d".format(htg_lt3.GetName()), htg_lt3.GetTitle(), htg_lt3.GetNbinsX(), htg_lt3.GetXaxis().GetXmin(), htg_lt3.GetXaxis().GetXmax(), htg_lt3.GetNbinsY(), htg_lt3.GetYaxis().GetXmin(), htg_lt3.GetYaxis().GetXmax(), 7, 4.35, 5.75)
        for ix in range(1, htg_ext3d.GetNbinsX()+1):
            for iy in range(1, htg_ext3d.GetNbinsY()+1):
                cont = htg_lt3.GetBinContent(ix, iy)
                for iz in range(1, htg_ext3d.GetNbinsZ()+1):
                    htg_ext3d.SetBinContent(ix, iy, iz, cont)
        htg_lt3 = htg_ext3d
       
    pathfileexp = pathfilelt.replace('Livetime', 'Exposure')
    file_exp = ROOT.TFile(pathfileexp, 'RECREATE')

    # Performance file
    file_perf = ROOT.TFile(pathfileperf, 'READ')
    htg_acc = file_perf.Get("acc_cth_hist")
    print htg_acc.GetName(), "is found."

    print "===================="
    file_exp.cd()
    name_htg_exp3 = namehtglt2.replace('Lt', 'Exp')
    title_htg_lt3 = htg_lt3.GetTitle()
    title_htg_exp3 = title_htg_lt3.replace('Livetime', 'Exposure')
    NBIN_LT2_X = htg_lt3.GetXaxis().GetNbins()
    NBIN_LT2_Y = htg_lt3.GetYaxis().GetNbins()
    NBIN_ENERGY = htg_lt3.GetZaxis().GetNbins() #7
    EDGE_ENERGY_LOW =htg_lt3.GetZaxis().GetBinLowEdge(1) # 4.35
    EDGE_ENERGY_UP = htg_lt3.GetZaxis().GetBinUpEdge(NBIN_ENERGY) #5.75
    if suffix!="":
        name_htg_exp3 = name_htg_exp3 + '_' + suffix
        title_htg_exp3 = title_htg_exp3 + ' (' + suffix + ')'
    htg_exp3 = ROOT.TH3D(name_htg_exp3, title_htg_exp3, NBIN_LT2_X, htg_lt3.GetXaxis().GetBinLowEdge(1), htg_lt3.GetXaxis().GetBinUpEdge(NBIN_LT2_X), NBIN_LT2_Y, htg_lt3.GetYaxis().GetBinLowEdge(1), htg_lt3.GetYaxis().GetBinUpEdge(NBIN_LT2_Y), NBIN_ENERGY, EDGE_ENERGY_LOW, EDGE_ENERGY_UP)

    for iz in range(1, NBIN_ENERGY+1):
        for ix in range(1, NBIN_LT2_X+1):
            ix_acc = htg_acc.GetXaxis().FindBin(htg_exp3.GetZaxis().GetBinCenter(iz)) #Energy
            iy_acc = htg_acc.GetYaxis().FindBin(htg_exp3.GetXaxis().GetBinCenter(ix)) #Cos(Inclination)
            acc = htg_acc.GetBinContent(ix_acc, iy_acc)
            for iy in range(1, NBIN_LT2_Y+1):
                htg_exp3.SetBinContent( ix, iy, iz,  htg_lt3.GetBinContent(ix, iy, iz)*acc*htg_exp3.GetXaxis().GetBinWidth(ix)/htg_acc.GetYaxis().GetBinWidth(iy_acc)/4./math.pi )
    htg_exp3.Write()


if __name__ == '__main__':
    main()
