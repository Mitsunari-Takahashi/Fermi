#!/usr/bin/env python

import sys
import os.path
import ROOT
from ROOT import TTree, TChain, TH1, TH2, TH3
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
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
from pMETandMJD import *
from pColor import *


@click.command()
@click.argument('pathfilelt', type=str)
@click.argument('pathfileroc', type=str)
@click.argument('namehtglt', type=str)
@click.option('--suffix', '-s', default='')
@click.option('--outdir', '-o', default='.')
@click.option('--zmax', '-z', default=90, type=float)
#@click.option('--upward', '-u', is_flag=True)
@click.option('--exclude', default=None)
def main(pathfilelt, pathfileroc, namehtglt, suffix, outdir, zmax, exclude):
    file_lt = ROOT.TFile(pathfilelt, 'READ')
    print file_lt.GetName(), "is opened."
    htg_lt = file_lt.Get(namehtglt)
    print htg_lt.GetName(), "is found."

    if exclude is not None:
        print 'Total livetime:',htg_lt.Integral()
        file_lt_exclude = ROOT.TFile(exclude, 'READ')
        print file_lt_exclude.GetName(), "is opened."
        htg_lt_exclude = file_lt_exclude.Get(namehtglt)
        print htg_lt_exclude.GetName(), "is found."
        htg_lt.Add(htg_lt_exclude, -1)
        pathfilelt_excluded = pathfilelt.replace('.root', '_excluded.root')
        file_lt_excluded = ROOT.TFile(pathfilelt_excluded, 'UPDATE')
        htg_lt.Write()
        print 'Livetime after subtraction:',htg_lt.Integral()

    # ROC file
    file_roc = ROOT.TFile(pathfileroc, 'READ')
    htg_acc = file_roc.Get("sig_acc_cth_cut")
    print htg_acc.GetName(), "is found."

    # Exposure file
    namefileexp = os.path.basename(pathfilelt)
    namefileexp = namefileexp.replace('Livetime', 'Exposure').replace('livetime', 'exposure').replace('.root', '_zmax{0:0>3.0f}.root'.format(zmax))
    if suffix!='':
        namefileexp = namefileexp.replace('.root', '_{0}.root'.format(suffix))
    #if upward is True:
    #    namefileexp = namefileexp.replace('.root', '_negativeZDIR.root')
    file_exp = ROOT.TFile('{0}/{1}'.format(outdir,namefileexp), 'RECREATE')
    file_exp.cd()
    
#    if upward is True:
    htg_exp = htg_acc.Clone("exp_cth_cut")
    # else:
    #     htg_exp = TH3F("exp_cth_cut", htg_acc.GetTitle(), htg_acc.GetXaxis().GetNbins(), htg_acc.GetXaxis().GetBinLowEdge(1), htg_acc.GetXaxis().GetBinLowEdge(htg_acc.GetXaxis().GetNbins()+1), htg_acc.GetYaxis().GetNbins(), htg_acc.GetYaxis().GetBinLowEdge(1), htg_acc.GetYaxis().GetBinLowEdge(htg_acc.GetYaxis().GetNbins()+1), htg_acc.GetZaxis().GetNbins(), -htg_acc.GetZaxis().GetBinLowEdge(htg_acc.GetYaxis().GetNbins()+1), -htg_acc.GetZaxis().GetBinLowEdge(1))
    htg_exp.SetTitle("Exposure [m^2 s];log10(Energy);BDT;cos#theta")
    for iz in range(1, htg_exp.GetNbinsZ()+1): # Inclination
        cthbin_lo = htg_exp.GetZaxis().GetBinLowEdge(iz)
        cthbin_hi = htg_exp.GetZaxis().GetBinLowEdge(iz+1)
        #if upward is False:
        lv = htg_lt.Integral(htg_lt.GetXaxis().FindBin(cthbin_lo), htg_lt.GetXaxis().FindBin(cthbin_hi), 1, htg_lt.GetYaxis().FindBin(zmax))
        #else:
         #   lv = htg_lt.Integral(htg_lt.GetXaxis().FindBin(-cthbin_hi), htg_lt.GetXaxis().FindBin(-cthbin_lo), htg_lt.GetYaxis().FindBin(180.-zmax), htg_lt.GetYaxis().GetNbins())
        for ix in range(1, htg_exp.GetNbinsX()+1): # Energy
            for iy in range(1, htg_exp.GetNbinsY()+1): # BDT Cut
                acc = htg_acc.GetBinContent(ix, iy, iz)
                aeff = acc / (2.*pi*htg_acc.GetZaxis().GetBinWidth(iz))
                accerr = htg_acc.GetBinError(ix, iy, iz)
                aefferr = accerr / (2.*pi*htg_acc.GetZaxis().GetBinWidth(iz))
                htg_exp.SetBinContent(ix, iy, iz, aeff*lv)
                htg_exp.SetBinError(ix, iy, iz, aefferr*lv)
    htg_exp.Write()


if __name__ == '__main__':
    main()
