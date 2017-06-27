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
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pMETandMJD import *
from pColor import *


@click.command()
@click.argument('pathfilelt', type=str)
@click.argument('pathfileperf', type=str)
@click.argument('namehtglt', type=str)
@click.option('--suffix', '-s', default='')
@click.option('--outdir', '-o', default='.')
@click.option('--extend3d', '-e', is_flag=True)
@click.option('--exclude', default=None)
@click.option('--limpsf', default=None, type=float)
def main(pathfilelt, pathfileperf, namehtglt, suffix, outdir, extend3d, exclude, limpsf):
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

    if extend3d==True:
        htg_ext3d = ROOT.TH3D("{0}_ext3d".format(htg_lt.GetName()), htg_lt.GetTitle(), htg_lt.GetNbinsX(), htg_lt.GetXaxis().GetXmin(), htg_lt.GetXaxis().GetXmax(), htg_lt.GetNbinsY(), htg_lt.GetYaxis().GetXmin(), htg_lt.GetYaxis().GetXmax(), 28, 4.35, 5.75)
        for ix in range(1, htg_ext3d.GetNbinsX()+1):
            for iy in range(1, htg_ext3d.GetNbinsY()+1):
                cont = htg_lt.GetBinContent(ix, iy)
                for iz in range(1, htg_ext3d.GetNbinsZ()+1):
                    htg_ext3d.SetBinContent(ix, iy, iz, cont)
        htg_lt = htg_ext3d

    for ix in range(1, htg_lt.GetNbinsX()+1):
        for iy in range(1, htg_lt.GetNbinsY()+1):
            for iz in range(1, htg_lt.GetNbinsZ()+1):
                htg_lt.SetBinError(ix, iy, iz, 0)
       
    namefileexp = os.path.basename(pathfilelt)
    namefileexp = namefileexp.replace('Livetime', 'Exposure')
    if limpsf is not None:
        if suffix=='':
            suffix = 'limpsf{1}deg'.format(suffix, int(limpsf+0.5))
        else:
            suffix = '{0}_limpsf{1}deg'.format(suffix, int(limpsf+0.5))
    if suffix!='':
        namefileexp = namefileexp.replace('.root', '_{0}.root'.format(suffix))
    #pathfileexp = pathfilelt.replace('Livetime', 'Exposure')
    #pathfileexp = pathfileexp.replace('.root', '_{0}.root'.format(suffix))
    file_exp = ROOT.TFile('{0}/{1}'.format(outdir,namefileexp), 'RECREATE')

    # Performance file
    file_perf = ROOT.TFile(pathfileperf, 'READ')
    htg_acc = file_perf.Get("acc_cth_hist")
    print htg_acc.GetName(), "is found."
    if limpsf is not None:
        print 'PSF < {0} deg'.format(limpsf)
        htg_psf = file_perf.Get("psf_cth_q68_hist")
        print htg_psf.GetName(), "is found."
        for ix in range(1, htg_acc.GetXaxis().GetNbins()+1):
            for iy in range(1, htg_acc.GetYaxis().GetNbins()+1):
                psf_local = htg_psf.GetBinContent(htg_psf.GetXaxis().FindBin(htg_acc.GetXaxis().GetBinCenter(ix)), htg_psf.GetYaxis().FindBin(htg_acc.GetYaxis().GetBinCenter(iy)))
                if psf_local>limpsf:
                    htg_acc.SetBinContent(ix, iy, 0)

    print "===================="
    file_exp.cd()
    name_htg_exp3 = namehtglt.replace('Lt', 'Exp')
    title_htg_lt = htg_lt.GetTitle()
    title_htg_exp3 = title_htg_lt.replace('Livetime', 'Exposure')
    NBIN_LT_X = htg_lt.GetXaxis().GetNbins()
    NBIN_LT_Y = htg_lt.GetYaxis().GetNbins()
    NBIN_ENERGY = htg_lt.GetZaxis().GetNbins() #28
    EDGE_ENERGY_LOW =htg_lt.GetZaxis().GetBinLowEdge(1) # 4.35
    EDGE_ENERGY_UP = htg_lt.GetZaxis().GetBinUpEdge(NBIN_ENERGY) #5.75
    if suffix!="":
        name_htg_exp3 = name_htg_exp3 + '_' + suffix
        title_htg_exp3 = title_htg_exp3 + ' (' + suffix + ')'
    htg_exp3 = ROOT.TH3D(name_htg_exp3, title_htg_exp3, NBIN_LT_X, htg_lt.GetXaxis().GetBinLowEdge(1), htg_lt.GetXaxis().GetBinUpEdge(NBIN_LT_X), NBIN_LT_Y, htg_lt.GetYaxis().GetBinLowEdge(1), htg_lt.GetYaxis().GetBinUpEdge(NBIN_LT_Y), NBIN_ENERGY, EDGE_ENERGY_LOW, EDGE_ENERGY_UP)

    for iz in range(1, NBIN_ENERGY+1):
        for ix in range(1, NBIN_LT_X+1):
            ix_acc = htg_acc.GetXaxis().FindBin(htg_exp3.GetZaxis().GetBinCenter(iz)) #Energy
            iy_acc = htg_acc.GetYaxis().FindBin(htg_exp3.GetXaxis().GetBinCenter(ix)) #Cos(Inclination)
            acc = htg_acc.GetBinContent(ix_acc, iy_acc)
            accerr = htg_acc.GetBinError(ix_acc, iy_acc)
            for iy in range(1, NBIN_LT_Y+1):
                htg_exp3.SetBinContent( ix, iy, iz,  htg_lt.GetBinContent(ix, iy, iz)*acc*htg_exp3.GetXaxis().GetBinWidth(ix)/htg_acc.GetYaxis().GetBinWidth(iy_acc)/4./math.pi )
                htg_exp3.SetBinError( ix, iy, iz,  htg_lt.GetBinContent(ix, iy, iz)*accerr*htg_exp3.GetXaxis().GetBinWidth(ix)/htg_acc.GetYaxis().GetBinWidth(iy_acc)/4./math.pi )
    htg_exp3.Write()


if __name__ == '__main__':
    main()
