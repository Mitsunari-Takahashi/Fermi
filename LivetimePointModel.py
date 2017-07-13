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
#import commands
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
#from pColor import *
from pAnalysisConfig import *
from pFindHEALPix import *
from pLivetime import *
import click


@click.command()
@click.argument('name', type=str)
@click.argument('ra', type=float)
@click.argument('dec', type=float)
#@click.argument('dec', type=str)
@click.option('--perf', type=str, default='/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_R100_perf.root')
@click.option('--start', type=float, default=239557417.)
@click.option('--stop', type=float, default=536457605.)
@click.option('--exstart', type=float, default=0)
@click.option('--exstop', type=float, default=0)
@click.option('--torigin', type=float, default=0.)
@click.option('--inclin', type=float, default=0.)
@click.option('--energy', type=float, default=0.)
@click.option('--errrad', type=float, default=0.)
@click.option('--suffix', type=str, default='')
@click.option('--truepoint', is_flag=True)
def main(name, ra, dec, perf, inclin, energy, suffix, start, stop, torigin, truepoint, errrad, exstart, exstop):

    # Spacecraft data
    pathFileScAll = "/disk/gamma/cta/store/takhsm/FermiData/spacecraft/mtakahas-AstroServer-00011-ft2-30s.fits"

    print "===================="
    
    # PSF histogram
    file_perf = ROOT.TFile(perf, 'READ')
    htg2_psf = file_perf.Get('psf_cth_q68_hist')
    FIXED_PSF_ENERGY = energy
    if FIXED_PSF_ENERGY<=0:
        print 'Scaling is variable for energy'
    else:
        print 'Scaling is fixed for energy at', FIXED_PSF_ENERGY
    FIXED_PSF_INCLIN = inclin
    if FIXED_PSF_INCLIN<=0:
        print 'Scaling is variable for inclination'
    else:
        print 'Scaling is fixed for inclination at', FIXED_PSF_INCLIN


    # ON/OFF regions
    nOff = 0;
#    NHPSIDE_ON = 512
    NHPSIDE_ON = 32
    ANG_CUT = 5.
    ANG_CUT_RAD = radians(5.)
    aHpx_array = [find_pointon_healpxs(ra, dec, ANG_CUT, nhpside=NHPSIDE_ON)]
    aCoordsRegion = [SkyCoord(ra, dec, unit="deg")]
    aCoordsPix_array = []
    aAreaPix_array = []
    aAreaPix_sum = []
    aStrRegion = []
    str_lit_unit = []
    for (iRegion, coordsRegion) in enumerate(aCoordsRegion):
        if iRegion==0:
            aStrRegion.append("ON")
        else:
            aStrRegion.append("OFF{0}".format(iRegion))
        if truepoint==True:
            aCoordsPix_array.append(coordsRegion)
            aAreaPix_array.append(0)
            print aCoordsPix_array[-1]
            str_lit_unit.append('sec')
            aAreaPix_sum.append(0)
        else:
            aCoordsPix_array.append([])
            aAreaPix_array.append([])
            for npix in aHpx_array[iRegion]:
                aAngPix = hppf.pix2ang(NHPSIDE_ON, npix)
                aCoordsPix_array[-1].append(SkyCoord(aAngPix[1], pi/2.-aAngPix[0], unit="rad"))
                area_pix = hppf.nside2pixarea(NHPSIDE_ON)
                aAreaPix_array[-1].append(area_pix)
            str_lit_unit.append('sec sr')
            print aCoordsPix_array[-1]
            aAreaPix_sum.append(sum(aAreaPix_array[-1]))
            print 'Solid angle =', aAreaPix_sum[-1], 'sr'

    # Output objects
    aFileToI = []
    if suffix!="":
        suffix = "_" + suffix
    if truepoint==True:
        suffix = suffix + "_TruePoint"
    
    fileRoot = ROOT.TFile("Livetime_{0}_MET{1}-{2}{3}.root".format(name, int(start), int(stop), suffix), "update")
    aHtgLt = []
    NBIN_CTH = 50
    EDGE_CTH_LOW = 0.0
    EDGE_CTH_UP = 1.0
    NBIN_ZEN = 180
    EDGE_ZEN_LOW = 0
    EDGE_ZEN_UP = 180
    NBIN_ENE = htg2_psf.GetNbinsX()
    EDGE_ENE_LOW =  htg2_psf.GetXaxis().GetBinLowEdge(1)
    EDGE_ENE_UP =  htg2_psf.GetXaxis().GetBinUpEdge(htg2_psf.GetNbinsX())
    if torigin==0:
        torigin = start 
    tmin_htg = 239557417
    tmax_htg = 536457605
        
    for hRegion in range(nOff+1):
        aHtgLt.append(ROOT.TH3D("htgLt_{0}".format(hRegion), "Livetime for {0} [{1}];Cos(Inclination angle);Zenith angle [deg];Time - {2} [sec]".format(name, str_lit_unit[hRegion], torigin), NBIN_CTH, EDGE_CTH_LOW, EDGE_CTH_UP, NBIN_ZEN, EDGE_ZEN_LOW, EDGE_ZEN_UP, max(10, int(tmax_htg-tmin_htg)/54000), tmin_htg, tmax_htg))#tPro, tPost))
    # make_livetime_histogram(aHtgLt, nOff+1 ,pathFileScAll, start+torigin, stop+torigin, aFileToI, aCoordsPix_array, aAreaPix_array, torigin) !!!BUG!!!
    make_livetime_histogram(aHtgLt, nOff+1 ,pathFileScAll, start, stop, aFileToI, aCoordsPix_array, aAreaPix_array, torigin, exstart, exstop)
    aHtgLt_projYX = []
    aHtgLt_scaled = []
    print 'Making output products...'
    fileRoot.cd()
    for jR in range(nOff+1):
        print 'Region', jR
        aHtgLt[jR].Write()
        aHtgLt_projYX.append(aHtgLt[jR].Project3D("yx"))
        aHtgLt_projYX[jR].Write()
        if truepoint==False:
            aHtgLt_scaled.append(ROOT.TH3D('{0}_scaled'.format(aHtgLt[jR].GetName()), '{0} scaled;{1};{2};log_{{10}}Energy [MeV]'.format(aHtgLt[jR].GetTitle(), aHtgLt_projYX[jR].GetXaxis().GetTitle(), aHtgLt_projYX[jR].GetYaxis().GetTitle()), NBIN_CTH, EDGE_CTH_LOW, EDGE_CTH_UP, NBIN_ZEN, EDGE_ZEN_LOW, EDGE_ZEN_UP, NBIN_ENE, EDGE_ENE_LOW, EDGE_ENE_UP))
            for iz in range(1, 1+NBIN_ENE):
                print '  Energy {0} - {1}'.format(aHtgLt_scaled[-1].GetZaxis().GetBinLowEdge(iz), aHtgLt_scaled[-1].GetZaxis().GetBinUpEdge(iz))
                if FIXED_PSF_ENERGY<=0:
                    nbin_ene_psf = htg2_psf.GetXaxis().FindBin(aHtgLt_scaled[-1].GetZaxis().GetBinCenter(iz))
                else:
                    nbin_ene_psf = htg2_psf.GetXaxis().FindBin(FIXED_PSF_ENERGY)
                for ix in range(1, 1+NBIN_CTH):
                    if FIXED_PSF_INCLIN<=0:
                        nbin_inc_psf = htg2_psf.GetYaxis().FindBin(aHtgLt_scaled[-1].GetXaxis().GetBinCenter(ix))
                    else:
                        nbin_inc_psf = htg2_psf.GetYaxis().FindBin(FIXED_PSF_INCLIN)
                    psf_cut_rad = radians(htg2_psf.GetBinContent(nbin_ene_psf, nbin_inc_psf) + errrad)
                    area_ratio = 2.*pi*(1.0-cos(psf_cut_rad)) / aAreaPix_sum[jR]
                    print '    Inclination {0} - {1} : Scaling factor = {2}'.format(aHtgLt_scaled[-1].GetXaxis().GetBinLowEdge(ix), aHtgLt_scaled[-1].GetXaxis().GetBinUpEdge(ix), area_ratio)
                    for iy in range(1, 1+NBIN_ZEN):
                        aHtgLt_scaled[-1].SetBinContent(ix, iy, iz, aHtgLt_projYX[jR].GetBinContent(ix, iy)*area_ratio)
            aHtgLt_scaled[-1].Write()
    print 'Livetime calculation finished.'


if __name__ == '__main__':
    main()
