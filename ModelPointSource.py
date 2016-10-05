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
ROOT.gROOT.SetBatch()
from array import array
import math
from math import pi, cos, sin, tan, acos, asin, atan, radians, degrees
#from pColor import *
from pAnalysisConfig import *
from pFindHEALPix import *
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")


class TrueSource: 
    def __init__(self, name_src, spatial_type="Point-like", htg_sed=""):
        self.NAME = name_src
        self.SPATIAL_TYPE = spatial_type
        self.NBIN_ENERGY = 7
        self.EDGE_ENERGY_LOW = 4.35
        self.EDGE_ENERGY_UP = 5.75
        if htg_sed=="":
            htg_sed = ROOT.TH1D('htg_sed_test', 'SED just for test', self.NBIN_ENERGY, self.EDGE_ENERGY_LOW, self.EDGE_ENERGY_UP)
            for hbin in range(1, htg_sed.GetNbinsX()+1):
                htg_sed.SetBinContent(hbin, hbin)
        self.HTG_SED = htg_sed # [ Counts / m^2 / sec ]
        self.dict_htg2_model = {}


class TruePointSource(TrueSource):
    def __init__(self, name_src, htg_sed, ra, dec):
        TrueSource.__init__(self, name_src, "Point-like", htg_sed)
        self.RA = ra
        self.DEC = dec
        self.RA_RAD = radians(self.RA)
        self.DEC_RAD = radians(self.DEC)
        

    def model(self, TP_HTG_KING, HTG2_LIVETIME, HTG2_ACCEPTANCE, NHPSIDE=512, THRESHOLD_ANGDIST=15):
        """Model the point source with the PSF of King function. Energy dispersion is ignored currently.
        """
        PIX_TRUE = hppf.ang2pix(NHPSIDE, pi/2.-self.DEC_RAD, self.RA_RAD)
        NPIX = nside2npix(NHPSIDE)
        coo_src = SkyCoord(self.RA, self.DEC, unit="deg")
        HTG1_LT = HTG2_LIVETIME.ProjectionX("{0}_projTheta".format(HTG2_LIVETIME.GetName()))
        htg2_model = ROOT.TH2D("htg2_model_{0}".format(NHPSIDE), "Model of {0} (NSIDE={1})".format(self.NAME, NHPSIDE), NPIX, 0, NPIX, self.NBIN_ENERGY, self.EDGE_NEERGY_LOW, self.EDGE_ENERGY_UP)
        if TP_HTG_KING[0].GetNbinsX()!=self.NBIN_ENERGY:
            print "Number of energy bins is not matched."
            return 1
        # PSF
        fc_King = ROOT.TF1("fc_King", "[0]*(1.-1./[2])*pow(1.+(x/[1])^2/2./[2],-[2])/2./TMath::Pi()/[1]^2", 0, pi)

        for iEne in range(1, TP_HTG_KING[0].GetNbinsX()+1):
            flux_true_diff = self.HTG_SED.GetBinContent(iEne)
            for iTh in range(1, HTG1_LT.GetNbinsX()+1):
                for ipar in range(3): # Setting the parameters of King function
                    # PSF
                    fc_King.FixParameter(ipar, TP_HTG_KING[ipar].GetBinContent(iEne, TP_HTG_KING[ipar].GetYaxis().FindBin(HTG1_LT.GetBinCenter(iTh))))
                # Acceptance
                scale_acc = HTG2_ACCEPTANCE.GetBinContent(HTG2_ACCEPTANCE.GetXaxis().FindBin(TP_HTG_KING[0].GetXaxis().GetBinCenter(iEne), HTG2_ACCEPTANCE.GetYaxis().FindBin(HTG1_LT.GetXaxis().GetBinCenter(iTh))))
                scale_exp = scale_acc * HTG1_LT.GetBinContent(iTh) # Integrated exposure value for a certain energy and inclination angle
                for ipix in range(NPIX):
                    tp_pix = hppf.pix2ang(NHPSIDE, ipix)
                    coo_pix = SkyCoord(tp_pix[1], pi/2.-tp_pix[0], unit="rad")
                    ang_pix2src = coo_src.separation(coo_pix)
                    deg_pix2src = float(ang_pix2src.to_string(unit=u.deg, decimal=True))
                    if deg_pix2src < THRESHOLD_ANGDIST:
                        rad_pix2src = float(ang_pix2src.to_string(unit=u.rad, decimal=True)) # Distance between the source and each pixel in radians
                        scale_psf = fc_King.Eval(rad_pix2src)*hppf.nside2pixarea(NHPSIDE)
                        htg2_model.Fill(ipix, htg2_model.GetBinCenter(iEne), flux_true_diff*scale_psf*scale_exp)
        self.dict_htg2_model[NHPSIDE] = htg2_model
        return htg2_model

    def write(self, path_file_out):
        file_out = ROOT.TFile(path_file_out, 'UPDATE')
        file_out.cd()
        self.HTG_SED.Write()
        for model in self.array_htg2_model.values():
            model.Write()


@click.command()
@click.argument('name', type=str)
@click.argument('ra', type=float)
@click.argument('dec', type=float)
@click.option('sed', type=str, default='')
@click.option('king', type=str, default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/Dispersion/AG_dispersion.root')
@click.option('acceptance', type=str, default='')
@click.argument('livetime', type=float)
def main(name, sed, ra, dec, king, acceptance, livetime):

    if sed!='':
        FILE_SED = ROOT.TFile(sed, 'READ')
        HTG_SED = FILE_SED.Get('hSED')
    else:
        HTG_SED = ''

    FILE_KING = ROOT.TFile(king, 'READ')
    TP_KING = (file_king.Get('htgKingN'), file_king.Get('htgKingS'), file_king.Get('htgKingG'))

    FILE_ACC = ROOT.TFile(acceptance, 'READ')
    HTG_ACC = FILR_ACC.Get('acc_cth_hist')

    FILE_LT = ROOT.TFile(livetime, 'READ')
    HTG_LT = FILR_LT.Get('acc_cth_hist')

    NSIDE_healpy = 512
    THRESHOLD_ROI = 15

    src_true = TruePointSource(name, HTG_SED, ra, dec)
    src_model = src_true.model(TP_KING, HTG_LT, HTG_ACC, NSIDE_healpy, THRESHOLD_ROI)

    src_true.write('ModelingPointSource_TEST.root')


if __name__ == '__main__':
    main()
