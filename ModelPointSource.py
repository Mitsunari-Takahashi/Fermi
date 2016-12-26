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
import click
import matplotlib.pyplot as plt
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")


class TrueSource: 
    def __init__(self, name_src, spatial_type="Point-like", htg_sed=""):
        self.NAME = name_src
        self.SPATIAL_TYPE = spatial_type
        #self.NBIN_ZENITH = 28
        #self.EDGE_ZENITH_LOW = 4.35
        self.ZENITH_CUT = 90.
        if htg_sed=="":
            li_flux = [  7.38720813e-10,   5.42467533e-10,   4.28598094e-10,   4.48707388e-10,
                         3.34952271e-10,   2.71831631e-10,   1.92920102e-10,   1.21634262e-10,
                         1.06100454e-10,   1.15034194e-10,   6.71599557e-11,   4.14993127e-11,
                         4.11690073e-11,   3.57418942e-11,   1.52070785e-11,   1.27542759e-11]
            li_flux_err = [  4.69432924e-11,   4.02066133e-11,   3.52432951e-11,   3.58411223e-11,
                             3.10366595e-11,   2.79193559e-11,   2.35158326e-11,   1.86931808e-11,
                             1.75856652e-11,   1.83589107e-11,   1.39442576e-11,   1.10684638e-11,
                             1.09408549e-11,   1.02402336e-11,   6.82109584e-12,   6.39910338e-12]
            htg_sed = ROOT.TH1D('htg_sed_Mkn421', 'Mkn 421 SED for test', 14, 4.35, 5.75)
            for hbin in range(1, htg_sed.GetNbinsX()+1):
                htg_sed.SetBinContent(hbin, li_flux[hbin])
                htg_sed.SetBinError(hbin, li_flux_err[hbin])
            htg_sed.Rebin(2)
        self.HTG_SED = htg_sed # [ photons / cm^2 / sec ]
        self.NBIN_ENERGY = htg_sed.GetXaxis().GetNbins()
        self.EDGE_ENERGY_LOW = htg_sed.GetXaxis().GetXmin()
        self.EDGE_ENERGY_UP = htg_sed.GetXaxis().GetXmax()
        self.dict_htg2_model = {} # Key is NSIDE
        self.dict_map = {} # Key is NSIDE
        self.dict_map_energy = {} # Key is NSIDE


class TruePointSource(TrueSource):
    def __init__(self, name_src, htg_sed, ra, dec):
        TrueSource.__init__(self, name_src, "Point-like", htg_sed)
        self.RA = ra
        self.DEC = dec
        self.RA_RAD = radians(self.RA)
        self.DEC_RAD = radians(self.DEC)
        self.dict_arr_map = {}        
        self.dict_arr_map_energy = {}


    def model(self, TP_HTG_KING, HTG2_LIVETIME, HTG2_ACCEPTANCE, NHPSIDE=256, THRESHOLD_ANGDIST=15):
        """Model the point source with the PSF of King function. Energy dispersion is ignored currently.
        """
        PIX_TRUE = hppf.ang2pix(NHPSIDE, pi/2.-self.DEC_RAD, self.RA_RAD) # #Pixel the true position of the source locates in
        TP_DIR_TRUE = (pi/2.-self.DEC_RAD, self.RA_RAD)
        print 'NSIDE:', NHPSIDE
        NPIX = hppf.nside2npix(NHPSIDE)
        print NPIX, 'pixels.'
        print 'Pixel resolution:', degrees(hppf.nside2resol(NHPSIDE)), 'deg'
        sa_pix = hppf.nside2pixarea(NHPSIDE) # Solid angle of a pixel [sr]
        print 'Solid angle of each pixel:', sa_pix, 'sr'
        coo_src = SkyCoord(self.RA, self.DEC, unit="deg") # True coordinate of the source
        THRESHOLD_ANGDIST_RADIANS = radians(THRESHOLD_ANGDIST)

        # ROI
        set_pix_roi = set([])
        dict_angdist = {}
#        li_lon = [self.RA_RAD-pi, self.RA_RAD-pi] #Min, Max
#        li_lat = [self.DEC_RAD, self.DEC_RAD] #Min, Max
        for ipix in range(NPIX):
            tp_pix = hppf.pix2ang(NHPSIDE, ipix)
            #coo_pix = SkyCoord(tp_pix[1], pi/2.-tp_pix[0], unit="rad")
            #ang_pix2src = coo_src.separation(coo_pix) # Angular distance between the souece and the pixel center
            dict_angdist[ipix] = hp.rotator.angdist(tp_pix, TP_DIR_TRUE) #float(ang_pix2src.to_string(unit=u.rad, decimal=True)) # Angular distance between the source and each pixel in radians
            if dict_angdist[ipix] < THRESHOLD_ANGDIST_RADIANS:
                set_pix_roi.add(ipix)
                # arr_ang_pix = hppf.pix2ang(NHPSIDE, ipix) #Theta, Phi
                # lat_pix = pi/2.-arr_ang_pix[0] #Theta -> Lat
                # lon_pix = arr_ang_pix[1]-pi #Theta -> Lat
                # if lat_pix<li_lat[0]:
                #     li_lat[0] = lat_pix
                # if lat_pix>li_lat[1]:
                #     li_lat[1] = lat_pix
                # if lon_pix<li_lon[0]:
                #     li_lon[0] = lon_pix
                # if lon_pix>li_lon[1]:
                #     li_lon[1] = lon_pix

        print 'ROI pixels:', set_pix_roi
        # li_lon_deg = (degrees(li_lon[0]), degrees(li_lon[1]))
        # li_lat_deg = (degrees(li_lat[0]), degrees(li_lat[1]))
        # print 'Longitude range:', li_lon_deg, '[deg]'
        # print 'Latitude range:', li_lat_deg, '[deg]'

        HTG1_LT = HTG2_LIVETIME.ProjectionX("{0}_projTheta".format(HTG2_LIVETIME.GetName()), 1, HTG2_LIVETIME.GetYaxis().FindBin(self.ZENITH_CUT)-1)
        htg2_model = ROOT.TH2D("htg2_model_{0}".format(NHPSIDE), "Model of {0} (NSIDE={1})".format(self.NAME, NHPSIDE), NPIX, 0, NPIX, self.NBIN_ENERGY, self.EDGE_ENERGY_LOW, self.EDGE_ENERGY_UP)
        if self.HTG_SED.GetNbinsX()!=self.NBIN_ENERGY:
            print "Number of energy bins is not matched between {0} and {1}.".format(self.HTG_SED.GetName(), htg2_model.GetName())
        # PSF
        fc_King_annulus = ROOT.TF1("fc_King_annulus", "TMath::Sin(x)*[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/[1]**2", 0, pi)
        fc_King = ROOT.TF1("fc_King", "[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/2./TMath::Pi()/[1]**2", 0, pi)

        for iEne in range(1, self.HTG_SED.GetNbinsX()+1):
            print 'Energy:', self.HTG_SED.GetXaxis().GetBinLowEdge(iEne), '-', self.HTG_SED.GetXaxis().GetBinUpEdge(iEne)
            flux_true_diff = self.HTG_SED.GetBinContent(iEne)
            print '  Flux:', flux_true_diff
            for iTh in range(1, HTG1_LT.GetNbinsX()+1):
                print '  cos(Inclination angle):', HTG1_LT.GetXaxis().GetBinLowEdge(iTh), '-', HTG1_LT.GetXaxis().GetBinUpEdge(iTh)
                scale_aeff = HTG2_ACCEPTANCE.GetBinContent(HTG2_ACCEPTANCE.GetXaxis().FindBin(self.HTG_SED.GetXaxis().GetBinCenter(iEne)), HTG2_ACCEPTANCE.GetYaxis().FindBin(HTG1_LT.GetXaxis().GetBinCenter(iTh))) / (2*math.pi*HTG2_ACCEPTANCE.GetYaxis().GetBinWidth(HTG2_ACCEPTANCE.GetYaxis().FindBin(HTG1_LT.GetXaxis().GetBinCenter(iTh)))) # Effective area [m^2]
                print '    Effective area:', scale_aeff, 'm^2'
                scale_exp = scale_aeff * 100.**2 * HTG1_LT.GetBinContent(iTh) # Integrated exposure value for a certain energy and inclination angle [cm^2 s]
                print '    Exposure:', scale_exp, 'cm^2 s'
                if scale_exp>0:
                    for ipar in range(3): # Setting the parameters of King function
                        # PSF
                        par_value = TP_HTG_KING[ipar].GetBinContent(TP_HTG_KING[ipar].GetXaxis().FindBin(self.HTG_SED.GetXaxis().GetBinCenter(iEne)), TP_HTG_KING[ipar].GetYaxis().FindBin(HTG1_LT.GetBinCenter(iTh)))
                        print '    Parameter No.{0}:'.format(ipar), par_value
                        fc_King_annulus.FixParameter(ipar, par_value)
                        fc_King.FixParameter(ipar, par_value)
                    factor_norm = 1.0/fc_King_annulus.Integral(0, pi)
                    print '    Normalization factor:', factor_norm
                    for ipix in set_pix_roi:
                        scale_psf = fc_King.Eval(dict_angdist[ipix])*hppf.nside2pixarea(NHPSIDE)
                        htg2_model.Fill(ipix+0.5, htg2_model.GetYaxis().GetBinCenter(iEne), flux_true_diff*scale_psf*scale_exp*factor_norm*sa_pix)
            print ''

        print 'Making map...'
        arr_map = []
        arr_map_energy = []
        tp_rotate = (self.RA, self.DEC, 0)
        for iEne in range(1, self.HTG_SED.GetNbinsX()+1):
            arr_map_energy.append((self.HTG_SED.GetXaxis().GetBinLowEdge(iEne), self.HTG_SED.GetXaxis().GetBinUpEdge(iEne)))
            print '  Energy:', arr_map_energy[-1][0], '-', arr_map_energy[-1][1]
            arr_map.append([])
            for ipix in range(htg2_model.GetXaxis().GetNbins()):
                if ipix in set_pix_roi:
                    arr_map[-1].append(htg2_model.GetBinContent(ipix+1, iEne))
                else:
                    arr_map[-1].append(hppf.UNSEEN)
            #print arr_map[-1]
            nparr_map = np.array(arr_map[-1])
            hp.visufunc.cartview(nparr_map, iEne, tp_rotate, unit='photons cm^-2 s^-1', lonra=[-THRESHOLD_ANGDIST, THRESHOLD_ANGDIST], latra=[-THRESHOLD_ANGDIST, THRESHOLD_ANGDIST], title='{0} ({1:.3f} - {2:.3f} GeV)'.format(self.NAME, pow(10, arr_map_energy[-1][0]-3), pow(10, arr_map_energy[-1][1]-3)), min=0)#, flip='astro')
            #plt.show()
            plt.savefig("{0}_NSIDE{1}_{2}-{3}.png".format(self.NAME, NHPSIDE, int(100*arr_map_energy[-1][0]+0.5), int(100*arr_map_energy[-1][1]+0.5)))
            hp.fitsfunc.write_map("{0}_NSIDE{1}_{2}-{3}.fits".format(self.NAME, NHPSIDE, int(100*arr_map_energy[-1][0]+0.5), int(100*arr_map_energy[-1][1]+0.5)), nparr_map)
            
        self.dict_htg2_model[NHPSIDE] = htg2_model
        self.dict_arr_map[NHPSIDE] = arr_map
        self.dict_arr_map_energy[NHPSIDE] = arr_map_energy
        return htg2_model

    def write(self, path_file_out):
        file_out = ROOT.TFile(path_file_out, 'UPDATE')
        file_out.cd()
        self.HTG_SED.Write()
        for model in self.dict_htg2_model.values():
            model.Write()


@click.command()
@click.argument('name', type=str)
@click.argument('ra', type=float)
@click.argument('dec', type=float)
@click.option('--sed', type=str, default='')
@click.option('--king', type=str, default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/Dispersion/AG_dispersion.root')
@click.option('--acceptance', type=str, default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R100_perf.root')
@click.argument('livetime', type=str)
@click.option('--suffix', type=str, default='')
@click.option('--nside', '-n', type=int, default='')
def main(name, sed, ra, dec, king, acceptance, livetime, suffix, nside):
    if math.log(nside,2)!=int(math.log(nside,2)) or nside>2**30:
        raise click.BadParameter('nside must be a power of 2, less than 2**30!!!')

    if sed!='':
        FILE_SED = ROOT.TFile(sed, 'READ')
        HTG_SED = FILE_SED.Get('hSED')
        print HTG_SED.GetName(), 'has been found.'
    else:
        HTG_SED= ''
    #     HTG_SED = ROOT.TH1D("hSED", "SED created automatically", 28, 4.35, 5.75)
    #     for ibin in range(1, 1+HTG_SED.GetNbinsX()):
    #         HTG_SED.SetBinContent(ibin, 1)
    #     print HTG_SED.GetName(), 'has been created.'

    FILE_KING = ROOT.TFile(king, 'READ')
    TP_KING = (FILE_KING.Get('htgKingN'), FILE_KING.Get('htgKingS'), FILE_KING.Get('htgKingG'))

    FILE_ACC = ROOT.TFile(acceptance, 'READ')
    HTG_ACC = FILE_ACC.Get('acc_cth_hist')
    print HTG_ACC.GetName(), 'has been found.'

    FILE_LT = ROOT.TFile(livetime, 'READ')
    HTG_LT = FILE_LT.Get('htgLt_GRB160509374_0_yx') #TBD
    print HTG_LT.GetName(), 'has been found.'

#    NSIDE_healpy = 128
    THRESHOLD_ROI = 15

    src_true = TruePointSource(name, HTG_SED, ra, dec)
    src_model = src_true.model(TP_KING, HTG_LT, HTG_ACC, nside, THRESHOLD_ROI)

    if suffix=='':
        suffix = "_" + suffix
    src_true.write('ModelingPointSource{0}.root'.format(suffix))


if __name__ == '__main__':
    main()
