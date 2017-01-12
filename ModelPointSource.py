#!/usr/bin/env python

import sys
import os
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
from pCandleCatalogue import aObjectDict
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")


class TrueSource: 
    def __init__(self, name_src, path_file_out, spatial_type="Point-like", htg_sed=""):
        self.NAME = name_src
        self.PATH_OUT = path_file_out
        self.FILE_OUT= ROOT.TFile(self.PATH_OUT, 'UPDATE')
        self.FILE_OUT.cd()

        self.SPATIAL_TYPE = spatial_type
        #self.NBIN_ZENITH = 28
        #self.EDGE_ZENITH_LOW = 4.35
        self.ZENITH_CUT = 90.

        # Make dedicated directories
        TPL_DIR_DEDICATED = ['png']
        for dir_ded in TPL_DIR_DEDICATED:
            if os.path.isdir(dir_ded) is False:
                os.mkdir(dir_ded)

        if htg_sed=="":
            li_flux = aObjectDict[self.NAME]["flux"]
            li_flux_err = aObjectDict[self.NAME]["flux_err"]
            
            htg_sed = ROOT.TH1D('htg_sed_Mkn421', 'Mkn 421 SED for test', 14, 4.35, 5.75)
            for hbin in range(1, htg_sed.GetNbinsX()+1):
                htg_sed.SetBinContent(hbin, li_flux[hbin])
                htg_sed.SetBinError(hbin, li_flux_err[hbin])
            htg_sed.Rebin(2)
        self.TPL_STR_CLASS = ('CalOnlyR100', 'CalOnlyR30', 'CalOnlyR10')
        self.HTG_SED = htg_sed # [ photons / cm^2 / sec ]
        self.NBIN_ENERGY = htg_sed.GetXaxis().GetNbins()
        self.EDGE_ENERGY_LOW = htg_sed.GetXaxis().GetXmin()
        self.EDGE_ENERGY_UP = htg_sed.GetXaxis().GetXmax()
        # self.dctt_htg2_model = {} # Key is NSIDE, CLASS
        # self.dctt_htg2_on = {} # Key is NSIDE, CLASS
        # self.dctt_map = {} # Key is NSIDE, CLASS
        # self.dctt_map_energy = {} # Key is NSIDE, CLASS
        # self.dct_htg_sed_model = {} # Key is CLASS
        # self.dct_htg_sederrSq_model = {} # Key is CLASS
        # self.dctt_htg_sed_model_on = {} # Key is NSIDE, CLASS

        # for cla in self.TPL_STR_CLASS:

        #     self.dctt_htg2_model[cla] = {}
        #     self.dctt_htg2_on[cla] = {}
        #     self.dctt_map[cla] = {}
            #self.dctt_map_energy[cla] = {}

            #self.dct_htg_sed_model[cla] = self.HTG_SED.Clone('{0}_modeled_{1}'.format(htg_sed.GetName(), cla))
            #self.dct_htg_sed_model[cla] = self.dct_htg_sed_model[cla].Clone('{0}_errSq'.format(self.dct_htg_sed_model[cla].GetName()))
            #self.dctt_htg_sed_model_on[cla] = {} # Key is NSIDE, CLASS

            #for hbin in range(0, self.dct_htg_sed_model[cla].GetNbinsX()+2):
                #self.dct_htg_sed_model.SetBinContent(hbin, 0)
                #self.dct_htg_sed_model.SetBinError(hbin, 0)
                #self.dct_htg_sederrSq_model.SetBinContent(hbin, 0)
                #self.dct_htg_sederrSq_model.SetBinError(hbin, 0)

        #for clas in self.TPL_STR_CLASS:
         #   self.dct_htg_sed_model_on[clas] = self.dct_htg_sed_model[cla].Clone('{0}_on_{1}'.format(self.dct_htg_sed_model[cla].GetName(), clas))
            #self.htg_sederrSq_model_on = self.htg_sederrSq_model.Clone('{0}_on_{1}'.format(self.htg_sederrSq_model.GetName(), clas))


class TruePointSource(TrueSource):
    def __init__(self, name_src, path_file_out, htg_sed, ra, dec):
        TrueSource.__init__(self, name_src, path_file_out, "Point-like", htg_sed)
        self.RA = ra
        self.DEC = dec
        self.RA_RAD = radians(self.RA)
        self.DEC_RAD = radians(self.DEC)
        #self.dctt_arr_map = {}        
        #self.dctt_arr_map_energy = {}
        self.tp_rotate = (self.RA, self.DEC, 0)

    def model(self, TP_HTG_KING, HTG2_LIVETIME, NHPSIDE=256, THRESHOLD_ANGDIST=15, tpl_path_perf=('/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R30_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R10_perf.root')):
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
        for ipix in range(NPIX):
            tp_pix = hppf.pix2ang(NHPSIDE, ipix)
            dict_angdist[ipix] = hp.rotator.angdist(tp_pix, TP_DIR_TRUE) #float(ang_pix2src.to_string(unit=u.rad, decimal=True)) # Angular distance between the source and each pixel in radians
            if dict_angdist[ipix] < THRESHOLD_ANGDIST_RADIANS:
                set_pix_roi.add(ipix)

        print 'ROI pixels:', set_pix_roi

        HTG1_LT = HTG2_LIVETIME.ProjectionX("{0}_projTheta".format(HTG2_LIVETIME.GetName()), 1, HTG2_LIVETIME.GetYaxis().FindBin(self.ZENITH_CUT)-1)
#        htg2_model = ROOT.TH2D("htg2_model_{0}".format(NHPSIDE), "Model of {0} (NSIDE={1})".format(self.NAME, NHPSIDE), NPIX, 0, NPIX, self.NBIN_ENERGY, self.EDGE_ENERGY_LOW, self.EDGE_ENERGY_UP)

        # Preparation of ON region setup
        #dct_htg2_model = {}
        #dct_htg2_model_on = {}
        dct_path_perf = {}
        dct_file_perf = {}
        dct_htg_psf95 = {}
        dct_htg_acc = {}
        #dct_htg_sed_model_on = {}
#        if self.HTG_SED.GetNbinsX()!=self.NBIN_ENERGY:
 #           print "Number of energy bins is not matched between {0} and {1}.".format(self.HTG_SED.GetName(), dct_htg2_model[cla].GetName())
        # PSF
        fc_King_annulus = ROOT.TF1("fc_King_annulus", "TMath::Sin(x)*[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/[1]**2", 0, pi)
        fc_King = ROOT.TF1("fc_King", "[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/2./TMath::Pi()/[1]**2", 0, pi)

        for (icla,cla) in enumerate(self.TPL_STR_CLASS):
            print '======================='
            print cla
            print '======================='

            htg2_model = ROOT.TH2D("htg2_model_{0}_NSIDE{1}".format(cla, NHPSIDE), "Model of {0} {1} (NSIDE={2})".format(cla, self.NAME, NHPSIDE), NPIX, 0, NPIX, self.NBIN_ENERGY, self.EDGE_ENERGY_LOW, self.EDGE_ENERGY_UP)
            htg2_model_on = htg2_model.Clone('{0}_ON'.format(htg2_model.GetName()))
            htg_sed_model = ROOT.TH1D('htg_sed_model_{0}'.format(cla), 'SED model of {0} {1}'.format(cla, self.NAME), self.HTG_SED.GetNbinsX(), self.HTG_SED.GetXaxis().GetXmin(), self.HTG_SED.GetXaxis().GetXmax())
            htg_sederrSq_model = htg_sed_model.Clone('{0}_errSq'.format(htg_sed_model.GetName()))
            #htg_sed_model_on = htg_sed_model.Clone('{0}_ON_{1}'.format(htg_sed_model.GetName(), cla))
            dct_path_perf[cla] = tpl_path_perf[icla]
            dct_file_perf[cla] = ROOT.TFile(dct_path_perf[cla], 'READ')
            dct_htg_psf95[cla] = dct_file_perf[cla].Get('psf_q95_hist')
            dct_htg_acc[cla] = dct_file_perf[cla].Get('acc_cth_hist')

            for iEne in range(1, self.HTG_SED.GetNbinsX()+1):
                print 'Energy:', self.HTG_SED.GetXaxis().GetBinLowEdge(iEne), '-', self.HTG_SED.GetXaxis().GetBinUpEdge(iEne)
                flux_true_itgl = self.HTG_SED.GetBinContent(iEne)
                fluxerr_true_itgl = self.HTG_SED.GetBinError(iEne)
                print '  Flux:', flux_true_itgl, '+/-', fluxerr_true_itgl, '[photons cm^-2 s^-1]'
                for iTh in range(1, HTG1_LT.GetNbinsX()+1):
                    print '  cos(Inclination angle):', HTG1_LT.GetXaxis().GetBinLowEdge(iTh), '-', HTG1_LT.GetXaxis().GetBinUpEdge(iTh)
                    scale_aeff = dct_htg_acc[cla].GetBinContent(dct_htg_acc[cla].GetXaxis().FindBin(self.HTG_SED.GetXaxis().GetBinCenter(iEne)), dct_htg_acc[cla].GetYaxis().FindBin(HTG1_LT.GetXaxis().GetBinCenter(iTh))) / (2*math.pi*dct_htg_acc[cla].GetYaxis().GetBinWidth(dct_htg_acc[cla].GetYaxis().FindBin(HTG1_LT.GetXaxis().GetBinCenter(iTh)))) # Effective area [m^2]
                    print '    Effective area:', scale_aeff, 'm^2'
                    scale_exp = scale_aeff * 100.**2 * HTG1_LT.GetBinContent(iTh) # Integrated exposure value for a certain energy and inclination angle [cm^2 s]
                    print '    Exposure:', scale_exp, 'cm^2 s'
                    nphoton = flux_true_itgl*scale_exp
                    nphotonerr = fluxerr_true_itgl*scale_exp
                    htg_sed_model.Fill(htg_sed_model.GetBinCenter(iEne), nphoton)
                #htg2_model_errSq.Fill(iEne, nphotonerr**2)
                    htg_sederrSq_model.Fill(htg_sederrSq_model.GetBinCenter(iEne), nphotonerr**2)
                    print '    Number of photons:', nphoton, '+/-', nphotonerr
                    kxbin = TP_HTG_KING[0].GetXaxis().FindBin(self.HTG_SED.GetXaxis().GetBinCenter(iEne))
                    kybin = TP_HTG_KING[0].GetYaxis().FindBin(HTG1_LT.GetBinCenter(iTh))
                    if nphoton>0 and kxbin>0 and kybin>0:
                        for ipar in range(3): # Setting the parameters of King function
                        # PSF
                            par_value = TP_HTG_KING[ipar].GetBinContent(kxbin, kybin)
                            print '    Parameter No.{0}:'.format(ipar), par_value
                            fc_King_annulus.FixParameter(ipar, par_value)
                            fc_King.FixParameter(ipar, par_value)
                        factor_norm = 1.0/fc_King_annulus.Integral(0, pi)
                        print '    Normalization factor:', factor_norm
                        for ipix in set_pix_roi:
                            scale_psf = fc_King.Eval(dict_angdist[ipix])
                            htg2_model.Fill(ipix+0.5, htg2_model.GetYaxis().GetBinCenter(iEne), nphoton*scale_psf*factor_norm*sa_pix)
#                        for cla in self.TPL_STR_CLASS:
                            deg_psf95 = dct_htg_psf95[cla].GetBinContent(dct_htg_psf95[cla].FindBin(self.HTG_SED.GetBinCenter(iEne)))
                            if radians(min(15, deg_psf95))>dict_angdist[ipix]:
                                htg2_model_on.SetBinContent(ipix, iEne, 1)
                            else:
                                htg2_model_on.SetBinContent(ipix, iEne, 0)
                #htg_sed_model.SetBinError(iEne, math.sqrt(htg_sederrSq_model.GetBinContent(iEne)))
                htg_sed_model.SetBinError(iEne, self.HTG_SED.GetBinError(iEne)/self.HTG_SED.GetBinContent(iEne)*htg_sed_model.GetBinContent(iEne))
                print '  Observable photon number:', htg_sed_model.GetBinContent(iEne), '+/-', htg_sed_model.GetBinError(iEne), 'photons'
                print ''

            print htg2_model_on.Integral()
            htg2_model_on.Multiply(htg2_model)
            for ienr in range(1, htg2_model_on.GetYaxis().GetNbins()+1):
                for ipix in range(1, htg2_model_on.GetXaxis().GetNbins()+1):
                    content = htg2_model.GetBinContent(ipix, ienr)
                    content_err = htg2_model.GetBinError(ipix, ienr)
                    content_on = htg2_model_on.GetBinContent(ipix, ienr)
                    if content>0:
                        htg2_model_on.SetBinError(ipix, ienr, content_err*content_on/content)
                    else:
                        htg2_model_on.SetBinError(ipix, ienr, 0)
            print htg2_model_on.Integral()
            htg_sed_model_on = htg2_model_on.ProjectionY('{0}_ON'.format(htg_sed_model.GetName()))
            print 'ON photon number:', htg_sed_model_on.Integral()

            for ienr in range(1, htg_sed_model_on.GetXaxis().GetNbins()+1):
                content = htg_sed_model.GetBinContent(ienr)
                content_err = htg_sed_model.GetBinError(ienr)
                content_on = htg_sed_model_on.GetBinContent(ienr)
                if content>0:
                    htg_sed_model_on.SetBinError(ienr, content_err*content_on/content)
                else:
                    htg_sed_model_on.SetBinError(ienr, 0)

            print 'Making map...'
            arr_map = []
            #nparr_map = np.zeros((self.HTG_SED.GetNbinsX()+1, htg2_model.GetXaxis().GetNbins())) # Energy vs. pixel
            arr_map_energy = []
            for iEne in range(0, self.HTG_SED.GetNbinsX()+1):
                arr_map.append([])
                if iEne>0:
                    arr_map_energy.append((self.HTG_SED.GetXaxis().GetBinLowEdge(iEne), self.HTG_SED.GetXaxis().GetBinUpEdge(iEne)))
                    for ipix in range(htg2_model.GetXaxis().GetNbins()):
                        if ipix in set_pix_roi:
                            arr_map[-1].append(htg2_model.GetBinContent(ipix+1, iEne))
                        else:
                            arr_map[-1].append(hppf.UNSEEN)
                else:
                    arr_map_energy.append((self.HTG_SED.GetXaxis().GetBinLowEdge(1), self.HTG_SED.GetXaxis().GetBinUpEdge(self.HTG_SED.GetNbinsX())))
                    for ipix in range(htg2_model.GetXaxis().GetNbins()):
                        if ipix in set_pix_roi:
                            arr_map[-1].append(htg2_model.Integral(ipix+1, ipix+1, 1, htg2_model.GetNbinsY()))
                        else:
                            arr_map[-1].append(hppf.UNSEEN)

                print '  Energy:', arr_map_energy[-1][0], '-', arr_map_energy[-1][1]
                nparr_map = np.array(arr_map[-1]) # Take this energy bin
                hp.visufunc.cartview(nparr_map, iEne, self.tp_rotate, unit='cm^-2 s^-1 / {0:.1E} sr'.format(sa_pix), lonra=[-THRESHOLD_ANGDIST, THRESHOLD_ANGDIST], latra=[-THRESHOLD_ANGDIST, THRESHOLD_ANGDIST], title='{0} ({1:.3f} - {2:.3f} GeV)'.format(self.NAME, pow(10, arr_map_energy[-1][0]-3), pow(10, arr_map_energy[-1][1]-3)), min=0, flip='astro')
                plt.savefig("png/Model-{0}_{1}_NSIDE{2}_{3}-{4}.png".format(self.NAME, cla, NHPSIDE, int(100*arr_map_energy[-1][0]+0.5), int(100*arr_map_energy[-1][1]+0.5)))
                plt.close()
            #hp.fitsfunc.write_map("{0}_NSIDE{1}_{2}-{3}.fits".format(self.NAME, NHPSIDE, int(100*arr_map_energy[-1][0]+0.5), int(100*arr_map_energy[-1][1]+0.5)), nparr_map)
                #htg1_model_px = htg2_model.ProjectionY('{0}_px{1}'.format(htg2_model.GetName(), iEne), 1, htg2_model.GetXaxis().GetNbins())
                #print '  ', htg1_model_px.GetBinContent(iEne), 'photons'
            
            self.FILE_OUT.cd()
            htg_sed_model.Write()
            htg2_model.Write()
            htg_sed_model_on.Write()
        
            #self.dctt_htg2_model[NHPSIDE][cla] = htg2_model
            #self.dctt_arr_map[NHPSIDE][cla] = arr_map
            #self.dctt_arr_map_energy[NHPSIDE][cla] = arr_map_energy
        self.HTG_SED.Write()        


@click.command()
@click.argument('name', type=str)
@click.argument('ra', type=float)
@click.argument('dec', type=float)
@click.option('--sed', type=str, default='')
@click.option('--king', type=str, default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/Dispersion/AG_dispersion.root')
#@click.option('--acceptance', type=str, default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R100_perf.root')
@click.argument('livetime', type=str)
@click.option('--suffix', type=str, default='')
@click.option('--nside', '-n', type=int, default=256)
def main(name, sed, ra, dec, king, livetime, suffix, nside):
    if math.log(nside,2)!=int(math.log(nside,2)) or nside>2**30:
        raise click.BadParameter('nside must be a power of 2, less than 2**30!!!')

    if sed!='':
        FILE_SED = ROOT.TFile(sed, 'READ')
        HTG_SED = FILE_SED.Get('hSED')
        print HTG_SED.GetName(), 'has been found.'
    else:
        HTG_SED= ''
    FILE_KING = ROOT.TFile(king, 'READ')
    TP_KING = (FILE_KING.Get('htgKingN'), FILE_KING.Get('htgKingS'), FILE_KING.Get('htgKingG'))
    FILE_LT = ROOT.TFile(livetime, 'READ')
    HTG_LT = FILE_LT.Get('htgLt_0_yx') #TBD
    print HTG_LT.GetName(), 'has been found.'

    THRESHOLD_ROI = 20

    if suffix!='':
        suffix = "_" + suffix

    src_true = TruePointSource(name, 'ModelingPointSource_{0}{1}.root'.format(name, suffix), HTG_SED, ra, dec)
    src_model = src_true.model(TP_KING, HTG_LT, nside, THRESHOLD_ROI)


if __name__ == '__main__':
    main()
