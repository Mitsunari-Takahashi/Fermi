import ROOT
from astropy.io import fits
import numpy as np
from pAnalysisConfig import *
import pColor
import pCommon
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp
from healpy import pixelfunc as hppf
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TH1, TH2


class Healcube:
    def __init__(self, str_gamsrc, str_evtclass, htg_cube, lon_cntr=0, lat_cntr=0, char_coord='E', deg_radius=None, lst_validpix=None):

        self.htg = htg_cube
        self.eregion=EnergyLogRegion(self.htg.GetXaxis().GetNbins(), self.htg.GetXaxis().GetBinLowEdge(1), self.htg.GetXaxis().GetBinWidth(1))
        print 'Energy region:', self.eregion.printR()
        self.cthregion=EnergyLogRegion(self.htg.GetYaxis().GetNbins(), self.htg.GetYaxis().GetBinLowEdge(1), self.htg.GetYaxis().GetBinWidth(1))
        print 'cos(theta) region:', self.cthregion.printR()
        self.NPIX=self.htg.GetZaxis().GetNbins()
        self.NSIDE=hppf.npix2nside(self.NPIX)

        if lst_validpix is None:
            self.validpix = range(self.NPIX)
        else:
            self.validpix = lst_validpix

        self.name=str_gamsrc
        self.evtclass = str_evtclass
        self.loncntr=lon_cntr
        self.latcntr=lat_cntr
        self.coord=char_coord
        self.radius=deg_radius 

    def setmap(self):
        self.map = []
        for ienr in range(self.eregion.nBin+1):
            self.map.append([])
            enrlow = ienr if ienr>0 else 1
            enrup = ienr if ienr>0 else self.htg.GetXaxis().GetNbins()
            for icth in range(self.cthregion.nBin+1):
                self.map[-1].append(np.zeros(self.NPIX))
                cthlow = icth if icth>0 else 1
                cthup = icth if icth>0 else self.htg.GetYaxis().GetNbins()
                htg_pZ = self.htg.ProjectionZ('{0}_pz{1}{2}'.format(self.htg.GetName(), ienr, icth), enrlow, enrup, cthlow, cthup)
                for ipix in range(self.NPIX):
                    if ipix in self.validpix:
                        self.map[-1][-1][ipix] = htg_pZ.GetBinContent(ipix+1)
                    else:
                        self.map[-1][-1][ipix] = hppf.UNSEEN

    def draw(self, mollweide=False):
        for ienr in range(self.eregion.nBin+1):
            enrlow = ienr-1 if ienr>0 else 0
            enrup = ienr-1 if ienr>0 else self.eregion.nBin-1
            if mollweide==False:
                hp.visufunc.cartview(self.map[ienr][0], rot=(self.loncntr, self.latcntr, 0), coord=self.coord, lonra=[-self.radius, self.radius], latra=[-self.radius, self.radius], min=0, flip='astro', title="{0} {1} ({2:.1f} - {3:.1f} GeV)".format(self.evtclass, self.name, 10**(self.eregion.getBin(enrlow)[0]-3), 10**(self.eregion.getBin(enrup)[1]-3)), unit='counts')
            else:
                hp.visufunc.mollview(self.map[ienr][0], rot=(self.loncntr, self.latcntr, 0), coord=self.coord, min=0, flip='astro', title="{0} {1} ({2:.1f} - {3:.1f} GeV)".format(self.evtclass, self.name, 10**(self.eregion.getBin(enrlow)[0]-3), 10**(self.eregion.getBin(enrup)[1]-3)), unit='counts')
            plt.savefig("{0}_{1}_E{2}-{3}.png".format(self.name, self.evtclass, int(100*self.eregion.getBin(enrlow)[0]+0.5), int(100*self.eregion.getBin(enrup)[1]+0.5)))
            plt.clf()


    def smear(self, nparr_dist, path_king, deg_threshold=None):

        sa_pix = hppf.nside2pixarea(self.NSIDE) # Solid angle of a pixel [sr]
        htg_smr = self.htg.Clone("{0}_smeared".format(self.GetName()))
        for hx in range(htg_smr.GetXaxis().GetNbins()+2):
            for hy in range(htg_smr.GetYaxis().GetNbins()+2):
                for hz in range(htg_smr.GetZaxis().GetNbins()+2):
                    htg_smr.SetBinContent(hx, hy, hz, 0)
                    htg_smr.SetBinError(hx, hy, hz, 0)
        FILE_KING = ROOT.TFile(path_king, 'READ')
        TP_HTG_KING = (FILE_KING.Get('htgKingN'), FILE_KING.Get('htgKingS'), FILE_KING.Get('htgKingG'))
        fc_King_annulus = ROOT.TF1("fc_King_annulus", "TMath::Sin(x)*[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/[1]**2", 0, pi)
        fc_King = ROOT.TF1("fc_King", "[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/2./TMath::Pi()/[1]**2", 0, pi)

        for ienr in range(1, self.eregion.nBin+1):
            kxbin = TP_HTG_KING[0].GetXaxis().FindBin(self.eregion.getBinCenter(ienr-1))
            for icth in range(1, self.cthregion.nBin+1):
                kybin = TP_HTG_KING[0].GetYaxis().FindBin(self.cthregion.getBinCenter(icth-1))
                if kxbin>0 and kybin>0:
                    for ipar in range(3): # Setting the parameters of King function
                        # PSF
                        par_value = TP_HTG_KING[ipar].GetBinContent(kxbin, kybin)
                        #print '    Parameter No.{0}:'.format(ipar), par_value
                        fc_King_annulus.FixParameter(ipar, par_value)
                        fc_King.FixParameter(ipar, par_value)
                    factor_norm = 1.0/fc_King_annulus.Integral(0, pi)
                    for ipix in range(1, self.NPIX+1):
                        cnt = self.htg.GetBinContent(ienr, icth, ipix)
                        if cnt>0:
                            for jpix in range(1, self.NPIX+1):
                                angdist = nparr_dist[ipix-1][jpix-1]
                                htg_smr.Fill(ienr, icth, jpix, cnt*fc_King.Eval(angdist)*factor_norm*sa_pix)
        return htg_smr

        

def Setdistance(nside):
    """Returns a table of distance between each HEALPix pixel in radians.
"""
    nparr_dist = np.zeros((nside, nside))
    for (ipix, npar_pix) in enumerate(nparr_dist[0]):
        for (jpix, rad_pix) in enumerate(npar_pix):
            rad_pix = hp.rotator.angdist(hppf.pix2ang(nside, ipix), hppf.pix2ang(nside, jpix))
            nparr_dist[ipix][jpix] = rad_pix
    np.save('PixelsDistance_{0}.npy'.format(nside), nparr_dist)
    return nparr_dist # In radians
        
