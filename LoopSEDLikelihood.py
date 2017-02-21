#!/usr/bin/env python

import sys
#import os
import numpy as np
#import yaml
#import datetime
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TH1, TH2
ROOT.gROOT.SetBatch()
from pColor import *
#from pMETandMJD import *
import ModelPointSource
import ExtrapolateFlux


def LoopSEDLikelihood(name, observed, tpl_dnde_log, tpl_index, eref, ra, dec, king, livetime, suffix, nside, redshift, addregular, dnde_pri_min, dnde_pri_max, idx_pri_min, idx_pri_max, binw=0.294, enr_min=4.00, nbin=4, t_start=770., t_stop=8233.):
    #t_stop will be included in the time window.

    TPL_CLASS = ('CalOnlyR100',) #, 'CalOnlyR30', 'CalOnlyR10')
    NREBIN = 1
    #SCALE_FLUX = 1.0e-13

    FILE_IN = ROOT.TFile(observed, 'READ')
    HTG_OBS = FILE_IN.Get('spectrum_observed')
    print HTG_OBS, 'has been found.'
    HTG_OBS_ENR = HTG_OBS.ProjectionY('{0}_projEnr'.format(HTG_OBS.GetName()), HTG_OBS.GetXaxis().FindBin(t_start), HTG_OBS.GetXaxis().FindBin(t_stop))
    print HTG_OBS_ENR.Integral(HTG_OBS_ENR.GetXaxis().FindBin(enr_min), HTG_OBS_ENR.GetXaxis().FindBin(enr_min+binw*nbin)), 'observed events.'
    HTG_OBS_ENR.Rebin(NREBIN)
    HTG_OBS_ENR.SetLineWidth(0)
    HTG_OBS_ENR.SetLineColor(ROOT.kRed)
    HTG_OBS_ENR.SetMarkerColor(ROOT.kRed)
    HTG_OBS_ENR.SetMarkerStyle(20)
    HTG_OBS_ENR.SetFillStyle(0)

    PATH_FILE_OUT = 'LoopLikelihood_{0}{1}'.format(name, suffix)
    FILE_OUT = ROOT.TFile('{0}.root'.format(PATH_FILE_OUT), 'RECREATE')
    FILE_OUT.cd()

    # Histogram for results
    #xaxis = array('d', tpl_dnde)
    xaxis = np.array(tpl_dnde_log+(2.*tpl_dnde_log[-1]-tpl_dnde_log[-2],), dtype=float)
    #xaxis_scaled = xaxis/SCALE_FLUX
    #yaxis = array('d', tpl_index)
    yaxis = np.array(tpl_index+(2.*tpl_index[-1]-tpl_index[-2],), dtype=float)
    dct_htg_likeresult = {}
    dct_htg_likeratio = {}
    dct_cvs_likeresult = {}
    dct_cvs_likeratio = {}
    for cla in TPL_CLASS:
        dct_htg_likeresult[cla] = ROOT.TH2D('htg_likeresult', 'Likelihood;log_{{10}}dN/dE at {0:1.2e} MeV;PWL index'.format(eref), len(tpl_dnde_log), xaxis, len(tpl_index), yaxis)
        dct_htg_likeratio[cla] = ROOT.TH2D('htg_likeratio', 'Likelihood Ratio;log_{{10}}dN/dE at {0:1.2e} MeV;PWL index'.format(eref), len(tpl_dnde_log), xaxis, len(tpl_index), yaxis)
        dct_cvs_likeresult[cla] = ROOT.TCanvas('cvs_likeresult_{0}'.format(cla), '{0} Likelihood'.format(cla), 750, 750)
        dct_cvs_likeratio[cla] = ROOT.TCanvas('cvs_likeratio_{0}'.format(cla), '{0} Likelihood Ratio'.format(cla), 750, 750)

    likelihood_ceil = math.exp(-HTG_OBS_ENR.Integral())
    for ienr in range(1, HTG_OBS_ENR.GetNbinsX()+1):
        ni = HTG_OBS_ENR.GetBinContent(ienr)
        likelihood_ceil = likelihood_ceil * math.pow(ni, ni)/math.factorial(ni)
    print 'Ideal maximum likelihood =', likelihood_ceil

    # Loop over dN/dE and PWL-index
    for (ix, dnde_log) in enumerate(tpl_dnde_log):
        dnde = 10**dnde_log
        print '===================='
        print 'dN/dE = {0:1.2e} at {1:1.1e} MeV'.format(dnde, eref)
        for (iy, idx_pl) in enumerate(tpl_index):
            print '--------------------'
            print 'PWL index = {0}'.format(idx_pl)
            lst_flux_itgl = ExtrapolateFlux.ExtrapolateFlux(eref, dnde, idx_pl, binw, enr_min, nbin, redshift)
            htg_flux = ROOT.TH1D('htg_flux', 'Integral flux', nbin, enr_min, enr_min+nbin*binw)
            for ibin in range(1, htg_flux.GetNbinsX()+1):
                htg_flux.SetBinContent(ibin, lst_flux_itgl[ibin-1])
                htg_flux.SetBinError(ibin, 0)
            str_fp = 'dNdE{0:0>12d}_PWL{1}{2:0>3d}'.format(int(dnde*1e20+0.5), "n" if idx_pl<0 else "p", int(idx_pl*100+0.5))
            suffix_fp = suffix + str_fp
            dct_htg_model = ModelPointSource.ModelPointSource(name, htg_flux, ra, dec, king, livetime, suffix_fp, nside, addregular)
            print dct_htg_model
            hs = ROOT.THStack('spectrum_{0}'.format(str_fp), 'log_{{10}}dN/dE={0:.2f} at {1} MeV, PWL-index={2:+f};log_{{10}}Energy [MeV];[counts]'.format(dnde_log, eref, idx_pl))
            hs.Add(HTG_OBS_ENR)
            for (icla,cla) in enumerate(TPL_CLASS):
                print cla
                htg_model = dct_htg_model[cla]
                htg_model.Rebin(NREBIN)
                htg_model.SetLineWidth(2)
                htg_model.SetLineColor(ROOT.kGray)
                htg_model.SetLineStyle(icla+1)
                htg_model.SetMarkerColor(ROOT.kGray)
                hs.Add(htg_model)
                likelihood = math.exp(-htg_model.Integral())
                for ienr in range(1, htg_model.GetNbinsX()+1):
                    print 'Energy range (model): 10^{0} - 10^{1}'.format(htg_model.GetXaxis().GetBinLowEdge(ienr), htg_model.GetXaxis().GetBinUpEdge(ienr))
                    print 'Energy range (observed): 10^{0} - 10^{1}'.format(HTG_OBS_ENR.GetXaxis().GetBinLowEdge(ienr), HTG_OBS_ENR.GetXaxis().GetBinUpEdge(ienr))
                    mi = htg_model.GetBinContent(ienr)
                    ni = HTG_OBS_ENR.GetBinContent(ienr)
                    likelihood = likelihood * math.pow(mi, ni)/math.factorial(ni)
                dct_htg_likeresult[cla].SetBinContent(ix+1, iy+1, likelihood)
                dct_htg_likeratio[cla].SetBinContent(ix+1, iy+1, likelihood/likelihood_ceil)
            FILE_OUT.cd()
            hs.Write()
            del dct_htg_model
            del htg_flux
            
    FILE_OUT.cd()
    likelihood_max = 0.0
    likelihood_temp = 0.0
    xlocmax = ROOT.Long()
    ylocmax = ROOT.Long()
    zlocmax = ROOT.Long()
    for cla in TPL_CLASS:
        dct_htg_likeresult[cla].Write()
        dct_htg_likeratio[cla].GetZaxis().SetRangeUser(0.05, 1.0)
        dct_htg_likeratio[cla].Write()
        dct_cvs_likeresult[cla].cd()
        dct_cvs_likeresult[cla].SetLogz()
        dct_htg_likeresult[cla].Draw("colz")
        likelihood_max = dct_htg_likeresult[cla].GetMaximum()
        dct_htg_likeresult[cla].GetMaximumBin(xlocmax, ylocmax, zlocmax)
        print '===== Maximum likelihood ====='
        print 'dNdE =', dct_htg_likeresult[cla].GetXaxis().GetBinLowEdge(xlocmax), 'at', eref, 'MeV'
        print 'PWL-index =', dct_htg_likeresult[cla].GetYaxis().GetBinLowEdge(ylocmax)
        dct_htg_likeresult[cla].GetZaxis().SetRangeUser(0.001*dct_htg_likeresult[cla].GetMaximum(), dct_htg_likeresult[cla].GetMaximum())
        dct_cvs_likeresult[cla].Write()
        dct_cvs_likeratio[cla].cd()
        dct_cvs_likeratio[cla].SetLogz()
        dct_htg_likeratio[cla].Draw("colz")
        dct_cvs_likeratio[cla].Write()
    return dct_htg_likeresult


@click.command()
@click.argument('name', type=str)
@click.argument('ra', type=float)
@click.argument('dec', type=float)
@click.argument('observed', type=str)
@click.argument('livetime', type=str)
@click.option('--king', type=str, default='/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/Dispersion/AG_dispersion.root')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--nside', '-n', type=int, default=256)
@click.option('--addregular', type=str, default=None, help='Set a file containing htgEaRegular')
@click.option('--redshift', '-z', type=float, default=0.)
@click.option('--dndeprimin', type=float, default=-20.)
@click.option('--dndeprimax', type=float, default=-10.)
@click.option('--indexprimin', type=float, default=-3.)
@click.option('--indexprimax', type=float, default=6.)
def main(name, ra, dec, observed, king, livetime, suffix, nside, redshift, addregular, dndeprimin, dndeprimax, indexprimin, indexprimax):

    FACTOR_RANGE_INDEX = 10
    INDEX_MIN = int(-3.0*FACTOR_RANGE_INDEX-0.5)
    INDEX_MAX = int(6.0*FACTOR_RANGE_INDEX+0.5)
    INDEX_STEP = int(0.2*FACTOR_RANGE_INDEX+0.5)
    TPL_INDEX = tuple([ float(x)/FACTOR_RANGE_INDEX for x in range(INDEX_MIN, INDEX_MAX+1, INDEX_STEP)])

    FACTOR_RANGE_FLUX = 10
    FLUX_LOG_MIN = int(-20.0*FACTOR_RANGE_FLUX-0.5)
    FLUX_LOG_MAX = int(-10.0*FACTOR_RANGE_FLUX-0.5)
    FLUX_LOG_STEP = int(0.2*FACTOR_RANGE_FLUX+0.5)
    #TPL_FLUX = tuple([ 10**(float(x)/FACTOR_RANGE_FLUX) for x in range(FLUX_LOG_MIN, FLUX_LOG_MAX+1, FLUX_LOG_STEP)])
    TPL_FLUX = tuple([ float(x)/FACTOR_RANGE_FLUX for x in range(FLUX_LOG_MIN, FLUX_LOG_MAX+1, FLUX_LOG_STEP)]) # In log scale
    if suffix!="":
        suffix = "_" + suffix
    LoopSEDLikelihood(name, observed, TPL_FLUX, TPL_INDEX, 10000., ra, dec, king, livetime, suffix, nside, redshift, addregular, dndeprimin, dndeprimax, indexprimin, indexprimax)


if __name__ == '__main__':
    main()
