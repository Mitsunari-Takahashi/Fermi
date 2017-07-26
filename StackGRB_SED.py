#!/usr/bin/env python

import os
import os.path
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import fermipy.utils
from fermipy.gtanalysis import GTAnalysis
#import GtApp
#import FluxDensity
#from LikelihoodState import LikelihoodState
#from fermipy.gtutils import BinnedAnalysis, SummedLikelihood
#import BinnedAnalysis as ba
#import pyLikelihood as pyLike
import ROOT
ROOT.gROOT.SetBatch()
from ROOT import gStyle
from array import array
import numpy as np
import math
from math import log10
from scipy import interpolate
from pLsList import ls_list
#from pReadGBMCatalogueInfo import ReadGBMCatalogueOneLine
import click


def stack_grb_sed(lst_sed_npy, path_outdir, suffix):
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)
    if suffix != '':
        suffix = '_' + suffix
    pathout_base = '{0}/StackedGRBs{1}'.format(path_outdir, suffix)
    fileout = ROOT.TFile(pathout_base+'.root', 'RECREATE')
    gStyle.SetPalette(51)
    gStyle.SetOptStat(0)
    lst_htg_loglike = []
    lst_htg_dlikelihood = []
    lst_htg_nlikelihood = []
    lst_htg_nlike_scaled = []
    logf_min = -8
    logf_max = 0
    logf_binw = 0.01
    ax_logf = np.arange(logf_min, logf_max, logf_binw) 
    nbinf = len(ax_logf)
    ax_logf_edges = np.arange(logf_min-logf_binw/2., logf_max+logf_binw/2., logf_binw)

    ax_logf_scaled_edges = np.arange(logf_min-logf_binw/2.+4, logf_max+logf_binw/2.+4, logf_binw) 

    ax_logf_edges_display = ax_logf_scaled_edges[49:-49]
    print ax_logf_edges_display[0], ax_logf_edges_display[1], ax_logf_edges_display[-2], ax_logf_edges_display[-1]
    ax_logf_edges_display[0] = ax_logf_scaled_edges[0]
    ax_logf_edges_display[-1] = ax_logf_scaled_edges[-1]
    print ax_logf_edges_display[0], ax_logf_edges_display[1], ax_logf_edges_display[-2], ax_logf_edges_display[-1]

    ax_loge_intrinsic_edges = np.array([log10(x) for x in [316.228, 1778.28, 5623.41, 17782.8, 56234.1, 177828.0]])
    htg_nlike_scaled_stacked = ROOT.TH2D('nlike_scaled_stacked', '{0} GRBs stacked;log_{{10}}Energy [MeV];Relative flux to one in {1:.1f}-{2:.1f} MeV'.format(len(lst_sed_npy), pow(10, ax_loge_intrinsic_edges[0]), pow(10, ax_loge_intrinsic_edges[1])), len(ax_loge_intrinsic_edges)-1, ax_loge_intrinsic_edges, len(ax_logf_scaled_edges)-1, ax_logf_scaled_edges)

    for (ised, path_sed) in enumerate(lst_sed_npy):
        grbname = os.path.basename(path_sed)[7:16]
        print '=====', grbname, '====='
        print path_sed
        fileout.mkdir(grbname)
        fileout.cd(grbname)
        sed = np.load(path_sed).flat[0]
        ax_loge_edges = sed['loge_min']
        ax_loge_edges = np.append(ax_loge_edges, sed['loge_max'][-1])
        #ax_loge_edges = array('d', nda_loge_edges)
        print '.'
        lst_htg_loglike.append(ROOT.TH2D('loglike_{0:0>3d}'.format(ised), grbname+';log_{10}Energy [MeV];E^{2}dN/dE [MeV cm^{-2} s^{-1}]', len(ax_loge_edges)-1, ax_loge_edges, len(ax_logf_edges)-1, ax_logf_edges))
        print '.'
        lst_htg_dlikelihood.append(ROOT.TH2D('dlikelihood_{0:0>3d}'.format(ised), grbname+';log_{10}Energy [MeV];E^{2}dN/dE [MeV cm^{-2} s^{-1}]', len(ax_loge_edges)-1, ax_loge_edges, len(ax_logf_edges)-1, ax_logf_edges))
        print '.'
        lst_htg_nlikelihood.append(ROOT.TH2D('nlikelihood_{0:0>3d}'.format(ised), grbname+';log_{10}Energy [MeV];E^{2}dN/dE [MeV cm^{-2} s^{-1}]', len(ax_loge_edges)-1, ax_loge_edges, len(ax_logf_edges)-1, ax_logf_edges))
        print '.'
        lst_htg_nlike_scaled.append(ROOT.TH2D('nlike_scaled_{0:0>3d}'.format(ised), grbname+';log_{{10}}Energy [MeV] in rest frames;Relative flux to one in {0:.1f}-{1:.1f} MeV'.format(pow(10, ax_loge_edges[0]), pow(10, ax_loge_edges[1])), len(ax_loge_intrinsic_edges)-1, ax_loge_intrinsic_edges, len(ax_logf_scaled_edges)-1, ax_logf_scaled_edges))
        print '.'
        for ie in range(lst_htg_loglike[-1].GetXaxis().GetNbins()):
            e_ref = sed['e_ref'][ie]
            ref_e2dnde = sed['ref_e2dnde'][ie]
            e2dnde_best = sed['e2dnde'][ie]
            norm_scan = sed['norm_scan'][ie]
            dloglike_scan = sed['dloglike_scan'][ie]
            #loglike_scan = sed['loglike_scan'][ie]
            loglike_best = sed['loglike'][ie]

            m = norm_scan > 0
            e2dnde_scan = norm_scan[m] * ref_e2dnde
            loge2dnde_scan = np.log10(e2dnde_scan)
            logl = dloglike_scan[m]
            logl -= np.max(logl)
            try:
                fn = interpolate.interp1d(loge2dnde_scan, logl, fill_value='extrapolate')
                logli = fn(ax_logf)
            except:
                logli = np.interp(ax_logf, loge2dnde_scan, logl)
            for (iflux, flux) in enumerate(ax_logf):
                lst_htg_loglike[-1].Fill(lst_htg_loglike[-1].GetXaxis().GetBinCenter(ie+1), flux, logli[iflux])
                lst_htg_dlikelihood[-1].Fill(lst_htg_loglike[-1].GetXaxis().GetBinCenter(ie+1), flux, pow(10, logli[iflux]))
                #conf = 1.-fermipy.utils.twosided_dlnl_to_cl(-logli[iflux])
                #lst_htg_conf[-1].Fill(lst_htg_conf[-1].GetXaxis().GetBinCenter(ie+1), flux, conf)
            underlimit = 1.-fermipy.utils.onesided_dlnl_to_cl(-logli[0])
            print 'Under limit', underlimit
            overlimit = 1.-fermipy.utils.onesided_dlnl_to_cl(-logli[-1])
            print 'Over limit', overlimit
            lst_htg_dlikelihood[-1].SetBinContent(ie+1, 0, underlimit*lst_htg_dlikelihood[-1].Integral(ie+1, ie+1, 1, nbinf+1)/(1-underlimit-overlimit))
            lst_htg_dlikelihood[-1].SetBinContent(ie+1, lst_htg_dlikelihood[-1].GetYaxis().GetNbins()+1, overlimit*lst_htg_dlikelihood[-1].Integral(ie+1, ie+1, 1, nbinf+1)/(1-overlimit-underlimit))
            like_scale = lst_htg_dlikelihood[-1].Integral(ie+1, ie+1, 0, lst_htg_dlikelihood[-1].GetYaxis().GetNbins()+1)
            for jflux in range(len(ax_logf)+2):
                lst_htg_nlikelihood[-1].SetBinContent(ie+1, jflux, lst_htg_dlikelihood[-1].GetBinContent(ie+1, jflux)/like_scale)
            lst_htg_nlike_scaled[-1].SetBinContent(ie+1, 0, lst_htg_nlikelihood[-1].GetBinContent(ie+1, 0))
            lst_htg_nlike_scaled[-1].SetBinContent(ie+1, nbinf+1, lst_htg_nlikelihood[-1].GetBinContent(ie+1, nbinf+1))
            for (kflux, logf) in enumerate(ax_logf):
                if ie>1:
                    for (lflux, logfl) in enumerate(ax_logf):
                        logf_scaled = logf - logfl
                        wlike = lst_htg_nlikelihood[-1].GetBinContent(1, lflux+1)
                        lst_htg_nlike_scaled[-1].Fill(lst_htg_nlike_scaled[-1].GetXaxis().GetBinCenter(ie+1), logf_scaled, lst_htg_nlikelihood[-1].GetBinContent(ie+1, kflux+1)*wlike)
                else:
                    lst_htg_nlike_scaled[-1].Fill(lst_htg_nlike_scaled[-1].GetXaxis().GetBinCenter(ie+1), logf-np.log10(e2dnde_best), lst_htg_nlikelihood[-1].GetBinContent(ie+1, kflux+1))

            # nflux_best = lst_htg_prob[-1].GetYaxis().FindBin(np.log10(e2dnde_best))
            # print 'Best flux bin:', nflux_best
            # for jflux in range(nflux_best+1, lst_htg_prob[-1].GetYaxis().GetNbins()+1): # Above best flux value
            #     p_deriv = lst_htg_conf[-1].GetBinContent(ie+1, jflux) - lst_htg_conf[-1].GetBinContent(ie+1, jflux+1)
            #     lst_htg_prob[-1].SetBinContent(ie+1, jflux, p_deriv)
            # for kflux in range(1, nflux_best): # Below best flux value
            #     p_deriv = lst_htg_conf[-1].GetBinContent(ie+1, kflux) - lst_htg_conf[-1].GetBinContent(ie+1, kflux-1)
            #     lst_htg_prob[-1].SetBinContent(ie+1, kflux, p_deriv)
            # p_deriv = lst_htg_conf[-1].GetBinContent(ie+1, nflux_best) - lst_htg_prob[-1].Integral(ie+1, ie+1, 0, -1)
                
        lst_htg_loglike[-1].Write()
        lst_htg_dlikelihood[-1].Write()
        lst_htg_nlikelihood[-1].Write()
        lst_htg_nlike_scaled[-1].Write()
        htg_nlike_scaled_stacked.Add(lst_htg_nlike_scaled[-1])
        #lst_htg_conf[-1].Write()
        #lst_htg_prob[-1].Write()
            #llhMatrix[i, :] = logli
    fileout.cd()
    htg_nlike_scaled_stacked.Write()
    htg_nlike_scaled_stacked_forDisplay = ROOT.TH2D(htg_nlike_scaled_stacked.GetName()+'_forDisplay', '{0} GRBs stacked;log_{{10}}Energy [MeV] in rest frames;Relative flux to one in {1:.1f}-{2:.1f} MeV'.format(len(lst_sed_npy), pow(10, ax_loge_intrinsic_edges[0]), pow(10, ax_loge_intrinsic_edges[1])), len(ax_loge_intrinsic_edges)-1, ax_loge_intrinsic_edges, len(ax_logf_edges_display)-1, ax_logf_edges_display) #int((logf_max-logf_min)/logf_binw), logf_min+4, logf_max+4)
    #htg_nlike_scaled_stacked_forDisplay = htg_nlike_scaled_stacked.Clone(htg_nlike_scaled_stacked.GetName()+'_forDisplay')
    for le in range(1, htg_nlike_scaled_stacked.GetXaxis().GetNbins()+1):
        for lf in range(1, htg_nlike_scaled_stacked.GetYaxis().GetNbins()+1):
            htg_nlike_scaled_stacked_forDisplay.Fill(htg_nlike_scaled_stacked.GetXaxis().GetBinCenter(le), htg_nlike_scaled_stacked.GetYaxis().GetBinCenter(lf), htg_nlike_scaled_stacked.GetBinContent(le, lf))
        htg_nlike_scaled_stacked_forDisplay.Fill(htg_nlike_scaled_stacked.GetXaxis().GetBinCenter(le), htg_nlike_scaled_stacked.GetYaxis().GetBinCenter(1), htg_nlike_scaled_stacked.GetBinContent(le, 0))
        htg_nlike_scaled_stacked_forDisplay.Fill(htg_nlike_scaled_stacked.GetXaxis().GetBinCenter(le), htg_nlike_scaled_stacked.GetYaxis().GetBinCenter(nbinf+1), htg_nlike_scaled_stacked.GetBinContent(le, nbinf+1))

    htg_nlike_scaled_stacked_forDisplay.Write()
    cdisplay = ROOT.TCanvas('cdisplay', 'Stacked relative SED', 800, 800)
    cdisplay.SetLogz()
    cdisplay.cd()
    htg_nlike_scaled_stacked_forDisplay.Draw("colz")
    htg_nlike_scaled_stacked_forDisplay.SetMinimum(1e-6)
    cdisplay.Write()
    cdisplay.SaveAs(pathout_base+'.pdf')
        

@click.command()
@click.argument('inputs', type=str)
@click.option('-s', '--suffix', type=str, default='')
@click.option('--outpath', '-o', default='.')
def main(inputs, outpath, suffix):
    with open(inputs, "r") as filein:
        str_paths = filein.read()
        input_paths = str_paths.split('\n')[:-1]
        print input_paths
        stack_grb_sed(input_paths, outpath, suffix)
        #stack_grb_sed(ls_list(inputs), outpath, suffix)


if __name__ == '__main__':
    main()
