#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi, sqrt
import matplotlib as mpl
import matplotlib.pyplot as plt
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TFile, TDirectory, TH1F, TH2F, TH3F
ROOT.gROOT.SetBatch()
from pColor import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL

##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


class EventClass:
    def __init__(self, name, cut_string, performance_file):
        self.cut = cut_string
        self.perf_hists = {}
        for q in ('acc_cth', 'acc_psf', 'psf_q68', 'psf_q95', 'psf_cth_q68', 'psf_cth_q95', 'acc_edisp', 'edisp_q68', 'edisp_q95', 'edisp_cth_q68', 'edisp_cth_q95'):
            self.perf_hists[q] = performance_file.Get('_'.join([q, 'hist']))


    def lookup_psf(self, loge, costh, quantile=68):
        return self.perf_hists['psf_cth_q{q}'.format(q=quantile)].Interpolate(loge, costh)


    def lookup_edisp(self, loge, costh, quantile=68):
        return self.perf_hists['edisp_cth_q{q}'.format(q=quantile)].Interpolate(loge, costh)
        

    def calc_count(self, source_model_hist, livetime_hist, emin, emax):
        # source_model: vs. energy
        # livetime: vs. costheta vs. energy
        # acceptance: vs. costheta vs. energy
        # exposure: vs. costheta vs. energy

        exp_hist = livetime_hist.Clone('exposure')
        exp_hist.Multiply(self.perf_hists['acc_cth'])

        count_cth_hist = TH2F('count_cth', 'Predicted count', source_model_hist.GetXaxis(), exp_hist.GetYaxis())
        for ix in range(1, count_cth_hist.GetXaxis().GetNbins()+1):
            nebin_lo = exp_hist.GetXaxis().FindBin(count_cth_hist.GetXaxis().GetBinLowEdge(ix))
            nebin_up = exp_hist.GetXaxis().FindBin(count_cth_hist.GetXaxis().GetBinUpEdge(ix))
            flux = source_model_hist.GetBinContent(ix)
            flux_err = source_model_hist.GetBinError(ix)

            for jy in range(1, count_cth_hist.GetYaxis().GetNbins()+1):
                # Exposure except for the edge bins
                exp_err_nonedge = ROOT.Double(0)
                exp_nonedge = exp_hist.IntegralAndError(nebin_lo+1, nebin_up-1, jy, jy, exp_err)

                # Exposure of the lowest bin
                wfactor_loedge = (exp_hist.GetXaxis().GetBinUpEdge(nebin_lo)-count_cth_hist.GetXaxis().GetBinLowEdge(ix)) / exp_hist.GetXaxis().GetBinWidth(nebin_lo)
                exp_loedge = exp_hist.GetBinContent(nebin_lo, jy) * wfactor_loedge
                exp_err_loedge = exp_hist.GetBinError(nebin_lo, jy) * wfactor_loedge

                # Exposure of the highest bin
                wfactor_hiedge = (count_cth_hist.GetXaxis().GetBinUpEdge(ix) - exp_hist.GetXaxis().GetBinLowEdge(nebin_hi)) / exp_hist.GetXaxis().GetBinWidth(nebin_hi)
                exp_hiedge = exp_hist.GetBinContent(nebin_hi, jy) * wfactor_hiedge
                exp_err_hiedge = exp_hist.GetBinError(nebin_hi, jy) * wfactor_hiedge

                # Total exposure
                exp = exp_nonedge + exp_loedge + exp_hiedge
                exp_err = sqrt(exp_err_nonedge**2 + exp_loedge**2 + exp_err_hiedge**2)

                cnt = flux * exp
                cnt_err = 
                count_cth_hist.SetBinContent(ix, jy, cnt)
                count_cth_hist.SetBinError(ix, jy, cnt_err)
        return count_cth


@click.command()
@click.argument('name', type=str)
@click.argument('dst', nargs=-1)
@click.option('--suffix', type=str, default='')
@click.option('--values', type=(str, int))
@click.option('--values', multiple=True)
@click.option('--language', type=click.Choice(['Japanese', 'English']))
@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(name, sed, ra, dec, king, acceptance, livetime, suffix, nside):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

if __name__ == '__main__':
    main()
