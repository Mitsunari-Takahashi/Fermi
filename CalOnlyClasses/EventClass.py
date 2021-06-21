#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import matplotlib as mpl
import matplotlib.pyplot as plt
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
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
        

    def calc_count(self, source_model, livetime_hist, emin, emax):
        # source_model: vs. energy
        # livetime: vs. costheta vs. energy
        # acceptance: vs. costheta vs. energy
        # exposure: vs. costheta vs. energy

        exp_hist = livetime_hist.Clone('exposure')
        exp_hist.Multiply(self.perf_hists['acc_cth'])

        eaxis = self.perf_hists['acc_cth'].Clone('count_hist').GetXaxis()
        

        count_hist = self.perf_hists['acc_cth'].Clone('count_hist')
        


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
