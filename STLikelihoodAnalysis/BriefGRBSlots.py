#!/usr/bin/env python
"""Module for making light curves of LAT data.
The main class LightCurve is a chain of another module pLATLikelihoodConfig.py.
The authour: Mitsunari Takahashi
* Version: 0.0 (2017.10.31)
"""
import sys
import os
import os.path
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import logging
import pickle
import subprocess
import datetime
import numpy as np
import itertools
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
from sympy import *
from scipy import integrate
from astropy.io import fits
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import STLikelihoodAnalysis.pLATLikelihoodConfig as pLATLikelihoodConfig
from FindGoodstatPeriods import find_goodstat_periods, get_entries, get_event_time_and_energy
from FindCrossEarthlimb import find_cross_earthlimb
from STLikelihoodAnalysis import get_module_logger
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
import pickle_utilities
import pMatplot
from LightCurve import LightCurve, LightCurveGRB


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

##### VERSION OF THIS MACRO #####
VERSION = 0.1 # 2017.12.09
#VERSION = 0.0 # 2017.10.31


##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6

##### Catalogue #####
CATALOGUE_LAT = ReadLATCatalogueInfo.open_table(pLATLikelihoodConfig.GRB_CATALOGUE_LAT)
CATALOGUE_GBM = ReadGBMCatalogueInfo.open_table()


def brief_grb_slots(name, emin, emax, roi, suffix, refit, force, outdir, index, normanchor, tmin=0, tmax=100000):
    tb_one = ReadLATCatalogueInfo.select_one_by_name(CATALOGUE_LAT, name, CATALOGUE_GBM)
    lc = LightCurveGRB(name=name, wholephase='afterglow', tmin=tmin, tmax=tmax, emin=emin, emax=emax, deg_roi=roi, ngoodstat=0, rlim_goodstat=0, suffix=suffix, grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, refit=refit, force=force, outdir=None, phase='briefslots', spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2, 'Scale':emin}) #norm=1e-10, index=index, scalefactor=tb_one['GBM']['FLUENCE'])

    lognorm_anc = 2.5-2.*np.log10(emin)
    norms = 10**np.linspace(lognorm_anc-1., lognorm_anc+1., 101)*normanchor #-7, 0, 701)
    lcindices = np.linspace(-0.655, -1.555, 91) #np.array([-1.0, -1.3, -1.6])
    lc.setup()
    like_results = lc.scan_parameters(norms, lcindices, tnorm=10., rescaler=tb_one['GBM']['FLUENCE'])
    lc.pickle(like_results)


@click.command()
@click.option('--namemin', type=str, default='000000000')
@click.option('--namemax', type=str, default='999999999')
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--roi', type=float, default=12.)
@click.option('--specindex', type=float, default=-2.0)
@click.option('--normanchor', '-a', type=float, default=1.0)
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--outdir', '-o', type=str, default='')
@click.option('--bsub', '-b', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(namemin, namemax, emin, emax, roi, refit, force, suffix, outdir, specindex, normanchor, bsub, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    if bsub==True:
        tb_lat = ReadLATCatalogueInfo.select_by_name(CATALOGUE_LAT, namemin, namemax, CATALOGUE_GBM)
        tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
        tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
        tb_lat = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)

        if len(tb_lat)<1:
            print 'No GRBs.'
            return 1
        for (irow , row) in enumerate(tb_lat):
            name = row['GRBNAME']
            print '##### No.{0} GRB{1} #####'.format(irow, name)
            if not os.path.exists(name):
                os.mkdir(name)
            acmd = ['bsub', '-o','{0}/GRB{0}_briefslots{1}_E{2:0>7.0f}-{3:0>7.0f}MeV.log'.format(name, suffix if suffix=='' else '_'+suffix, emin, emax), '-J','bs{0}'.format(name[:-3]), '-W','900', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/BriefGRBSlots.py', '--emin', str(emin), '--emax', str(emax), '--normanchor', str(normanchor), '-s', suffix, '--roi', str(roi), '--namemin', name, '--outdir', outdir, '--specindex', str(specindex)]
            if force==True:
                acmd.append('--force')
            if refit==True:
                acmd.append('--refit')
            print acmd
            subprocess.call(acmd)
        return 0
    else:
        ##### Find GRBs #####
        tb_lat = ReadLATCatalogueInfo.select_one_by_name(CATALOGUE_LAT, namemin, CATALOGUE_GBM)

        ##### Main function #####
        brief_grb_slots(name=namemin, emin=emin, emax=emax, roi=roi, suffix=suffix, refit=refit, force=force, outdir=outdir, index=specindex, normanchor=normanchor, tmin=0, tmax=100000)


if __name__ == '__main__':
    main()
