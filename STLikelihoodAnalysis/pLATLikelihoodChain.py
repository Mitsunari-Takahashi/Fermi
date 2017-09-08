#!/usr/bin/env python
"""Module for LAT likelihood analysis chain of pLATLikelihoodConfig.py.
"""
import sys
import os
import os.path
import logging
import numpy as np
#import pickle
#import datetime
#from array import array
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
#import subprocess
import gt_apps as my_apps
import pyLikelihood
from UnbinnedAnalysis import *
from BinnedAnalysis import *
import matplotlib as mpl
import matplotlib.pyplot as plt
#from astropy.io import fits
#from sympy import *
#from scipy import integrate
#from fermipy.utils import get_parameter_limits
import copy
sys.path.append('/nfs/farm/g/glast/u/mtakahas/python_packages')
from make3FGLxml import *
import pLATLikelihoodConfig
from pLsList import ls_list
import pMETandMJD
from FindCrossEarthlimb import find_cross_earthlimb
from FindGoodstatPeriods import find_goodstat_periods, get_entries
from DownloadFermiData import download_fermi_data_grb
import ReadLTFCatalogueInfo
from STLikelihoodAnalysis import get_module_logger


##### Logger #####
logger = get_module_logger(__name__)
      

##### Analysis Chain Classes #####
class AnalysisChain:
    def __init__(self, lst_analysis_config):
        self.configs = lst_analysis_config
        self.nanalyses = len(self.configs)


class GRBExtrapolateSpectrum():
    def __init__(self, name, phase, emin_fitted, emax_fitted, emin_extrapolated, emax_extrapolated, tstop=10000., deg_roi=12., zmax=100., suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LTF):

        # Target GRB
        grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue)

        # Energy setups
        self.emin_fitted = emin_fitted
        self.emax_fitted = emax_fitted
        self.emin_extrapolated = emin_extrapolated
        self.emax_extrapolated = emax_extrapolated

        self.energybins = []

        # Analysis instance
        self.analysis_fit = pLATLikelihoodConfig.GRBConfig(target=grb, phase=phase, tstop=tstop, emin=self.emin_fitted, emax=self.emax_fitted, deg_roi=deg_roi, zmax=zmax, suffix=suffix)

        self.analysis_extrapolated = pLATLikelihoodConfig.GRBConfig(target=grb, phase=phase, tstop=tstop, emin=min(self.emin_fitted, self.emin_extrapolated), emax=max(self.emax_fitted, self.emax_extrapolated), deg_roi=deg_roi, zmax=zmax, suffix=suffix)

        AnalysisChain.__init__(self, [analysis_fit, analysis_extrapolated])


    def setup_fit(self):
        self.analysis_fit.setup()


    def setup_extrapolate(self):
        self.analysis_extrapolated.set_directories()
        self.analysis_extrapolated.download()
        self.analysis_extrapolated.filter()
        self.analysis_extrapolated.maketime()
        self.analysis_extrapolated.livetime()
        self.analysis_extrapolated.exposure()
        logger.debug('Path of reffered model: {0}'.format(self.analysis_fit.path_model_xml_new))
        self.analysis_extrapolated.use_external_model(self.analysis_fit.path_model_xml_new)
        self.analysis_extrapolated.diffuse_responses()
        

    def fit(self, redo=True):
        self.analysis_fit.fit(bredo=redo)


    def set_likelihood_extrapolate(self):
        self.analysis_extrapolated.set_likelihood()
        for i in range(self.analysis_extrapolated.like.nFreeParams()):
            self.analysis_extrapolated.like.freeze(i)


    def plot_extrapolated_count_spectrum(self):
        self.analysis_extrapolated.plot_countspectra_fitted()


    def eval_deviation(self):
        # Energy bins
        ebins = self.analysis_extrapolated.like.energies
        nebins = len(ebins)-1

        # Model count
        x_cspec_fit, y_model_all, y_model_target, y_model_others = self.analysis_extrapolated.count_axes()

        # Eval error
        y_model_err = np.zeros_like(y_model_all)
        for ie in range(nebins):
            flux, flux_err = self.analysis_fit.eval_flux_and_error_total(emin=ebins[ie], emax=ebins[ie+1])
            flux_frac_err = flux_err / flux
            y_model_err[ie] = y_model_all[ie] * flux_frac_err

        fig, axes = self.analysis_extrapolated.plot_countspectra_fitted()
        logger.debug(fig)
        axes[0].fill_between(x_cspec_fit, y_model_all+y_model_err, y_model_all-y_model_err, alpha=0.2, color='b', label='Fitting uncertainty')
        fig.savefig("{0}/Extrapolated_count_spectrum_{1}{2}.png".format(self.analysis_extrapolated.dir_work, self.analysis_extrapolated.target.name, self.analysis_extrapolated.suffix))


# class LightCurves(AnalysisChain):
#     def __init__(self, met0, emin, emax, deg_roi=12., zmax=100.,):
