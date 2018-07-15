#!/Usr/bin/env python
"""Modules for LAT likelihood analysis in python style.
"""
import sys
import os
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import os.path
import shutil
import logging
import itertools
import numpy as np
import scipy.misc
from scipy import interpolate, integrate
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
import ROOT
from collections import OrderedDict
import gt_apps as my_apps
import pyLikelihood
from UnbinnedAnalysis import *
from BinnedAnalysis import *
import FluxDensity
from LikelihoodState import LikelihoodState
import SummedLikelihood
from astropy.io import fits
from astropy.table import Table, Column
from fermipy.utils import get_parameter_limits
import matplotlib as mpl
import matplotlib.pyplot as plt
#from ROOT import TMath.Prob as Prob
import copy
sys.path.append('/nfs/farm/g/glast/u/mtakahas/python_packages')
from make3FGLxml import *
from pLsList import ls_list
import pMETandMJD
from FindCrossEarthlimb import find_cross_earthlimb
from FindGoodstatPeriods import get_entries, get_entries_roi
from DownloadFermiData import download_fermi_data_grb
import ReadLTFCatalogueInfo
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
from STLikelihoodAnalysis import get_module_logger
import GetGTI
from pCommon import MEVtoERG
import pMatplot


##### Logger #####
logger = get_module_logger(__name__)


##### Matplotlib #####
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams["font.size"] = 15


##### d-loglike value for certain signicicance cuts #####
# Significance(sigma): [ Doubled d-loglike for NDF=1,2,3 ]
TABLE_DLOGLIKE_SIGNIFICANCE = Table({'1.0': [1.00, 2.30, 3.53],
                                     '2.0': [4.00, 6.18, 8.03],
                                     '3.0': [9.00, 11.83, 14.16]
                                     })

##### PATH of Catalogue #####
GRB_CATALOGUE_LTF = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'
GRB_CATALOGUE_LAT = "/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LATBurstCatalogue.xml"
GRB_CATALOGUE_GBM = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/GBMcatalogue20171005.fits'


##### LAT IRFs #####
DCT_EVCLASSES = {8:'P8R2_TRANSIENT020E_V6', 
                 16:'P8R2_TRANSIENT020_V6',
                 32:'P8R2_TRANSIENT010E_V6',
                 64:'P8R2_TRANSIENT010E_V6',
                 128:'P8R2_SOURCE_V6',
                 256:'P8R2_CLEAN_V6',
                 512:'P8R2_ULTRACLEAN_V6',
                 1024:'P8R2_ULTRACLEANVETO_V6'}

DICT_EVCLASSES_BIT_CALONLY = {'R100':(4096, 8192, 16384, 32768),
                              'R030':(8192, 16384, 32768),
                              'R010':(16384, 32768),
                              'R010':(32768,)}

DICT_SPECTRALPARS_SET = {'PowerLaw': set(['Prefactor', 'Index', 'Scale']),
                         'PowerLaw2': set(['Integral', 'Index', 'LowerLimit', 'UpperLimit']),
                         'EblAtten::PowerLaw2': set(['Integral', 'Index', 'LowerLimit', 'UpperLimit', "tau_norm", "redshift", "ebl_model"]),
                         'ScaleFactor::PowerLaw2': set(['Integral', 'Index', 'LowerLimit', 'UpperLimit', 'ScaleFactor']),
                         'ExpCutoff': set(['Prefactor', 'Index', 'Scale', 'Ebreak', 'P1', 'P2', 'P3']),
                         'EblAtten::ExpCutoff': set(['Prefactor', 'Index', 'Scale', 'Ebreak', 'P1', 'P2', 'P3', "tau_norm", "redshift", "ebl_model"]),
                         'BrokenPowerLaw': set(['Prefactor', 'Index1', 'Index2', 'BreakValue', 'LowerLimit', 'UpperLimit']),
                         'BrokenPowerLaw2': set(['Integral', 'Index1', 'Index2', 'BreakValue', 'LowerLimit', 'UpperLimit'])
                         }

DICT_SPECTRALPARS_SCALE = {'PowerLaw':{'Prefactor':1, 'Index':1, 'Scale':1},
                           'PowerLaw2': {'Integral':1, 'Index':1, 'LowerLimit':1, 'UpperLimit':1},
                           'EblAtten::PowerLaw2': {'Integral':1, 'Index':1, 'LowerLimit':1, 'UpperLimit':1, 'tau_norm':1, 'redshift':1, 'ebl_model':1},
                           'ScaleFactor::PowerLaw2': {'Integral':1, 'Index':1, 'LowerLimit':1, 'UpperLimit':1, 'ScaleFactor':1},
                           'ExpCutoff': {'Prefactor':1, 'Index':1, 'Scale':1, 'Ebreak':1, 'P1':1, 'P2':1, 'P3':1},
                           'EblAtten::ExpCutoff': {'Prefactor':1, 'Index':1, 'Scale':1, 'Ebreak':1, 'P1':1, 'P2':1, 'P3':1, 'tau_norm':1, 'redshift':1, 'ebl_model':1},
                           'BrokenPowerLaw': {'Prefactor':1, 'Index1':1, 'Index2':1, 'BreakValue':1},
                           'BrokenPowerLaw2': {'Integral':1, 'Index1':1, 'Index2':1, 'BreakValue':1, 'LowerLimit':1, 'UpperLimit':1}
                           }

DICT_SPECTRALPARS_BOUNDS = {'PowerLaw':{'Prefactor':(0.,1.), 'Index':(-9.,3.), 'Scale':(100.,100000.)},
                           'PowerLaw2': {'Integral':(0.,1.), 'Index':(-9.,3.), 'LowerLimit':(30.,100000.), 'UpperLimit':(100.,1000000.)},
                           'EblAtten::PowerLaw2': {'Integral':(0.,1.), 'Index':(-9.,3.), 'LowerLimit':(30.,100000.), 'UpperLimit':(100.,1000000.), 'tau_norm':(0.,10.), 'redshift':(0.,10.), 'ebl_model':(0,7)},
                           'ScaleFactor::PowerLaw2': {'Integral':(0.,1.), 'Index':(-9.,3.), 'LowerLimit':(30.,100000.), 'UpperLimit':(100.,1000000.), 'ScaleFactor':(0.,100.)},
                           'ExpCutoff': {'Prefactor':(0.,1.), 'Index':(-9.,3.), 'Scale':(100.,100000.), 'Ebreak':(100.,100000.), 'P1':(100.,100000.), 'P2':(100.,100000.), 'P3':(100.,100000.)},
                           'EblAtten::ExpCutoff': {'Prefactor':(0.,1.), 'Index':(-9.,3.), 'Scale':(100.,100000.), 'Ebreak':(10.,100000.), 'P1':(10.,100000.), 'P2':(0.,100000.), 'P3':(0.,100000.), 'tau_norm':(0.,10.), 'redshift':(0.,10.), 'ebl_model':(0,7)},
                           'BrokenPowerLaw': {'Prefactor':(0.,1.), 'Index1':(-9.,3.), 'Index2':(-9.,3.), 'BreakValue':(100.,100000.)},
                           'BrokenPowerLaw2': {'Integral':(0.,1.), 'Index1':(-9.,3.), 'Index2':(-9.,3.), 'BreakValue':1, 'LowerLimit':(30.,100000.), 'UpperLimit':(100.,1000000.)}
                           }

DICT_SPECTRALPARS_FIXED = {'PowerLaw': ['Scale'],
                           'PowerLaw2': ['LowerLimit', 'UpperLimit'],
                           'EblAtten::PowerLaw2': ['LowerLimit', 'UpperLimit', 'tau_norm', 'redshift', 'ebl_model'],
                           'ScaleFactor::PowerLaw2': ['LowerLimit', 'UpperLimit', 'ScaleFactor'],
                           'ExpCutoff': ['Scale', 'P2', 'P3'],
                           'EblAtten::ExpCutoff': ['Scale', 'P2', 'P3', 'tau_norm', 'redshift', 'ebl_model'],
                           'BrokenPowerLaw': ['BreakValue'],
                           'BrokenPowerLaw2': ['LowerLimit', 'UpperLimit']
                           }


##### Target Classes #####
class _AstroTarget:
    def __init__(self, name, spatialmodel, spectraltype, spectralpars, spectralpars_scale, spectralpars_fixed, spectralpars_bounds=None, met0=0, redshift=np.nan):
        self.name = name
        self.spatialmodel = spatialmodel
        self.spectraltype = spectraltype
        if not set(spectralpars.keys()) == DICT_SPECTRALPARS_SET[self.spectraltype]:
            logger.error('Spectral parameters are NOT correct for {0}!!!'.format(self.spectraltype))
            logger.error(spectralpars.keys())
            logger.error('Correct parameters:')
            logger.error(DICT_SPECTRALPARS_SET[self.spectraltype])
            sys.exit(1)

        # Name of Normalization parameter
        if spectraltype in ('PowerLaw', 'BrokenPowerLaw', 'ExpCutoff', 'EblAtten::ExpCutoff'):
            self.norm_name = 'Prefactor' 
        elif spectraltype in ('PowerLaw2', 'BrokenPowerLaw2', 'ScaleFactor::PowerLaw2', 'EblAtten::PowerLaw2'):
            self.norm_name = 'Integral' 
        else:
            logger.critical('No normalization name is found.')
            sys.exit(1)

        self.spectralpars = spectralpars
        self.spectralpars_scale = spectralpars_scale
        self.spectralpars_fixed = spectralpars_fixed
        self.spectralpars_bounds = spectralpars_bounds
        #self.path_pickle = path_pickle
        self.met0 = met0
        self.redshift = redshift
        logger.debug('Target name: {0}'.format(self.name))


class PointTarget(_AstroTarget):
    def __init__(self, name, ra, dec, spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000.}, spectralpars_scale={'Prefactor':1, 'Index':1, 'Scale':1}, spectralpars_fixed=['Scale'], spectralpars_bounds=None, loc_err=0., met0=0, redshift=np.nan):
        _AstroTarget.__init__(self, name, 'PointSource', spectraltype, spectralpars, spectralpars_scale, spectralpars_fixed, spectralpars_bounds=spectralpars_bounds, met0=met0, redshift=redshift)
        logger.debug('MET0 = {0}'.format(self.met0))
        logger.debug('Target name: {0}'.format(self.name))
        self.ra = ra
        self.dec = dec
        self.loc_err = loc_err


class GRBTarget(PointTarget):
    def __init__(self, name, path_catalogue=GRB_CATALOGUE_LAT, spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000.}, spectralpars_fixed=None, spectralpars_bounds=None, redshift=0., synccutoff=None):

        tb_gbm = ReadGBMCatalogueInfo.open_table(1, GRB_CATALOGUE_GBM)
        self.path_catalogue = path_catalogue
        if path_catalogue[-5:] == '.fits':
            tb_lat = ReadLTFCatalogueInfo.open_table(1, self.path_catalogue)
            tb_one = ReadLTFCatalogueInfo.select_one_by_name(tb_lat, name)
        elif path_catalogue[-4:] == '.xml':
            tb_lat = ReadLATCatalogueInfo.open_table(GRB_CATALOGUE_LAT)
            tb_one = ReadLATCatalogueInfo.select_one_by_name(tb_lat, name, tb_gbm)
            self.table_grb_catalogue = tb_one
        else:
            logger.critical('Catalogue {0} is not in readable format!!!'.format(path_catalogue))
        self.synccutoff = synccutoff

        spectralpars_scale = DICT_SPECTRALPARS_SCALE[spectraltype]
        spectralpars_fixed = spectralpars_fixed if spectralpars_fixed is not None else DICT_SPECTRALPARS_FIXED[spectraltype]
        if self.synccutoff is not None and self.synccutoff is not False and self.synccutoff==self.synccutoff:
            for sfp in ['Ebreak', 'P1']:
                if not sfp in spectralpars_fixed:
                    spectralpars_fixed.append(sfp)
        spectralpars_bounds = spectralpars_bounds if spectralpars_bounds is not None else DICT_SPECTRALPARS_BOUNDS[spectraltype]

        PointTarget.__init__(self, name, float(tb_one['RA']), float(tb_one['DEC']), spectraltype=spectraltype, spectralpars=spectralpars, spectralpars_scale=spectralpars_scale, spectralpars_fixed=spectralpars_fixed, spectralpars_bounds=spectralpars_bounds, met0=tb_one['TRIGGER_TIME'], redshift=redshift) 

        if path_catalogue[-5:] == '.fits':
            tb_gbm_one = ReadGBMCatalogueInfo.select_one_by_name(tb_gbm, tb_one['GBM_assoc_key'])
            self.t05 = tb_one['T90_START']
            self.t90 = tb_one['T90']
            self.t25 = tb_gbm_one['T50_START']
            self.t50 = tb_gbm_one['T50']
        elif path_catalogue[-4:] == '.xml':
            if tb_one['GBM'] is None:
                self.t05 = 0
                self.t90 = 0
                self.t25 = 0
                self.t50 = 0
            else:
                self.t05 = tb_one['GBM']['T90_START']
                self.t90 = tb_one['GBM']['T90']
                self.t25 = tb_one['GBM']['T50_START']
                self.t50 = tb_one['GBM']['T50']
        self.met05 = self.met0 + self.t05
        self.t95 = self.t05 + self.t90
        self.met95 = self.met05 + self.t90
        self.met25 = self.met0 + self.t25
        logger.debug('Target name: {0}'.format(self.name))


    def synccutoff_time(self, t):
        if self.synccutoff is not None:
            return self.synccutoff*pow(max(t, self.t95)/1000., -3./8.)
        else:
            return None


##### Analysis Classes #####
PATH_BASEDIR = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone'

class AnalysisConfig:
    def __init__(self, target, emin, emax, tmin, tmax, evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., thetamax=65., index_fixed=None, suffix='', emin_fit=None, emax_fit=None, emin_eval=None, emax_eval=None, binned=False, psForce=False, roi_checked=None, gti_external=None):
        self.target = target
        logger.debug(self.target)
        self.binned = binned
        self.str_binned = 'Binned' if self.binned==True else 'Unbinned'
        self.emin = emin
        self.emax = emax
        self.emin_fit = emin_fit if emin_fit is not None else self.emin
        self.emax_fit = emax_fit if emax_fit is not None else self.emax
        self.emin_eval = emin_eval if emin_eval is not None else self.emin
        self.emax_eval = emax_eval if emax_eval is not None else self.emax
        logger.debug('Energy range: {emin} - {emax}'.format(emin=self.emin, emax=self.emax))
        self.tmin = tmin
        self.tmax = tmax
        self.metmin = (self.target.met0 + self.tmin) if self.tmin is not None else None
        self.metmax = (self.target.met0 + self.tmax) if self.tmax is not None else None
        self.evclass = evclass
        self.evtype = evtype
        self.ft2interval = ft2interval
        self.deg_roi = deg_roi
        self.zmax = zmax
        self.thetamax = thetamax
        self.index_fixed = index_fixed
        self.rad_margin = rad_margin
        self.suffix = suffix if suffix[0] in ('','_') else '_'+suffix
        self.roi_checked = self.deg_roi if roi_checked is None else roi_checked # Radius within which number of events is checked.

        self.str_time = 'T{0.tmin:0>6.0f}-{0.tmax:0>6.0f}'.format(self) if self.target.met0 >= 239557417 else 'MET{0.tmin:0>9.0f}-{0.tmax:0>9.0f}'.format(self) 
        self.str_energy = 'E{0.emin:0>7.0f}-{0.emax:0>7.0f}MeV'.format(self)
        self.str_roi = 'r{0:0>2.0f}deg'.format(self.deg_roi)
        self.str_index = 'Index{0}'.format(int(self.index_fixed*100) if (self.index_fixed==self.index_fixed and self.index_fixed is not None) else 'Free')
        
        self.dir_base = PATH_BASEDIR
        self.dir_work = '{base}/{target}/{energy}/{roi}/{time}/{spectype}/{binned}'.format(base=PATH_BASEDIR, target=self.target.name, energy=self.str_energy, roi=self.str_roi, time=self.str_time, spectype=self.target.spectraltype, binned=self.str_binned) #index=self.str_index)
        # Data
        self.path_dir_data = '{base}/{target}'.format(base=PATH_BASEDIR, target=self.target.name)

        # Initial values of products
        self.path_ft1 = None
        self.path_ft2 = None
        self.path_filtered = None
        self.path_filtered_gti = None
        self.path_ccube = None
        self.path_livetime = None
        self.path_exposure = None
        self.path_model_xml = None
        self.path_model_xml_new = None
        self.path_srcmap = None
        self.obs = None
        self.like = None
        self.likeobj = None
        self.loglike_inversed = None
        self.retcode = None

        # GTI
        #self.gti_filter = '(DATA_QUAL==1)&&(LAT_CONFIG==1)&&(ANGSEP(RA_ZENITH,DEC_ZENITH,{ra},{dec}) + {rad} < {zen})'.format(ra=self.target.ra, dec=self.target.dec, rad=12., zen=self.zmax)
        self.gti_filter = ' (DATA_QUAL==1)&&(LAT_CONFIG==1)'
        self.gti_external = gti_external
        self.roicut = False #True

        # CCube
        self.binsz = 0.2
        self.npix = int(self.deg_roi/self.binsz*2+0.5)

        # Livetime
        self.dcostheta = 0.05 #0.025
        self.binsz_lt = 1

        # Exposure
        self.irfs = 'CALDB'
        self.srcrad = self.deg_roi + self.rad_margin
        self.nlong_exp = int(self.srcrad * 2.+0.5)
        self.nlat_exp = int(self.srcrad * 2.+0.5)
        self.npix_exp = int(max(15+self.deg_roi, sqrt(2)*self.deg_roi)/self.binsz*2+0.5)
        self.nenergies_exp = int(log10(self.emax/self.emin)*10+0.5)
        self.nenergies_plot = int(log10(self.emax/self.emin)*4+0.5)
        self.energies = 10 ** np.linspace(log10(self.emin), log10(self.emax), self.nenergies_plot+1)
        logger.debug('Energy bins: {0}'.format(self.energies))

        # Modeling
        self.path_make3FGLxml = '/nfs/farm/g/glast/u/mtakahas/python_packages/make3FGLxml.py'
        self.radius_free_sources = 0.
        self.free_norms_only = True
        self.psForce = psForce

        # Catalogues
        self.path_catalogues = {'3FGL': '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/catalogProducts/v2r2/3FGL/gll_psc_v16.fit'}

        # Path of source files
        self.path_galdiff = '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/gll_iem_v06.fits'
        self.path_isodiff = '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/iso_{0}_v06.txt'.format(DCT_EVCLASSES[self.evclass])
        self.path_dir_extended = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/Catalogues/Extended_archive_v15/Templates'
        self.galdiff_name = os.path.basename(self.path_galdiff)[:-5]
        self.isodiff_name = os.path.basename(self.path_isodiff)[:-4]

        # Summary of results
        self.dct_summary_results = {}

        # CalOnly data
        self.calonly = None


    def set_directories(self):
        if not os.path.exists(self.dir_work):
            os.makedirs(self.dir_work)
            logger.info('Directory {0} is created.'.format(self.dir_work))
        os.chdir(self.dir_work)
        logger.info('Move to directory {0}.'.format(self.dir_work))


    def check_datafiles(self):
        path_ft1_exist = ls_list(self.path_dir_data+'/*_ft1*.fits')[0]
        path_ft2_exist = ls_list(self.path_dir_data+'/*_ft2-'+self.ft2interval+'.fits')[0]
        self.path_ft1 = path_ft1_exist if path_ft1_exist[0]=='/' else None
        logger.info('FT1 file: {0}'.format(self.path_ft1))
        self.path_ft2 = path_ft2_exist if path_ft2_exist[0]=='/' else None
        logger.info('FT2 file: {0}'.format(self.path_ft2))
        return (self.path_ft1, self.path_ft2)


    def filter(self, bforce=False):
        check_required_files({'FT1':self.path_ft1})
        self.path_filtered = '{0}/{1}'.format(self.dir_work, os.path.basename(self.path_ft1).replace('.fits', '_filtered.fits'))
        if check_required_files({'Filtered events':self.path_filtered})==0 and bforce==False:
            logger.info("""Filtering events is skipped.
""")
            return 0
        else:
            logger.info("""Filtering events is starting...""")

        my_apps.filter['evclass'] = self.evclass
        my_apps.filter['evtype'] = self.evtype
        my_apps.filter['ra'] = self.target.ra
        my_apps.filter['dec'] = self.target.dec
        if self.binned==False:
            my_apps.filter['rad'] = self.deg_roi
        else:
            my_apps.filter['rad'] = sqrt(2)*self.deg_roi
        my_apps.filter['emin'] = self.emin
        my_apps.filter['emax'] = self.emax
        my_apps.filter['zmax'] = self.zmax
        my_apps.filter['tmin'] = self.metmin
        my_apps.filter['tmax'] = self.metmax
        my_apps.filter['infile'] = self.path_ft1
        my_apps.filter['outfile'] = self.path_filtered
        my_apps.filter.run()
        logger.info("""Filtering events finished.
""")


    def maketime(self, bforce=False):
        check_required_files({'FT2':self.path_ft2, 'Filtered event':self.path_filtered})
        self.path_filtered_gti = self.path_filtered.replace('_filtered.fits', '_filtered_gti.fits')

        if self.gti_external is not None:
            path_gti_external = self.path_filtered_gti.replace('E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg'.format(emin=self.emin, emax=self.emax, roi=self.deg_roi), self.gti_external)
            logger.debug(self.gti_external)
            logger.debug(path_gti_external)

            tb_external = GetGTI.get_gti_table(path_gti_external)
            if len(tb_external)<1:
                self.duration = 0
                self.nevt_rough = 0
                return self.nevt_rough
            filter_external = '('
            for tstart, tstop in zip(tb_external['START'], tb_external['STOP']):
                filter_external += '(START>={sta}-0.5&&STOP<={sto}+0.5)||'.format(sta=tstart, sto=tstop)
            filter_external = filter_external[:-2] + ')'
            self.gti_filter = self.gti_filter + '&&' + filter_external
            logger.info('New GTI cut with external constraints: {0}'.format(self.gti_filter))

        if check_required_files({'Filtered GTI':self.path_filtered_gti})==0 and bforce==False:
            logger.info("""Making GTI is skipped.
""")
            return 0
        else:
            logger.info("""Making GTI is starting...""")

        tb_filtered = GetGTI.get_gti_table(self.path_filtered)
        logger.info('Filtered time: {0} intervals.'.format(len(tb_filtered)))

        my_apps.maketime['scfile'] = self.path_ft2
        my_apps.maketime['filter'] = self.gti_filter
        my_apps.maketime['roicut'] = 'yes' if self.roicut==True else 'no'
        my_apps.maketime['evfile'] = self.path_filtered
        my_apps.maketime['outfile'] = self.path_filtered_gti
        my_apps.maketime.run()
        logger.info("""Making GTI finished.
""")
        self.dct_summary_results['data_path'] = self.path_filtered_gti
        self.nevt_rough = get_entries_roi(self.path_filtered_gti, self.tmin, self.tmax, self.roi_checked, self.target.ra, self.target.dec, self.target.met0)
        if self.nevt_rough>0:
            logger.info('{0} events within {1} deg.'.format(self.nevt_rough, self.roi_checked))
        else:
            logger.warning('No valid events within {0} deg!!'.format(self.roi_checked))
        self.duration = GetGTI.get_duration(self.path_filtered_gti)['conservative']
        if self.duration<=0:
            logger.warning('GTI duration is zero!')
        return self.nevt_rough


    def evtbin(self, bforce=False):
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti})
        self.path_ccube = '{0}/{1}'.format(self.dir_work, os.path.basename(self.path_ft1).replace('.fits', '_ccube.fits'))
        my_apps.evtbin['evfile'] = self.path_filtered_gti
        my_apps.evtbin['scfile']= self.path_ft2
        my_apps.evtbin['algorithm'] = "ccube"
        my_apps.evtbin['outfile'] = self.path_ccube
        my_apps.evtbin['ebinalg'] = 'LOG'
        my_apps.evtbin['ebinfile'] = None
        my_apps.evtbin['emin'] = self.emin
        my_apps.evtbin['emax'] = self.emax
        my_apps.evtbin['enumbins'] = self.nenergies_exp
        my_apps.evtbin['coordsys'] = "CEL" 
        my_apps.evtbin['xref'] = self.target.ra
        my_apps.evtbin['yref'] = self.target.dec
        my_apps.evtbin['nxpix'] = self.npix
        my_apps.evtbin['nypix'] = self.npix
        my_apps.evtbin['binsz'] = self.binsz
        my_apps.evtbin['axisrot'] = 0.0 
        my_apps.evtbin['rafield'] = "RA" 
        my_apps.evtbin['decfield'] = "DEC" 
        my_apps.evtbin['proj'] = "CAR"
        logger.info('Going to make CCUBE.')
        my_apps.evtbin.run()
        logger.info("""Creating ccube finished.
""")


    def livetime(self, bforce=False):
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti})
        self.path_livetime = self.path_filtered_gti.replace('_filtered_gti.fits', '_ltCube.fits')
        if check_required_files({'Livetime':self.path_livetime})==0 and bforce==False:
            logger.info("""Calculating livetime is skipped.
""")
            return 0
        else:
            logger.info("""Calculating livetime is starting...""")
        my_apps.expCube['evfile'] = self.path_filtered_gti
        my_apps.expCube['scfile'] = self.path_ft2
        my_apps.expCube['outfile'] = self.path_livetime
        my_apps.expCube['zmax'] = self.zmax
        my_apps.expCube['dcostheta'] = self.dcostheta
        my_apps.expCube['binsz'] = self.binsz_lt
        my_apps.expCube.run()
        logger.info("""Calculating livetime finished.
""")


    def exposure(self, bforce=False):
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti, 'Livetime':self.path_livetime})
        self.path_exposure = self.path_livetime.replace('_ltCube.fits', '_expMap.fits')
        if check_required_files({'Exposure':self.path_exposure})==0 and bforce==False:
            logger.info("""Calculating exposure is skipped.
""")
            return 0
        else:
            logger.info("""Calculating exposure is starting...""")
        if self.binned==False:
            logger.info('Unbinned exposure map...')
            my_apps.expMap['evfile'] = self.path_filtered_gti
            my_apps.expMap['scfile'] = self.path_ft2
            my_apps.expMap['expcube'] = self.path_livetime
            my_apps.expMap['outfile'] = self.path_exposure
            my_apps.expMap['irfs'] = self.irfs
            my_apps.expMap['srcrad'] = self.srcrad
            my_apps.expMap['nlong'] = self.nlong_exp
            my_apps.expMap['nlat'] = self.nlat_exp
            my_apps.expMap['nenergies'] = self.nenergies_exp
            my_apps.expMap.run()
        else:
            logger.info('Binned exposure map...')
            #my_apps.gtexpcube2['evfile'] = self.path_filtered_gti
            #my_apps.gtexpcube2['scfile'] = self.path_ft2
            my_apps.gtexpcube2['infile'] = self.path_livetime
            my_apps.gtexpcube2['outfile'] = self.path_exposure
            my_apps.gtexpcube2['irfs'] = self.irfs
            my_apps.gtexpcube2['cmap'] = self.path_ccube
            my_apps.gtexpcube2['evtype'] = self.evtype
            #npix = int((10+self.deg_roi)/self.binsz*2+0.5)
            my_apps.gtexpcube2['nxpix'] = self.npix_exp
            my_apps.gtexpcube2['nypix'] = self.npix_exp
            my_apps.gtexpcube2['binsz'] = self.binsz 
            my_apps.gtexpcube2['coordsys'] = "CEL" 
            my_apps.gtexpcube2['xref'] = self.target.ra
            my_apps.gtexpcube2['yref'] = self.target.dec
            my_apps.gtexpcube2['axisrot'] = 0.0 
            my_apps.gtexpcube2['proj'] = "CAR" 
            my_apps.gtexpcube2['ebinalg'] = "LOG" 
            my_apps.gtexpcube2['emin'] = self.emin
            my_apps.gtexpcube2['emax'] = self.emax
            my_apps.gtexpcube2['enumbins'] = self.nenergies_exp
            my_apps.gtexpcube2.run()
        logger.info("""Calculating exposure finished.
""")


    def use_external_model(self, path_external):
        if path_external==None:
            logger.critical("""Variable {0} is NOT assigned!!!""".format(path_external))            
        if not os.path.exists(path_external):
            logger.critical("""File {0} does not exist!!!""".format(path_external))
            sys.exit(1)
        #self.path_model_xml = path_external
        self.path_model_xml = '/'.join((self.dir_work, os.path.basename(path_external)))
        shutil.copy(path_external, self.path_model_xml)
        logger.info('External XML file {0} is used for modeling.'.format(self.path_model_xml))


    def model_3FGL_sources(self, bforce=False, rm_catalogue_srcs=False):
        check_required_files({'Filtered GTI event':self.path_filtered_gti, 'Galactic diffuse':self.path_galdiff, 'Isotoropic diffuse':self.path_isodiff, 'Extened source templates':self.path_dir_extended})

        self.path_model_xml = self.path_filtered_gti.replace('_filtered_gti.fits', '_model.xml')
        if check_required_files({'Model':self.path_model_xml})==0 and bforce==False:
            logger.info("""Making 3FGL source model is skipped.
""")
            return 0
        else:
            logger.info("""Making 3FGL source model is starting...""")
        logger.debug('Output XML path: '+self.path_model_xml)
        mymodel = srcList(self.path_catalogues['3FGL'], self.path_filtered_gti, self.path_model_xml)

        logger.debug('Galdiff:'+ self.galdiff_name)
        logger.debug('Isodiff:'+ self.isodiff_name)
        if rm_catalogue_srcs==False:
            mymodel.makeModel(self.path_galdiff, self.galdiff_name, self.path_isodiff, self.isodiff_name, extDir=self.path_dir_extended, ExtraRad=self.rad_margin, radLim=self.radius_free_sources, normsOnly=self.free_norms_only, psForce=self.psForce)
        else:
            mymodel.makeModel(self.path_galdiff, self.galdiff_name, self.path_isodiff, self.isodiff_name, extDir=self.path_dir_extended, ExtraRad=0, radLim=0, normsOnly=self.free_norms_only, psForce=self.psForce)
        logger.info("""Making 3FGL source model finished.
""")


    def srcmaps(self, bforce=False):
        if self.binned==False:
            logger.error('srcmaps is NOT available for Binned analysis!!')
            return 1
        check_required_files({'FT2':self.path_ft2, 'CountCube':self.path_ccube, 'Livetime':self.path_livetime, 'SourceModel':self.path_model_xml, 'ExposureMap':self.path_exposure}) 
        my_apps.srcMaps['scfile'] = self.path_ft2
        my_apps.srcMaps['sctable'] = "SC_DATA" 
        my_apps.srcMaps['expcube'] = self.path_livetime
        my_apps.srcMaps['cmap'] = self.path_ccube
        my_apps.srcMaps['srcmdl'] = self.path_model_xml
        my_apps.srcMaps['bexpmap'] = self.path_exposure
        #my_apps.srcMaps['wmap'] =none 
        self.path_srcmap = self.path_filtered_gti.replace('_filtered_gti.fits', '_srcmap.fits')
        my_apps.srcMaps['outfile'] = self.path_srcmap
        my_apps.srcMaps['irfs'] = self.irfs #"P8R2_SOURCE_V6" 
        my_apps.srcMaps['evtype'] = self.evtype
        #my_apps.srcMaps['emapbnds'] = 'no'
        my_apps.srcMaps.run()
        logger.info("""Mapping sources finished.
""")


    def diffuse_responses(self, bforce=True, path_model_xml=None):
        if bforce==False:
            logger.warning('Calculating diffuse source responses is skipped.')
            return 0
        logger.info("""Calculating diffuse source responses is starting...""")
        if path_model_xml==None:
            path_model_xml = self.path_model_xml
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti, 'Source model':path_model_xml})

        my_apps.diffResps['evfile'] = self.path_filtered_gti
        my_apps.diffResps['scfile'] = self.path_ft2
        my_apps.diffResps['srcmdl'] = path_model_xml
        my_apps.diffResps['irfs'] = self.irfs
        my_apps.diffResps.run()
        logger.info("""Calculating diffuse source responses finished.
""")


    def setup(self, force={'download':False, 'filter':False, 'maketime':True, 'evtbin':False, 'livetime':False, 'exposure':False, 'model_3FGL_sources':False, 'diffuse_responses':True, 'srcmaps':False}, skip_zero_data=False):
        """Perform all preparation before fitting. Namely, downloading, filtering, making GTI, calculating livetime and exposure, modeling 3FGL sources, making diffuse source responses.
"""
        self.set_directories()
        self.download(bforce=force['download'])
        self.filter(bforce=force['filter'])
        self.maketime(bforce=force['maketime'])
        #logger.info('Roughly {0} events.'.format(self.nevt_rough))
        if self.binned==True:
            self.evtbin(bforce=force['evtbin'])
        if skip_zero_data==False:
            self.livetime(bforce=force['livetime'])
            self.exposure(bforce=force['exposure'])
            self.model_3FGL_sources(bforce=force['model_3FGL_sources'], rm_catalogue_srcs=False if self.nevt_rough>0 else True)
            self.diffuse_responses(bforce=force['diffuse_responses'])
            if self.binned==True:
                self.srcmaps(bforce=force['srcmaps'])
                if check_required_files({'srcmap':self.path_srcmap})>0:
                    logger.error('Making source map failed!!')
                    self.psForce = True
                    logger.warning('Source model is being recreated forcing diffuse sources to be point-like!')
                    self.model_3FGL_sources(bforce=force['model_3FGL_sources'], rm_catalogue_srcs=False if self.nevt_rough>0 else True)
                    logger.warning('Souce map is being recreated...')
                    self.srcmaps(bforce=force['srcmaps'])
                    if check_required_files({'srcmap':self.path_srcmap})>0:
                        logger.critical('Making source map failed again!!!')
                    else:
                        logger.info('Source model has been created forcing diffuse sources to be point-like.')

        elif skip_zero_data==True:
            if self.duration>0:
                self.livetime(bforce=force['livetime'])
                self.exposure(bforce=force['exposure'])
                self.model_3FGL_sources(bforce=force['model_3FGL_sources'], rm_catalogue_srcs=False if self.nevt_rough>0 else True)
                self.diffuse_responses(bforce=force['diffuse_responses'])
                if self.binned==True:
                    self.srcmaps(bforce=force['srcmaps'])
            else:
                logger.warning('Skipping calculation of livetime, exposre and diffuse responses.')
                #self.model_3FGL_sources(bforce=force['model_3FGL_sources'], rm_catalogue_srcs=False if self.nevt_rough>0 else True)
        return self.nevt_rough


    def set_parameter_initial(self, plike, parname, bound_lo, bound_hi, fixed=False):
        par = plike.getParam(parname)
        par.setScale(self.target.spectralpars_scale[parname])
        valset = self.target.spectralpars[parname]
        if self.target.spectraltype[-9:]=='ExpCutoff' and parname in ['Ebreak', 'P1'] and self.target.synccutoff in locals():
            if self.target.synccutoff is not None and self.target.synccutoff is not False and self.target.synccutoff==self.target.synccutoff:
                valset = self.target.synccutoff_time(self.tmin)
                logger.debug('Exponential cutoff: {0} MeV'.format(valset))
        par.setValue(valset)
        par.setBounds(bound_lo, bound_hi)
        plike.setParam(par)
        if fixed==True:
            plike.setParamAlwaysFixed(parname)
        logger.info('{n}: {v} {x}'.format(n=parname, v=valset, x='(fixed)' if fixed==True else ''))


    def set_parameters_initial(self, plike):
        for parname, parvalue in self.target.spectralpars.items():
            self.set_parameter_initial(plike, parname, self.target.spectralpars_bounds[parname][0], self.target.spectralpars_bounds[parname][1], fixed=True if parname in self.target.spectralpars_fixed else False)
            # par = plike.getParam(parname)
            # par.setScale(self.target.spectralpars_scale[parname])
            # par.setValue(parvalue)
            # par.setBounds(self.target.spectralpars_bounds[parname][0], self.target.spectralpars_bounds[parname][1])
            # plike.setParam(par)
            # if parname in self.target.spectralpars_fixed:
            #     plike.setParamAlwaysFixed(parname)


    def set_likelihood(self):
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti, 'Livetime':self.path_livetime, 'Exposure':self.path_exposure, 'Source model':self.path_model_xml})

        if self.binned==True:
            logger.info("""Setting up Binned likelihood...""")
            self.obs = BinnedObs(self.path_srcmap, self.path_livetime, self.path_exposure, irfs=self.irfs)
            self.like = BinnedAnalysis(self.obs, self.path_model_xml, optimizer='NewMinuit')
        else:
            logger.info("""Setting up Unbinned likelihood...""")
            self.obs = UnbinnedObs(self.path_filtered_gti, self.path_ft2, expMap=self.path_exposure, expCube=self.path_livetime, irfs=self.irfs)
            self.like = UnbinnedAnalysis(self.obs, self.path_model_xml, optimizer='NewMinuit')
        if self.binned==False:
            self.like.reset_ebounds(self.energies)
            logger.info('Energy bound has changed to {0}'.format(self.like.energies))

        # Target source
        if not self.target.name in self.like.sourceNames():
            logger.info("""Adding {0}.""".format(self.target.name))
            target_src = pyLike.PointSource(0, 0, self.obs.observation) #self.like.observation.observation)
            pl = pyLike.SourceFactory_funcFactory().create(self.target.spectraltype)
            self.set_parameters_initial(plike=pl)

            target_src.setSpectrum(pl)
            target_src.setName('{0}'.format(self.target.name))
            logger.debug('RA: {0}, DEC: {1}'.format(self.target.ra, self.target.dec))
            target_src.setDir(self.target.ra, self.target.dec,True,False)
            self.like.addSource(target_src)

        iso_norm_idx = self.like.par_index(self.isodiff_name, 'Normalization')
        self.like.freeze(iso_norm_idx)
        logger.info("""{0} has been fixed.""".format(self.isodiff_name))

        # Fix some parameters of the target
        p_idxs = [self.like.par_index(self.target.name, x) for x in self.target.spectralpars_fixed] 
        for p_idx in p_idxs:
            self.like.freeze(p_idx)
            logger.info(self.like.params()[p_idx])
        logger.info('Free parameter of {src}: {free}'.format(src=self.target.name, free=self.like.freePars(self.target.name)))

        self.like.writeXml()
        return self.like


    def set_likelihood_external_model(self, path_model):
        logger.info("""Setting up likelihood by {0} ...""".format(path_model))
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti, 'Livetime':self.path_livetime, 'Exposure':self.path_exposure, 'Source model':path_model})
        self.diffuse_responses(bforce=True, path_model_xml=path_model)
        self.obs = UnbinnedObs(self.path_filtered_gti, self.path_ft2, expMap=self.path_exposure, expCube=self.path_livetime, irfs=self.irfs)
        self.like = UnbinnedAnalysis(self.obs, path_model, optimizer='NewMinuit')
        if self.target.synccutoff is not None and self.target.synccutoff is not False and self.target.synccutoff==self.target.synccutoff:
            self.set_target_synccutoff()

        self.path_model_xml = path_model
        self.path_model_xml_new = None
        self.likeobj = pyLike.NewMinuit(self.like.logLike)
        self.loglike_inversed = self.like()
        self.retcode = None
        if self.binned==False:
            self.like.reset_ebounds(self.energies)
            logger.info('Energy bound has changed to {0}'.format(self.like.energies))


    def fit(self, bredo=True):
        self.set_likelihood()
        self.likeobj = pyLike.NewMinuit(self.like.logLike)
        #likeobj = pyLike.NewMinuit(like.logLike)
        logger.info("""Likelihood fitting is starting...""")
        logger.debug('Normalization of {0}: {1}'.format(self.target.name, self.like.normPar(self.target.name).getValue()))
        if bredo==True:
            try:
                self.dofit(tol=self.like.tol)
                logger.debug('Normalization of {0}: {1}'.format(self.target.name, self.like.normPar(self.target.name).getValue()))
            except RuntimeError:
                logger.error('RuntimeError!!')
                logger.info('Normalization value: {0}'.format(self.like.normPar(self.target.name).getValue()))
                if self.like.normPar(self.target.name).getValue() != self.like.normPar(self.target.name).getValue():
                    self.remove_other_sources()
                    return 1

                logger.warning('Tolerance is relaxed from {tol0} to {tol1}'.format(tol0=self.like.tol, tol1=self.like.tol*10.))
                logger.info('Resetting likelihood.')
                self.set_likelihood()
                self.likeobj = pyLike.NewMinuit(self.like.logLike)
                logger.debug('Normalization of {0}: {1}'.format(self.target.name, self.like.normPar(self.target.name).getValue()))
                logger.info('Fixing normalization of other sources.')
                for source in self.like.sourceNames():
                    if source not in (self.target.name):
                        self.like.normPar(source).setFree(False)
                logger.debug('Normalization of {0}: {1}'.format(self.target.name, self.like.normPar(self.target.name).getValue()))
                logger.info('Fitting again...')
                #self.like.normPar(self.target.name).setTrueValue(1E-12)
                try:
                    self.dofit(tol=self.like.tol*10.)
                except RuntimeError:
                    self.remove_other_sources()
                    return 1
        else:
            self.dofit()

        if self.retcode>0 or bredo==True:
            logger.info("""Redoing fitting...""")
            sourceDetails = {}
            for source in self.like.sourceNames():
                sourceDetails[source] = self.like.Ts(source)
            logger.info('Deleting non-significant sources...')
            for source,ts in sourceDetails.iteritems():
                logger.info('{0} TS={1}'.format(source, ts))
                if (ts < 1 and source not in (self.target.name, self.isodiff_name, self.galdiff_name)):
                    logger.info("Deleting... ")
                    self.like.deleteSource(source)

            # Fix Ebreak
            # if self.target.spectraltype=='ExpCutoff':
            #     eb_idx = self.like.par_index(self.target.name, 'Ebreak')
            #     self.like.freeze(eb_idx)
            #     logger.info(self.like.params()[eb_idx])
            #     logger.info('Free parameter of {src}: {free}'.format(src=self.target.name, free=self.like.freePars(self.target.name)))

            try:
                self.dofit()
            except RuntimeError:
                logger.error('RuntimeError!!')
                logger.warning('Tolerance is relaxed from {tol0} to {tol1}'.format(tol0=self.like.tol, tol1=self.like.tol*10.))
                try: 
                    self.dofit(tol=self.like.tol*10.)
                except RuntimeError:
                    logger.error('RuntimeError!!')
                    logger.warning('Fixing normalization of other sources with TS < 4.')
                    sourceDetails = {}
                    for source in self.like.sourceNames():
                        sourceDetails[source] = self.like.Ts(source)
                    for source,ts in sourceDetails.iteritems():
                        logger.info('{0} TS={1}'.format(source, ts))
                        if source != self.target.name and ts<4.0:
                            self.like.normPar(source).setFree(False)
                    self.dofit(tol=self.like.tol*10.)
            #self.like.setEnergyRange(self.emin, self.emax)

        # Save new XML model file
        self.path_model_xml_new = self.path_model_xml.replace('.xml', '_new.xml')
        self.like.writeXml(self.path_model_xml_new)
        return (self.retcode, self.loglike_inversed)


    def like_calonly(self):
        #logger.info('CalOnly likelihood')

        nobs = self.calonly.onevt
        #ecalc_min, ecalc_max = log10(self.emin), log10(self.emax)
        #ecalc = 10 ** np.linspace(ecalc_min, ecalc_max, int((ecalc_max-ecalc_min)*20.+1.5))
        npred = np.zeros_like(nobs)
        #for ie,(line0,line1) in enumerate(zip(ecalc[:-1], ecalc[1:])):
        for ie,(loge0,loge1) in enumerate(zip(self.calonly.energies[:-1], self.calonly.energies[1:])):
            # Local gamma
            exposure = self.calonly.func_onexp((loge0+loge1)/2.) #self.calonly.onexp[ie] * 1e4
            #exposure = self.calonly.func_onexp((np.log10(line0)+np.log10(line1))/2.) #self.calonly.onexp[ie] * 1e4
            flux_all = self.like.flux(self.target.name, 10**loge0, 10**loge1) #0.
            npred[ie] = flux_all * 0.68 * exposure
            # if ie==0:
            #     logger.info('  Exposure: {0}'.format(exposure))
            #     logger.info('  Flux: {0}'.format(flux_all))
        npred += self.calonly.bkgevt
        npred_total = sum(npred)

        calonlylike = -np.log(np.exp(-npred_total) * np.prod(pow(npred, nobs)/scipy.misc.factorial(nobs)))
        return calonlylike


    def like_with_calonly(self):
        calonlylike = self.like_calonly()
        # logger.info('Nobs: {0}'.format(nobs))
        # logger.info('Npred: {0}'.format(npred))
        # logger.info('CalOnly log-likelihood inversed: {0}'.format(calonlylike))
        return self.like() + calonlylike


    def fit_fixed_index(self, index, path_xml_temp):
        index_idx = self.like.par_index(self.target.name, 'Index')
        self.like[index_idx] = index
        self.like.freeze(index_idx)
        try:
            self.dofit()
        except RuntimeError:
            logger.error('RuntimeError!!')
            logger.warning('Tolerance is relaxed from {tol0} to {tol1}'.format(tol0=self.like.tol, tol1=self.like.tol*10.))
            self.dofit(tol=self.like.tol*10.)

        # Save new XML model file
        self.like.writeXml(path_xml_temp)
        return (self.retcode, self.loglike_inversed)

            
    def dofit(self, tol=None):
        self.loglike_inversed = self.like.fit(tol=tol, verbosity=0,covar=True,optObject=self.likeobj)
        logger.info('Likelihood fit results:')
        logger.info('-Log(likelihood) = {0}'.format(self.loglike_inversed))
        logger.info(self.like.model[self.target.name])
        self.retcode = self.likeobj.getRetCode()
        if self.retcode>0:
            logger.warning("""Return code is {0}!!
""".format(self.retcode))
        else:
            logger.info("""Fitting successfully done.
""")


    def remove_other_sources(self):
        for source in self.like.sourceNames():
            if not source in (self.target.name, self.galdiff_name, self.isodiff_name):
                logger.info("Deleting {0}... ".format(source))
                self.like.deleteSource(source)
        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        self.like[norm_idx] = self.target.spectralpars[self.target.norm_name]
        self.path_model_xml_new = self.path_model_xml.replace('.xml', '_new.xml')
        self.like.writeXml(self.path_model_xml_new)
        logger.info('Simpified source model is saved as {0}'.format(self.path_model_xml_new))


    def reset_target_norm(self):
        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        self.like[norm_idx] = self.target.spectralpars[self.target.norm_name]


    def set_target_synccutoff(self):
        vcutoff = self.target.synccutoff_time(self.tmin)
        ebreak_idx = self.like.par_index(self.target.name, 'Ebreak')
        self.like[ebreak_idx] = vcutoff
        self.like.freeze(ebreak_idx)
        p1_idx = self.like.par_index(self.target.name, 'P1')
        self.like[p1_idx] = vcutoff
        self.like.freeze(p1_idx)
        logger.debug('Exponential cutoff: {0} MeV'.format(vcutoff))                


    def summarize_fit_results(self, use_calonly=False):
        """Summarize the results after fitting and return a dictonary. The contents are the model parameters, photon flux, TS, return code.
"""
        for name_param in self.target.spectralpars.keys(): #model_pars:
            param = self.like.model[self.target.name].funcs['Spectrum'].getParam(name_param)
            self.dct_summary_results[name_param] = {'value':param.value(), 'error':param.error()}

        # Current best loglike
        self.dct_summary_results['loglike_inversed'] = self.like()
        if use_calonly==True:
            self.dct_summary_results['loglike_calonly_inversed'] = self.like_calonly()
            self.dct_summary_results['loglike_with_calonly_inversed'] = self.dct_summary_results['loglike_inversed'] + self.dct_summary_results['loglike_calonly_inversed']

        # TS
        name = self.target.name
        logger.debug('TS of {0}:'.format(name))
        self.dct_summary_results['TS'] = self.like.Ts(str(name))

        # Flux
        if self.dct_summary_results['TS']>0:
            flux_and_err = self.eval_flux_and_error(self.target.name, self.emin_eval, self.emax_eval)

        # Return code 
        self.dct_summary_results['retcode'] = self.retcode

        logger.debug(self.dct_summary_results)
        return self.dct_summary_results


    def eval_flux_and_error(self, name=None, emin=None, emax=None):
        name_eval = name if name is not None else self.target.name
        e0 = emin if emin is not None else self.emin_eval
        e1 = emax if emax is not None else self.emax_eval
        flux = self.like.flux(name_eval, e0, e1)
        flux_err = self.like.fluxError(name_eval, e0, e1)
        self.dct_summary_results['flux'] = {'value':flux, 'error':flux_err}
        eflux = self.like.energyFlux(name_eval, e0, e1)
        eflux_err = self.like.energyFluxError(name_eval, e0, e1)
        self.dct_summary_results['eflux'] = {'value':eflux, 'error':eflux_err}
        return (flux, flux_err)
        

    def eval_flux_and_error_total(self, emin=None, emax=None):
        e0 = emin if emin is not None else self.emin_eval
        e1 = emax if emax is not None else self.emax_eval
        flux_total = 0.
        flux_err_total_sq = 0.
        for source in self.like.sourceNames():
            fsrc = self.eval_flux_and_error(source, e0, e1)
            flux_total += fsrc[0]
            flux_err_total_sq += pow(fsrc[1], 2)
        flux_err_total = sqrt(flux_err_total_sq)
        self.dct_summary_results['flux_total'] = {'value':flux_total, 'error':flux_err_total}
        return (flux_total, flux_err_total)


    def eval_limits_powerlaw(self, emin=None, emax=None, eref=None, str_index_fixed=['best', 'free']):
        e0 = emin if emin is not None else self.emin_eval
        e1 = emax if emax is not None else self.emax_eval
        if eref is not None:
            e2 = eref
        elif self.target.spectraltype[-8:]=='PowerLaw':
            e2 = self.target.spectralpars['Scale']
        else:
            logger.critical('eref is {0}!!!'.format(eref))
            sys.exit(1)

        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        norm_value = self.dct_summary_results[self.target.norm_name]['value'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).value()
        norm_error = self.dct_summary_results[self.target.norm_name]['error'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).error()
        logx_lowest = -4.0
        logx_highest = max(4.0, 1+np.log10(norm_error/norm_value))
        nx = min(100, 10 * (logx_highest-logx_lowest))
        xvals = max(norm_value, norm_error) * 10 ** np.linspace(logx_lowest, logx_highest, nx)
        logger.info('Normarization = {0} +/- {1}'.format(norm_value, norm_error))
        if np.inf in xvals:
            logger.warning('Infinite profile normalization value exists!!')
            xvals = 10 ** np.linspace(-20, 0, 100)
        #xvals = np.insert(xvals, 0, 0.0)
        if norm_value<min(xvals):
            xvals_inter = 10 ** np.linspace(np.log(norm_value), np.log10(min(xvals)), (np.log10(min(xvals))-np.log(norm_value))*5.+1 )
            xvals_infra = 10 ** np.linspace(np.log(norm_value)-2, np.log(norm_value), 11)
            xvals = np.insert(xvals, 0, xvals_inter[:-1])
            xvals = np.insert(xvals, 0, xvals_infra[:-1])
        if norm_value>max(xvals):
            xvals = np.insert(xvals, len(xvals), norm_value)
            xvals = np.insert(xvals, len(xvals), norm_value*10)
        self.like.normPar(self.target.name).setBounds(0, max(1e-2, xvals[-1])) #xvals[0], xvals[-1])
        logger.info("""Profile normalization factor: 
{0}
{1} points.""".format(xvals, len(xvals)))

        # Index parameter
        index_name = 'Index'
        index_idx = self.like.par_index(self.target.name, index_name)
        index_value = self.dct_summary_results[index_name]['value'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).value()
        index_error = self.dct_summary_results[index_name]['error'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).error()
        index_values = {'best':index_value, 'free':index_value}
        index_values['harder'] = index_value + index_error * (1 if index_values<0 else -1)
        index_values['softer'] = index_value + index_error * (-1 if index_values<0 else 1)

        # Current best loglike
        loglike0 = -self.like()
        # Current values
        v0 = {}
        v0['flux'] = self.like[self.target.name].flux(e0, e1)
        v0['eflux'] = self.like[self.target.name].energyFlux(e0, e1)

        # Profile results
        o = {}
        # Limit results
        limits = {}
        # Profile normalization
        failed_freeIndex = False
        for str_index_assumed in str_index_fixed:
            logger.info('Index = {idxa} ({st}) is assumed.'.format(idxa=index_values[str_index_assumed], st=str_index_assumed))

            cl_1sigma = 0.6827 #1.
            cl_2sigma = 0.9545#4.
            if str_index_assumed=='free':
                cl_1sigma = 0.87063
                cl_2sigma = 0.96821
                
            if self.target.spectraltype[-8:]=='PowerLaw' or self.target.spectraltype[-9:]=='PowerLaw2':
                if self.target.spectraltype[-8:]=='PowerLaw':
                    v0['dnde'] = norm_value * (e2 / self.target.spectralpars['Scale']) ** index_values[str_index_assumed]
                elif self.target.spectraltype[-9:]=='PowerLaw2':
                    v0['dnde'] =  norm_value * (index_values[str_index_assumed]+1.) * pow(e2, index_values[str_index_assumed]) / (pow(self.target.spectralpars['UpperLimit'], index_values[str_index_assumed]+1) - pow(self.target.spectralpars['LowerLimit'], index_values[str_index_assumed]+1))
                v0['e2dnde'] = v0['dnde'] * e2 * e2
                if not index_values[str_index_assumed]==index_values[str_index_assumed]:
                    logger.error('Index value is NOT valid!!! {0}'.format(index_values[str_index_assumed]))
                    sys.exit(1)
                self.like[index_idx] = index_values[str_index_assumed]
                if str_index_assumed == 'free':
                    self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[index_idx:index_idx+1], value=1)
                else:
                    self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[index_idx:index_idx+1], value=0)
            o[str_index_assumed] = {'xvals': xvals,
                                    'dloglike': np.zeros(len(xvals)),
                                    'loglike': np.zeros(len(xvals))
                                    }
            limits[str_index_assumed] = {}
            # Loop for profiling of normalization
            for i, x in enumerate(xvals):
                self.like[norm_idx] = x
                self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[norm_idx:norm_idx+1], value=0)
                if str_index_assumed == 'free':
                    try:
                        loglike1 = -self.like.fit(verbosity=0,covar=False,optObject=self.likeobj)
                    except RuntimeError:
                        logger.error('RuntimeError!!!')
                        logger.error('Fitting with free index failed!!')
                        failed_freeIndex = True
                        #sys.exit(1)
                else:
                    loglike1 = -self.like() 
                logger.debug('No.{ip} normalization factor = {xval:.3E}, expected count = {nph}, loglikelihood = {ll:.3E}'.format(ip=i, xval=x, nph=self.like._npredValues(), ll=loglike1)) #NpredValue(self.target.name)))
                o[str_index_assumed]['dloglike'][i] = loglike1 - loglike0
                o[str_index_assumed]['loglike'][i] = loglike1

            # limits
            if str_index_assumed=='free' and failed_freeIndex==True:
                break
            limits[str_index_assumed]['norm'] = get_parameter_limits(xval=xvals, loglike=o[str_index_assumed]['dloglike'], cl_limit=cl_2sigma, cl_err=cl_1sigma)
            #if self.target.spectraltype=='PowerLaw':
            for item in ('flux', 'eflux', 'dnde', 'e2dnde'):
                limits[str_index_assumed][item] = {}
                limits[str_index_assumed][item]['x0'] = v0[item] * limits[str_index_assumed]['norm']['x0'] / norm_value
                limits[str_index_assumed][item]['err_lo'] = v0[item] * limits[str_index_assumed]['norm']['err_lo'] / norm_value
                limits[str_index_assumed][item]['err_hi'] = v0[item] * limits[str_index_assumed]['norm']['err_hi'] / norm_value
                limits[str_index_assumed][item]['ul'] = v0[item] * limits[str_index_assumed]['norm']['ul'] / norm_value
                limits[str_index_assumed][item]['ll'] = v0[item] * limits[str_index_assumed]['norm']['ll'] / norm_value
                limits[str_index_assumed][item]['err'] = v0[item] * limits[str_index_assumed]['norm']['err'] / norm_value

        self.dct_summary_results['limits'] = limits
        return limits


    def eval_limits_powerlaw_index(self, emin=None, emax=None, eref=None):
        e0 = emin if emin is not None else self.emin_eval
        e1 = emax if emax is not None else self.emax_eval
        e2 = eref if (eref is not None or self.target.spectraltype!='PowerLaw') else self.target.spectralpars['Scale']

        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        norm_value = self.dct_summary_results[self.target.norm_name]['value'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).value()
        norm_error = self.dct_summary_results[self.target.norm_name]['error'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).error()
        norm_values = {'best':norm_value, 'free':norm_value}

        # Index parameter
        index_name = 'Index'
        index_idx = self.like.par_index(self.target.name, index_name)
        index_value = self.dct_summary_results[index_name]['value'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).value()
        index_error = self.dct_summary_results[index_name]['error'] #self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).error()
        logx_lowest = -4.0
        logx_highest = max(4.0, 1+np.log10(norm_error/norm_value))

        xll = max(min(index_value - 2.*index_error, -4), self.target.spectralpars_bounds[index_name][0])
        xul = min(max(index_value + 2.*index_error, -1), self.target.spectralpars_bounds[index_name][1])
        nx = min(300, int((xul-xll)*100.))
        xvals = np.linspace(xll, xul, nx+1)
        logger.info('Spectral index = {0} +/- {1}'.format(index_value, index_error))
        #self.like.normPar(self.target.name).setBounds(xvals[0], xvals[-1])
        logger.info("""Profile spectral index: 
{0}
{1} points.""".format(xvals, len(xvals)))

        # Current best loglike
        loglike0 = -self.like()

        # Profile results
        o = {}
        # Limit results
        limits = {}
        # Profile normalization
        failed_freeIndex = False
        for str_norm_assumed in ('free',):
            logger.info('Normalization = {na} ({st}) is assumed.'.format(na=norm_values[str_norm_assumed], st=str_norm_assumed))
            limits[str_norm_assumed] = {}

            cl_1sigma = 0.6827 #1.
            cl_2sigma = 0.9545#4.
            if str_norm_assumed=='free':
                cl_1sigma = 0.87063
                cl_2sigma = 0.96821
                
            #if self.target.spectraltype[-8:]=='PowerLaw' or self.target.spectraltype[-9:]=='PowerLaw2':
            self.like[norm_idx] = norm_value
            if str_norm_assumed == 'free':
                self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[norm_idx:norm_idx+1], value=1)
            else:
                self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[norm_idx:norm_idx+1], value=0)
            o[str_norm_assumed] = {'xvals': xvals,
                                   'dloglike': np.zeros(len(xvals)),
                                   'loglike': np.zeros(len(xvals))
                                   }

            # Loop for profiling of normalization
            for i, x in enumerate(xvals):
                self.like[index_idx] = x
                self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[index_idx:index_idx+1], value=0)
                if str_norm_assumed == 'free':
                    try:
                        loglike1 = -self.like.fit(verbosity=0,covar=False,optObject=self.likeobj)
                    except RuntimeError:
                        logger.error('RuntimeError!!!')
                        logger.error('Fitting with free normalization failed!!')
                        failed_freeIndex = True
                        #sys.exit(1)
                else:
                    loglike1 = -self.like()
                logger.debug('No.{ip} spectral index = {xval:.3E}, loglikelihood = {ll:.3E}'.format(ip=i, xval=x, ll=loglike1)) 
                #logger.debug('No.{ip} spectral index = {xval:.3E}, expected count = {nph}, loglikelihood = {ll:.3E}'.format(ip=i, xval=x, nph=self.like._npredValues(), ll=loglike1)) 
                o[str_norm_assumed]['dloglike'][i] = loglike1 - loglike0
                o[str_norm_assumed]['loglike'][i] = loglike1
            limits[str_norm_assumed]['index'] = get_parameter_limits(xval=xvals, loglike=o[str_norm_assumed]['dloglike'], cl_limit=cl_2sigma, cl_err=cl_1sigma)

        self.dct_summary_results['index_limit'] = limits
        return limits


    def map_norm_range(self, ndivperdec=10, insert_zero=False): #Let ndivperdec be back to 20
        norm_value = self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).value()
        norm_error = self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).error()
        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        if norm_error>0:
            logx_lowest = -3.0
            logx_highest = max(2.0, 4.*np.log10(norm_error/norm_value))
            nx = min(300, ndivperdec * (logx_highest-logx_lowest))
            xvals = norm_value * 10 ** np.linspace(logx_lowest, logx_highest, nx)
        else:
            logx_lowest = -8.0
            logx_highest = -2.0
            nx = min(400, ndivperdec * (logx_highest-logx_lowest))
            xvals = 10 ** np.linspace(logx_lowest, logx_highest, nx)

        if np.inf in xvals:
            logger.warning('Infinite profile normalization value exists!!')
            logx_lowest = -7.0
            logx_highest = -3.0
            nx = min(1000, ndivperdec * (logx_highest-logx_lowest))
            xvals = 10 ** np.linspace(logx_lowest, logx_highest, nx)
        # if norm_value<min(xvals):
        #     xvals_inter = 10 ** np.linspace(np.log(norm_value), np.log10(min(xvals)), (np.log10(min(xvals))-np.log(norm_value))*5.+1 )
        #     xvals_infra = 10 ** np.linspace(np.log(norm_value)-2, np.log(norm_value), 11)
        if insert_zero==True:
            xvals = np.insert(xvals, 0, 0.)
            #xvals = np.insert(xvals, 0, xvals_inter[:-1])
        #     xvals = np.insert(xvals, 0, xvals_infra[:-1])
        # if norm_value>max(xvals):
        #     xvals = np.insert(xvals, len(xvals), norm_value)
        #     xvals = np.insert(xvals, len(xvals), norm_value*10)
        return xvals


    def map_index_range(self, ndivperuni=20, index_range=(-6., 2.)): # Back to 40
        # Index parameter
        # index_name = 'Index'
        # index_idx = self.like.par_index(self.target.name, index_name)
        # index_value = self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).value()
        # index_error = self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).error()
        #xvals = np.linspace(index_value-3.*index_error, index_value+3.*index_error, ndiv+1)

        xvals_mid = np.linspace(-3, 0, 3.*ndivperuni*2.+1)
        xvals_soft = np.linspace(index_range[0], -3., (-3.-index_range[0])*ndivperuni+1)
        xvals_hard = np.linspace(0., index_range[1], (index_range[1]-0.)*ndivperuni+1)
        xvals = np.hstack((xvals_soft, xvals_mid[1:-1], xvals_hard)) #np.linspace(index_range[0], index_range[1], ndiv+1)
        return xvals


    def scan_norm_and_index(self, eref, use_calonly=False):
        logger.info('Parameter scan over normalization and spectral index.')
        if use_calonly==True:
            if self.calonly is None:
                logger.error('CalOnly data have NOT been added!!')
                use_calonly = False
            else:
                logger.info('CalOnly data are used...')
        else:
                logger.info('CalOnly data are not used...')
        # Mapping scanned values
        norms = self.map_norm_range(insert_zero=True)
        indices = self.map_index_range() #120) #, index_range=(-3, 3))
        self.like.normPar(self.target.name).setBounds(0, max(1e-2, norms[-1]))
        logger.info("""Profile normalization factor: 
{0}
{1} points.""".format(norms, len(norms)))
        indices_mesh, norms_mesh = np.meshgrid(indices, norms)
        loglikes = np.zeros_like(indices_mesh)
        efluxes = np.zeros_like(indices_mesh)
        fluxes = np.zeros_like(indices_mesh)
        efluences = np.zeros_like(indices_mesh)
        e2dnde = np.zeros_like(indices_mesh)
        npreds_caltkr = np.zeros_like(indices_mesh)
        npreds_calonly = np.zeros_like(indices_mesh)
        if use_calonly==True:
            loglikes_calonly = np.zeros_like(indices_mesh)
            loglikes_without_calonly = np.zeros_like(indices_mesh)

        # Indexes in likelihood object
        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        index_name = 'Index'
        index_idx = self.like.par_index(self.target.name, index_name)        

        loge_segmented = np.linspace(np.log10(self.like.energies[0]), np.log10(self.like.energies[-1]), int((-np.log10(self.like.energies[0])+np.log10(self.like.energies[-1]))*20+1.5))
        line_segmented = 10 ** loge_segmented

        for inorm, jindex in itertools.product(range(len(norms)), range(len(indices))):
            # Setting model
            logger.debug('Normalization: No.{0} {1}'.format(inorm, norms_mesh[inorm][jindex]))
            logger.debug('Spectral index: No.{0} {1}'.format(jindex, indices_mesh[inorm][jindex]))
            self.like[norm_idx] = norms[inorm]
            self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[norm_idx:norm_idx+1], value=0)
            self.like[index_idx] = indices[jindex]
            self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[index_idx:index_idx+1], value=0)
            # Npred(CalTkr)
            sum_model = np.zeros_like(self.like._srcCnts(self.like.sourceNames()[0]))
            for srcname in self.like.sourceNames():
                sum_model = sum_model + self.like._srcCnts(srcname)
            npreds_caltkr[inorm][jindex] = sum(sum_model)
            #Npred(CalOnly)
            if use_calonly==True:
                for je in range(len(line_segmented)-1):
                    npreds_calonly[inorm][jindex] += self.like.flux(self.target.name, line_segmented[je], line_segmented[je+1]) * self.calonly.func_onexp((loge_segmented[je]+loge_segmented[je+1])/2.) * 0.68
                npreds_calonly[inorm][jindex] += sum(self.calonly.bkgevt)
            # Evaluation
            loglikes[inorm][jindex] = -self.like() #if use_calonly==False else -self.like_with_calonly()
            if use_calonly==True:
                loglikes_without_calonly[inorm][jindex] = loglikes[inorm][jindex]
                loglikes_calonly[inorm][jindex] = -self.like_calonly()
                loglikes[inorm][jindex] += loglikes_calonly[inorm][jindex]
            fluxes[inorm][jindex] = self.like.flux(self.target.name, self.emin_eval, self.emax_eval)
            efluxes[inorm][jindex] = self.like.energyFlux(self.target.name, self.emin_eval, self.emax_eval)
            efluences[inorm][jindex] = efluxes[inorm][jindex] * self.duration
            if self.target.spectraltype[-8:]=='PowerLaw':
                e2dnde[inorm][jindex] = norms[inorm] * eref * eref * pow(eref/self.target.spectralpars['Scale'], indices[jindex])
            elif self.target.spectraltype[-9:]=='PowerLaw2':
                if indices[jindex]!=-1:
                    e2dnde[inorm][jindex] = norms[inorm] * (indices[jindex]+1.) * pow(eref, indices[jindex]) / (pow(self.target.spectralpars['UpperLimit'], indices[jindex]+1) - pow(self.target.spectralpars['LowerLimit'], indices[jindex]+1)) * eref * eref
                else:
                    e2dnde[inorm][jindex] = norms[inorm] / (np.log(self.target.spectralpars['UpperLimit']) - np.log(self.target.spectralpars['LowerLimit'])) * eref
            else:
                #logger.warning('Use PowerLaw or PowerLaw2 instead of {0}!!'.format(self.target.spectraltype))
                e2dnde = None
                #sys.exit(1)

        logger.info("""likelihood:
{0}""".format(loglikes))
        loglike_inv_best = -np.max(loglikes, axis=None) #min((self.dct_summary_results['loglike_inversed'] if use_calonly==False else self.dct_summary_results['loglike_calonly_inversed']), -np.max(loglikes, axis=None))
        logger.info("""Best likelihood: {0}""".format(loglike_inv_best))
        dloglike = -loglikes - loglike_inv_best 
        logger.info("""d-likelihood:
{0}""".format(dloglike))
        dloglike_doubled = 2.*dloglike
        # CalOnly
        if use_calonly==True:
            loglike_calonly_inv_best = -np.max(loglikes_calonly, axis=None) 
            loglike_without_calonly_inv_best = -np.max(loglikes_without_calonly, axis=None) 
            dloglike_calonly = -loglikes_calonly - loglike_calonly_inv_best 
            dloglike_without_calonly = -loglikes_without_calonly - loglike_without_calonly_inv_best 
            dloglike_calonly_doubled = 2.*dloglike_calonly
            dloglike_without_calonly_doubled = 2.*dloglike_without_calonly

        # Result storage
        self.dct_summary_results['dloglike'] = {}
        self.dct_summary_results['dloglike']['normalization'] = norms
        self.dct_summary_results['dloglike']['Index'] = indices
        self.dct_summary_results['dloglike']['loglike'] = loglikes
        self.dct_summary_results['dloglike']['npred'] = {'Regular':npreds_caltkr, 'CalOnly':npreds_calonly}
        self.dct_summary_results['dloglike']['dloglike'] = dloglike
        if use_calonly==True:
            self.dct_summary_results['dloglike']['loglike_calonly'] = loglikes_calonly
            self.dct_summary_results['dloglike']['loglike_without_calonly'] = loglikes_without_calonly
            self.dct_summary_results['dloglike']['dloglike_calonly'] = dloglike_calonly
            self.dct_summary_results['dloglike']['dloglike_without_calonly'] = dloglike_without_calonly
        self.dct_summary_results['dloglike']['flux'] = fluxes
        self.dct_summary_results['dloglike']['eflux'] = efluxes
        self.dct_summary_results['dloglike']['efluence'] = efluences
        self.dct_summary_results['dloglike']['e2dnde'] = e2dnde
        arg_bestlike = dloglike.argmin()
        args_bestlike = zip(np.where(dloglike==dloglike.min())[0], np.where(dloglike==dloglike.min())[1])[0]
        logger.info('Best likelihood point: {0}'.format(args_bestlike))
        self.dct_summary_results['dloglike']['best'] = {}
        self.dct_summary_results['dloglike']['best']['normalization'] = norms[args_bestlike[0]]
        logger.info('Normalization (best likelihood): {0}'.format(self.dct_summary_results['dloglike']['best']['normalization']))
        self.dct_summary_results['dloglike']['best']['Index'] = indices[args_bestlike[1]]
        logger.info('Index (best likelihood): {0}'.format(self.dct_summary_results['dloglike']['best']['Index']))
        for k,v in self.dct_summary_results['dloglike'].items():
            if k in ('loglike', 'dloglike', 'flux', 'eflux', 'efluence', 'e2dnde'):
                if v is not None:
                    self.dct_summary_results['dloglike']['best'][k] = v.flatten()[arg_bestlike]
                else:
                    self.dct_summary_results['dloglike']['best'][k] = None
        self.dct_summary_results['dloglike']['TS'] =2. * (self.dct_summary_results['dloglike']['loglike'][0,0] - self.dct_summary_results['dloglike']['best']['loglike'])
        logger.info('TS of {0}: {1}'.format(self.target.name, self.dct_summary_results['dloglike']['TS']))
        #logger.info('CalOnly best d-likelihood: {0}'.format(self.dct_summary_results['dloglike']['dloglike'][218,120]))

        self.plot_spectrum_scanned2D(name='dloglike', norms_mesh=norms_mesh, e2dnde=e2dnde, efluences=efluences, indices_mesh=indices_mesh, zvalues=dloglike_doubled, cont_levels=[2.30, 6.18, 11.83], shown_map=(dloglike_doubled<=11.83*1.5), eref=eref)
        if self.target.spectraltype[-8:]=='PowerLaw' or self.target.spectraltype[-9:]=='PowerLaw2':
            self.plot_sed_bowtie(name='dloglike', norms_mesh=norms_mesh, indices_mesh=indices_mesh, dict_meshes_shown={'1sigma':dloglike_doubled<=2.30, '2sigma':dloglike_doubled<=6.18})
        if use_calonly==True:
            self.plot_spectrum_scanned2D(name='dloglike_calonly', norms_mesh=norms_mesh, e2dnde=e2dnde, efluences=efluences, indices_mesh=indices_mesh, zvalues=dloglike_calonly_doubled, cont_levels=[2.30, 6.18, 11.83], shown_map=(dloglike_calonly_doubled<=11.83*1.5), eref=eref)
            self.plot_spectrum_scanned2D(name='dloglike_without_calonly', norms_mesh=norms_mesh, e2dnde=e2dnde, efluences=efluences, indices_mesh=indices_mesh, zvalues=dloglike_without_calonly_doubled, cont_levels=[2.30, 6.18, 11.83], shown_map=(dloglike_without_calonly_doubled<=11.83*1.5), eref=eref)


    def get_sed_bowtie(self, norms_mesh, indices_mesh, dict_meshes_shown, erange=None, neperdec=20):
        if erange is None:
            erange = (self.emin_eval, self.emax_eval)
        loge_min = np.log10(erange[0])
        loge_max = np.log10(erange[1])
        evals = 10 ** np.linspace(loge_min, loge_max, int((loge_max-loge_min)*neperdec+1))
        def e2dnde(e, norm, phindex):
            if self.target.spectraltype[-8:]=='PowerLaw':
                return e*e*norm*pow(e/self.target.spectralpars['Scale'], phindex)
            elif self.target.spectraltype[-9:]=='PowerLaw2':
                return (phindex!=-1)*(e*e*norm*(phindex+1)*pow(e, phindex)/(pow(self.target.spectralpars['UpperLimit'], phindex+1) - pow(self.target.spectralpars['LowerLimit'], phindex+1))) + (phindex==-1)*(e*norm/(np.log(self.target.spectralpars['UpperLimit']) - np.log(self.target.spectralpars['LowerLimit'])))
            else:
                #logger.error('Use PowerLaw or PowerLaw2 instead of {0}'.format(self.target.spectraltype))
                #sys.exit(1)
                return 0
        
        odict_curves_lo = OrderedDict()
        odict_curves_hi = OrderedDict()
        odict_energies = OrderedDict()

        list_e2dnde_maps = [e2dnde(e, norms_mesh, indices_mesh) for e in evals]
            
        for name_shown,shown in dict_meshes_shown.items():
            bound_lo = []
            bound_hi = []
            energies = []
            for ie,e in enumerate(evals):
                e2dnde_shown = list_e2dnde_maps[ie][shown]
                if len(e2dnde_shown)>0:
                    bound_lo.append(min(e2dnde_shown.flatten()))  #np.amin(e2dnde_map + sys.maxint*(1-shown))
                    bound_hi.append(max(e2dnde_shown.flatten())) #np.amax(e2dnde_map - sys.maxint*(1-shown))
                    energies.append(e)
            if len(energies)>0:
                odict_curves_lo[name_shown] = np.array(bound_lo)
                odict_curves_hi[name_shown] = np.array(bound_hi)
                odict_energies[name_shown] = np.array(energies)
            else:
                odict_curves_lo[name_shown] = None
                odict_curves_hi[name_shown] = None
                odict_energies[name_shown] = None
        return (odict_energies, (odict_curves_lo, odict_curves_hi))


    def plot_sed_bowtie(self, name, norms_mesh, indices_mesh, dict_meshes_shown):
        fig, ax = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(10, 5))
        odict_energies,odict_curves = self.get_sed_bowtie(norms_mesh=norms_mesh, indices_mesh=indices_mesh, dict_meshes_shown=dict_meshes_shown)
        if name in self.dct_summary_results:
            self.dct_summary_results[name]['sed_bowtie'] = {'energies':odict_energies, 'curve_lo':odict_curves[0], 'curve_hi':odict_curves[1]}
        #for energies, curve_lo, curve_hi in zip(odict_energies, list_curves[0], list_curves[1]):
        for name_shown in odict_energies.keys():
            energies = odict_energies[name_shown]
            curve_lo = odict_curves[0][name_shown]
            curve_hi = odict_curves[1][name_shown]
            if energies is not None:
                ax[0].fill_between(energies, curve_lo*MEVtoERG, curve_hi*MEVtoERG, alpha=0.2, label=name_shown)
                #ax[0].plot(energies, curve_lo*MEVtoERG, label=name_shown)
                #ax[0].plot(energies, curve_hi*MEVtoERG, label=name_shown)
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].set_ylim((1E-15, 1E-7))
        ax[0].set_xlabel('Energy [MeV]')
        ax[0].set_ylabel(r'$\nu F_{\nu} \ \rm{[erg/cm^2 \cdot s]}$')
        ax[0].grid()
        ax[0].legend(loc=0, fontsize=12, fancybox=True, framealpha=0.5)
        fig.tight_layout() 

        for ff in ('pdf', 'png'):
            path_save = "{0}/{1}_sed_bowtie_{2}{3}.{4}".format(self.dir_work, name, self.target.name, self.suffix, ff)
            fig.savefig(path_save)
            logger.info('{0} has been saved.'.format(path_save))
        fig.clf()


    def plot_spectrum_scanned2D(self, name, norms_mesh, e2dnde, efluences, indices_mesh, zvalues, cont_levels, shown_map, eref):
        indices_min = np.amin(indices_mesh[1:,:][shown_map[1:,:]])
        indices_max = np.amax(indices_mesh[1:,:][shown_map[1:,:]])
        norms_min = np.amin(norms_mesh[1:,:][shown_map[1:,:]])
        norms_max = np.amax(norms_mesh[1:,:][shown_map[1:,:]])
        if e2dnde is not None:
            e2dnde_min = min(e2dnde[1:,:].flatten()[shown_map[1:,:].flatten()]) #e2dnde[1:,:][shown_map[1:,:]])
            e2dnde_max = max(e2dnde[1:,:].flatten()[shown_map[1:,:].flatten()]) #e2dnde[1:,:][shown_map[1:,:]])
        efluences_min = np.amin(efluences[1:,:][shown_map[1:,:]])
        efluences_max = np.amax(efluences[1:,:][shown_map[1:,:]])

        fig, ax = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(15, 5))
        cont = ax[0].contour(indices_mesh, norms_mesh, zvalues, levels=cont_levels, colors='k', linestyles=pMatplot.TPL_LINE) 
        ax[0].set_xlim((indices_min, indices_max)) #(np.amin(indices_mesh + sys.maxint*unshown_map), np.amax(indices_mesh - sys.maxint*unshown_map)))
        ax[0].set_ylim((norms_min, norms_max)) #(np.amin(norms_mesh + sys.maxint*unshown_map), np.amax(norms_mesh - sys.maxint*unshown_map)))
        ax[0].set_yscale('log')
        ax[0].set_xlabel('Spectral index')
        ax[0].set_ylabel('Normalization factor [a.u.]')
        ax[0].grid()

        if e2dnde is not None:
            cont_eflux = ax[1].contour(indices_mesh, e2dnde*MEVtoERG, zvalues, levels=cont_levels, colors='k', linestyles=pMatplot.TPL_LINE)
            ax[1].set_xlim((indices_min, indices_max)) #(np.amin(indices_mesh + sys.maxint*unshown_map), np.amax(indices_mesh - sys.maxint*unshown_map)))
            ax[1].set_ylim((e2dnde_min*MEVtoERG, e2dnde_max*MEVtoERG)) #(np.amin(e2dnde*MEVtoERG + sys.maxint*unshown_map), np.amax(e2dnde*MEVtoERG - sys.maxint*unshown_map)))
            ax[1].set_yscale('log')
            ax[1].set_xlabel('Spectral index')
            ax[1].set_ylabel(r'$\nu F_{{\nu}} \ \rm{{[erg/cm^2 \cdot s]}} \ at \ {eref:.1f} GeV$'.format(eref=eref/1000.))
            ax[1].grid()

        cont_efluence = ax[2].contour(indices_mesh, efluences*MEVtoERG, zvalues, levels=cont_levels, colors='k', linestyles=pMatplot.TPL_LINE)
        ax[2].set_xlim((indices_min, indices_max)) #(np.amin(indices_mesh + sys.maxint*unshown_map), np.amax(indices_mesh - sys.maxint*unshown_map)))
        ax[2].set_ylim((efluences_min*MEVtoERG, efluences_max*MEVtoERG)) #(np.amin(efluences*MEVtoERG + sys.maxint*unshown_map), np.amax(efluences*MEVtoERG - sys.maxint*unshown_map)))
        ax[2].set_yscale('log')
        ax[2].set_xlabel('Spectral index')
        ax[2].set_ylabel(r'$Energy \ fluence \ \rm{[erg/cm^2]}$')
        ax[2].grid()
        #ax[2].legend(loc=0, fontsize=12)
        fig.tight_layout() 

        for ff in ('pdf', 'png'):
            path_save = "{0}/{1}_scanned2D_{2}{3}.{4}".format(self.dir_work, name, self.target.name, self.suffix, ff)
            fig.savefig(path_save)
            logger.info('{0} has been saved.'.format(path_save))
        fig.clf()


    def order_likeratio(self, nseries=0, eref=None, use_calonly=False):
        logger.info('Likelihood ratio ordering is starting...')
        nevt_observed = sum(self.like._Nobs())
        if use_calonly==True and self.calonly is None:
            logger.error('CalOnly data have NOT been added!!')
            use_calonly = False
        if use_calonly==True:
            nevt_observed += sum(self.calonly.onevt)
        if nevt_observed==0:
            nseries = 1
        else:
            nseries = len(self.like.energies)-1
        logger.info('Number of energy bins: {0}'.format(nseries))

        eref = eref if eref is not None else self.target.spectralpars['Scale']
        norms = self.map_norm_range(insert_zero=True)
        indices = self.map_index_range()
        self.like.normPar(self.target.name).setBounds(0, max(10, norms[-1]))
        logger.info("""Profile normalization factor: 
{0}
{1} points.""".format(norms, len(norms)))

        indices_mesh, norms_mesh = np.meshgrid(indices, norms)
        if nseries==1:
            npreds = np.zeros((len(norms), len(indices)))
            #npreds = np.zeros((len(norms), len(indices)))
        else:
            npreds = np.zeros((len(norms), len(indices), nseries)) 
            #npreds = np.zeros((len(norms), len(indices), len(self.like.energies)-1)) 
        efluxes = np.zeros_like(indices_mesh)
        efluences = np.zeros_like(indices_mesh)
        e2dnde = np.zeros_like(indices_mesh)
        prob_unliker = np.zeros_like(indices_mesh)
        prob_liker = np.zeros_like(indices_mesh)

        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
        index_name = 'Index'
        index_idx = self.like.par_index(self.target.name, index_name)        

        # CalOnly data
        if use_calonly==True:
            logger.info('CalOnly data have been loaded.')
            calonly_on, calonly_bkg = self.calonly.get_rebinned_events(np.log10(self.like.energies)) #self.calonly.onevt, self.calonly.bkgevt #
            logger.info('Energy bins: {0}'.format(self.like.energies))
            logger.info('ON events: {0}'.format(calonly_on))
            logger.info('Predicted backgrounds: {0}'.format(calonly_bkg))
            if nseries==1:
                calonly_on = sum(calonly_on) 
                calonly_bkg = sum(calonly_bkg)
        else:
            if nseries==1:
                calonly_on, calonly_bkg = 0, 0
            else:
                calonly_on, calonly_bkg = np.zeros_like(self.like._Nobs()), np.zeros_like(self.like._Nobs())

        if nseries==1:
            nobs = sum(self.like._Nobs())
            nobs += calonly_on
            NOBS_POSSIBLE_MAX = max(int(nobs+5.*sqrt(nobs)), 5)
        else:
            nobs = self.like._Nobs()
            nobs += calonly_on
            NOBS_POSSIBLE_MAX = max(int(sum(nobs)+5.*sqrt(sum(nobs))), 5)
        logger.info('Observed events: {0}'.format(nobs))
        nobs_possible = np.array(range(NOBS_POSSIBLE_MAX))
        nobs_possible_T = nobs_possible[:, np.newaxis]

        #np.log(np.exp(-nobs_possible))
        #pow(nobs_possible, nobs_possible)
        loglike_possible_ideal = np.log(np.exp(-nobs_possible) * ((nobs_possible>0)*pow(nobs_possible, nobs_possible) + (nobs_possible<=0)*1) / scipy.misc.factorial(nobs_possible))
        loglike_possible_ideal_T = loglike_possible_ideal[:, np.newaxis]

        for inorm, jindex in itertools.product(range(len(norms)), range(len(indices))):
            #logger.debug('Normalization: No.{0} {1}'.format(inorm, norms_mesh[inorm][jindex]))
            #logger.debug('Spectral index: No.{0} {1}'.format(jindex, indices_mesh[inorm][jindex]))
            self.like[norm_idx] = norms[inorm]
            self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[norm_idx:norm_idx+1], value=0)
            self.like[index_idx] = indices[jindex]
            self.like.setFreeFlag(srcName=self.target.name, pars=self.like.params()[index_idx:index_idx+1], value=0)
            sum_model = np.zeros_like(self.like._srcCnts(self.like.sourceNames()[0]))
            for srcname in self.like.sourceNames():
                sum_model = sum_model + self.like._srcCnts(srcname)

            # CalOnly
            model_calonly = np.zeros_like(self.like._Nobs()) if nseries>1 else 0
            if use_calonly==True:
                for ie, (e0,e1) in enumerate(zip(self.like.energies[:-1], self.like.energies[1:])):

                    loge_segmented = np.linspace(np.log10(e0), np.log10(e1), int((-np.log10(e0)+np.log10(e1))*20+1.5))
                    line_segmented = 10 ** loge_segmented
                    #if inorm==50:
                    #    logger.info('Energy: {0} - {1}'.format(e0, e1))
                    #    logger.info(loge_segmented)
                    for je in range(len(line_segmented)-1):
                        if nseries>1:
                            model_calonly[ie] = model_calonly[ie] + ( self.like.flux(self.target.name, line_segmented[je], line_segmented[je+1]) * self.calonly.func_onexp((loge_segmented[je]+loge_segmented[je+1])/2.) * 0.68 )
                        else:
                            model_calonly = model_calonly + ( self.like.flux(self.target.name, line_segmented[je], line_segmented[je+1]) * self.calonly.func_onexp((loge_segmented[je]+loge_segmented[je+1])/2.) * 0.68 )

                    # if nseries>1:
                    #     model_calonly[ie] = self.like.flux(self.target.name, e0, e1) * self.calonly.func_onexp((np.log10(e0)+np.log10(e1))/2.) * 0.68
                    # else:
                    #     model_calonly += self.like.flux(self.target.name, e0, e1) * self.calonly.func_onexp((np.log10(e0)+np.log10(e1))/2.) * 0.68

            model_calonly += calonly_bkg
            if nseries==1:
                npreds[inorm][jindex] = sum(sum_model) + model_calonly
            else:
                npreds[inorm][jindex] = sum_model + model_calonly
                #if inorm==50:
                #    logger.info('Regular: {0}'.format(sum_model))
                #    logger.info('CalOnly: {0}'.format(model_calonly))
            efluxes[inorm][jindex] = self.like.energyFlux(self.target.name, self.emin_eval, self.emax_eval)
            efluences[inorm][jindex] = efluxes[inorm][jindex] * self.duration
            if self.target.spectraltype[-8:]=='PowerLaw':
                e2dnde[inorm][jindex] = norms[inorm] * pow(eref/self.target.spectralpars['Scale'], indices[jindex]) * eref * eref
            elif self.target.spectraltype[-9:]=='PowerLaw2':
                e2dnde[inorm][jindex] = norms[inorm] * (indices[jindex]+1) * pow(eref, indices[jindex]) / (pow(self.target.spectralpars['UpperLimit'], indices[jindex]+1) - pow(self.target.spectralpars['LowerLimit'], indices[jindex]+1)) * eref * eref
            else:
                e2dnde[inorm][jindex] = 0
            
            loglike_possible = np.log(np.exp(-npreds[inorm][jindex]) * pow(npreds[inorm][jindex], nobs_possible_T) / scipy.misc.factorial(nobs_possible_T))
            loglikeratio = loglike_possible - loglike_possible_ideal_T
            if nseries==1:
                loglike_possible_meshs = loglike_possible
                loglikeratio_meshs = loglikeratio[:,0]
            if nseries==2:
                loglike_possible_meshs = np.meshgrid(loglike_possible[:,0], loglike_possible[:,1])
                loglikeratio_meshs = np.meshgrid(loglikeratio[:,0], loglikeratio[:,1])
            if nseries==3:
                loglike_possible_meshs = np.meshgrid(loglike_possible[:,0], loglike_possible[:,1], loglike_possible[:,2])
                loglikeratio_meshs = np.meshgrid(loglikeratio[:,0], loglikeratio[:,1], loglikeratio[:,2])
            if nseries==4:
                loglike_possible_meshs = np.meshgrid(loglike_possible[:,0], loglike_possible[:,1], loglike_possible[:,2], loglike_possible[:,3])
                loglikeratio_meshs = np.meshgrid(loglikeratio[:,0], loglikeratio[:,1], loglikeratio[:,2], loglikeratio[:,3])
            if nseries==5:
                loglike_possible_meshs = np.meshgrid(loglike_possible[:,0], loglike_possible[:,1], loglike_possible[:,2], loglike_possible[:,3], loglike_possible[:,4])
                loglikeratio_meshs = np.meshgrid(loglikeratio[:,0], loglikeratio[:,1], loglikeratio[:,2], loglikeratio[:,3], loglikeratio[:,4])

            if nseries==1:
                loglike_possible_sum = loglike_possible_meshs
                loglikeratio_sum = loglikeratio_meshs
            else:
                loglike_possible_sum = sum(loglike_possible_meshs)
                loglikeratio_sum = sum(loglikeratio_meshs)

            if nseries==1:
                loglikeratio_obs = loglikeratio_sum[int(nobs+0.5)]
            if nseries==2:
                loglikeratio_obs = loglikeratio_sum[int(nobs[0]+0.5)][int(nobs[1]+0.5)]
            if nseries==3:
                loglikeratio_obs = loglikeratio_sum[int(nobs[0]+0.5)][int(nobs[1]+0.5)][int(nobs[2]+0.5)]
            if nseries==4:
                loglikeratio_obs = loglikeratio_sum[int(nobs[0]+0.5)][int(nobs[1]+0.5)][int(nobs[2]+0.5)][int(nobs[3]+0.5)]
            if nseries==5:
                loglikeratio_obs = loglikeratio_sum[int(nobs[0]+0.5)][int(nobs[1]+0.5)][int(nobs[2]+0.5)][int(nobs[3]+0.5)][int(nobs[4]+0.5)]
            loglikeratio_best = np.amax(loglikeratio_sum) #self.dct_summary_results['loglike_inversed']
            
            bool_likesum_unliker = loglikeratio_sum<=loglikeratio_obs 
            bool_likesum_liker = loglikeratio_sum>loglikeratio_obs
            if nseries==1:
                prob_unliker[inorm][jindex] = np.sum(np.reshape(np.exp(loglike_possible_sum), (len(bool_likesum_unliker),)) * bool_likesum_unliker)
                prob_liker[inorm][jindex] = np.sum(np.reshape(np.exp(loglike_possible_sum), (len(bool_likesum_liker),)) * bool_likesum_liker)
            else:
                prob_unliker[inorm][jindex] = np.sum(np.exp(loglike_possible_sum) * bool_likesum_unliker)
                prob_liker[inorm][jindex] = np.sum(np.exp(loglike_possible_sum) * bool_likesum_liker)

            #if inorm==0 and jindex==0: #inorm==1 and jindex==5:
                #logger.info('----------')
                #logger.info('loglikeratio_sum: {0}'.format(loglikeratio_sum))
                #logger.info('loglikeratio_obs: {0}'.format(loglikeratio_obs))
                #logger.info('bool_likesum_liker: {0}'.format(bool_likesum_liker))
                #if nseries>1:
                #    logger.info('np.exp(loglike_possible_sum) * bool_likesum_liker: {0}'.format(np.exp(loglike_possible_sum) * bool_likesum_liker))
                #else:
                #    logger.info('np.reshape(np.exp(loglike_possible_sum), (len(bool_likesum_liker),)) * bool_likesum_liker: {0}'.format(np.reshape(np.exp(loglike_possible_sum), (len(bool_likesum_liker),)) * bool_likesum_liker))
                #logger.info('bool_likesum_unliker: {0}'.format(bool_likesum_unliker))
                #if nseries>1:
                #    logger.info('np.exp(loglike_possible_sum) * bool_likesum_unliker: {0}'.format(np.exp(loglike_possible_sum) * bool_likesum_unliker))
                #else:
                #    logger.info('np.reshape(np.exp(loglike_possible_sum), (len(bool_likesum_unliker),)) * bool_likesum_unliker: {0}'.format(np.reshape(np.exp(loglike_possible_sum), (len(bool_likesum_unliker),)) * bool_likesum_unliker))
                # logger.info('loglike_possible:')
                # logger.info(loglike_possible)
                # logger.info('loglike_possible_sum:')
                # logger.info(loglike_possible_sum)
                # logger.info('Sum of loglike_possible_sum:')
                # logger.info(np.sum(np.exp(loglike_possible_sum)))
                # logger.info('loglikeratio_sum:')
                # logger.info(loglikeratio_sum)
                # logger.info('loglikeratio_obs:')
                # logger.info(loglikeratio_obs)
                # logger.info('Unliker fraction: {0}'.format(prob_unliker[inorm][jindex]))
                # logger.info('Liker fraction: {0}'.format(prob_liker[inorm][jindex]))
                # logger.info('Sum of Liker and unliker fraction: {0}'.format(prob_liker[inorm][jindex]+prob_unliker[inorm][jindex]))

        logger.info('-----')
        logger.info("""Npreds: 
{0}""".format(npreds))
        logger.info("""Prob_unliker:
{0}""".format(prob_unliker))
#        logger.info("""Shown Prob_unliker:
#{0}""".format(prob_unliker[prob_unliker>=0.000063]))
        self.dct_summary_results['likeratioordering'] = {}
        self.dct_summary_results['likeratioordering']['normalization'] = norms
        self.dct_summary_results['likeratioordering']['index'] = indices
        self.dct_summary_results['likeratioordering']['eflux'] = efluxes
        self.dct_summary_results['likeratioordering']['efluence'] = efluences
        self.dct_summary_results['likeratioordering']['e2dnde'] = e2dnde 
        self.dct_summary_results['likeratioordering']['pliker'] = prob_liker
        self.dct_summary_results['likeratioordering']['punliker'] = prob_unliker

        self.plot_spectrum_scanned2D(name='likeratioordering', norms_mesh=norms_mesh, e2dnde=e2dnde, efluences=efluences, indices_mesh=indices_mesh, zvalues=prob_unliker, cont_levels=[2.7e-3, 4.55e-2, 3.173e-1], shown_map=(prob_unliker>=0.0000063), eref=eref)
        if self.target.spectraltype[-8:]=='PowerLaw' or self.target.spectraltype[-9:]=='PowerLaw2':
            self.plot_sed_bowtie(name='likeratioordering', norms_mesh=norms_mesh, indices_mesh=indices_mesh, dict_meshes_shown={'2sigma':prob_unliker>=4.55e-2, '1sigma':prob_unliker>=3.173e-1})


    def add_calonly(self, path_onevt, path_onexp, path_bkgevt, nhistogram=0, rclass='R100'):
        file_onevt = ROOT.TFile(path_onevt, "READ")
        tr_onevt = file_onevt.Get('EVENTS_GRB{0}'.format(self.target.name))
        file_onexp = ROOT.TFile(path_onexp, "READ")
        file_bkgevt = ROOT.TFile(path_bkgevt, "READ")
        dir_bkgevt = file_bkgevt.Get('TimeBin{0}'.format(nhistogram))
        #file_offexp = ROOT.TFile(path_offexp, "READ")
        self.calonly = CalOnlyData(tr_onevt=tr_onevt, htg_onexp=file_onexp.Get('htgExp_{0}_projYX'.format(nhistogram)), htg_bkgevt=dir_bkgevt.Get('htgExBKG_CalOnly_{rcla}_PSF68_projE'.format(rcla=rclass)), on_energy=(self.emin, self.emax), on_time=(self.target.met0+self.tmin, self.target.met0+self.tmax), on_classes=DICT_EVCLASSES_BIT_CALONLY[rclass], on_zenith=(0., self.zmax), on_theta=(0., self.thetamax)) #, ebins=np.log10(self.energies),
        if self.calonly.ne_bins<1:
            logger.warning('No CalOnly data in the energy range.')
            return 0
        logger.info('CalOnly data have been added.')
        logger.info('On event: {0}'.format(self.calonly.onevt))
        logger.info('On exposure: {0}'.format(self.calonly._onexp))
        logger.info('Predicted backgrounds: {0}'.format(self.calonly.bkgevt))


def get_allowed_intervals(dict_meshes, select_mesh):
    dict_intervals = {}
    for name_mesh, mesh in dict_meshes.items():
        mesh_selected = mesh[select_mesh]
        dict_intervals[name_mesh] = (np.amin(mesh_selected), np.amax(mesh_selected))
        
    return dict_intervals


class GRBConfig(AnalysisConfig):
    def __init__(self, target, phase, tstop=100000., emin=100., emax=100000., evclass=128, evtype=3, ft2interval=None, deg_roi=12., zmax=100., index_fixed=None, suffix='', tmin_special=None, tmax_special=None, emin_fit=None, emax_fit=None, emin_eval=None, emax_eval=None, binned=False, psForce=False, gti_external=None):

        self.phase = phase
        # Set FT2 interval
        if not ft2interval in ('1s', '30s'):
            if self.phase in ('prompt', 'lightcurve', 'primary', 'intermittent', 'briefslots', 'unified'):
                ft2interval = '1s'
            elif self.phase in ('afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', 'T95to03ks', '01ksto10ks', '03ksto10ks', '01ksto100ks', '03ksto100ks'):
                ft2interval = '30s'
            elif self.phase in ('special', 'lightcurve', 'briefslots'):
                if tmin_special==None or tmax_special==None:
                    logger.critical('Special time window is NOT assigned!!! Set tmin_special and tmax_special as arguments of GRBConfig.')
                    sys.exit(1)
                elif tmax_special-tmin_special>=1000.:
                    ft2interval = '30s'
                else:
                    ft2interval = '1s'
        
        # Define time window
        if self.phase in ('lightcurve', 'special', 'briefslots'):
            tmin = tmin_special
            tmax = tmax_special
        else:
            tphase = define_timephase(target, self.phase, tstop)
            if tphase is (None, None):
                logger.warning('Time phase is NOT valid!')
                tphase = (0,0)
            elif  self.phase=='intermittent' and tphase[0]+10.>=tphase[1]:
                logger.warning('Time phase ("intermittent") is shorter than 10s!')
                tphase = (0,0)
            logger.debug('Time window:{0}'.format(tphase))
            tmin, tmax = tphase

        AnalysisConfig.__init__(self, target, emin=emin, emax=emax, tmin=tmin, tmax=tmax, evclass=128, evtype=3, ft2interval=ft2interval, deg_roi=deg_roi, zmax=zmax, index_fixed=index_fixed, suffix=suffix, emin_fit=emin_fit, emax_fit=emax_fit, emin_eval=emin_eval, emax_eval=emax_eval, binned=binned, psForce=psForce, gti_external=gti_external)

        # Reassign work directory
        self.str_time = 'T{0:0>9.0f}-{1:0>9.0f}ms'.format(self.tmin*1000, self.tmax*1000)
        if self.phase in ('lightcurve', 'special', 'briefslots'):
            self.dir_work = '{base}/{target}/{energy}/{roi}/{phase}/{time}/{spectype}/{binned}'.format(base=PATH_BASEDIR, target=self.target.name, energy=self.str_energy, roi=self.str_roi, phase=self.phase, time=self.str_time, spectype=self.target.spectraltype, binned=self.str_binned) #index=self.str_index)
        else:
            self.dir_work = '{base}/{target}/{energy}/{roi}/{phase}/{spectype}/{binned}'.format(base=PATH_BASEDIR, target=self.target.name, energy=self.str_energy, roi=self.str_roi, phase=self.phase, spectype=self.target.spectraltype, binned=self.str_binned)#, index=self.str_index)


    def download(self, bforce=False):
        self.check_datafiles()
        if self.path_ft1 is not None and bforce==True:
            os.rename(self.path_ft1, self.path_ft1.replace('.fits', '_old.fits'))
            self.path_ft1 = None
        if self.path_ft2 is not None and  bforce==True:
            os.rename(self.path_ft2, self.path_ft2.replace('.fits', '_old.fits'))
            self.path_ft2 = None
        if self.path_ft1==None:
            logger.info('Downloading FT1 data...')
            os.chdir(self.path_dir_data)
            download_fermi_data_grb(self.target.name, lst_ft=[1], path_outdir=self.path_dir_data, reference='GBM' if self.target.table_grb_catalogue['GBM'] is not None else 'LAT')
            #download_fermi_data_grb(self.target.name, lst_ft=[1], path_catalogue=self.path_catalogue, path_outdir=self.path_dir_data)
        if self.path_ft2==None:
            logger.info('Downloading FT2 data...')
            os.chdir(self.path_dir_data)
            download_fermi_data_grb(self.target.name, lst_ft=[2], ft2_interval=self.ft2interval, path_outdir=self.path_dir_data, reference='GBM' if self.target.table_grb_catalogue['GBM'] is not None else 'LAT')
            #download_fermi_data_grb(self.target.name, lst_ft=[2], ft2_interval=self.ft2interval, path_catalogue=self.target.path_catalogue, path_outdir=self.path_dir_data)
        else:
            logger.info("""Downloading data of FT1 and FT2 is skipped.
""")
        # Update data file infomation
        self.check_datafiles()


    def count_axes(self):

        y_model_all = np.zeros_like(self.like._srcCnts(self.like.sourceNames()[0]))
        y_model_target = np.zeros_like(y_model_all)
        y_model_others = np.zeros_like(y_model_all)
        for sourceName in self.like.sourceNames():
            y_model_all = y_model_all + self.like._srcCnts(sourceName)
            if sourceName==self.target.name:
                y_model_target = y_model_target + self.like._srcCnts(sourceName)
            else:
                y_model_others = y_model_others + self.like._srcCnts(sourceName)
        return (y_model_all, y_model_target, y_model_others)

        
    def plot_countspectra_fitted(self):
        fig_cspec_fit, ax_cspec_fit = plt.subplots(2, 1, figsize=(16, 10))
        x_cspec_fit = (self.like.energies[:-1] + self.like.energies[1:])/2.
        y_model_all, y_model_target, y_model_others = self.count_axes()
        ax_cspec_fit[0].loglog(x_cspec_fit, y_model_all,label='Total model')
        ax_cspec_fit[0].loglog(x_cspec_fit, y_model_target,label=self.target.name)
        ax_cspec_fit[0].loglog(x_cspec_fit, y_model_others,label='Sum of other sources')
        if self.binned==True:
            fccube = fits.open(self.path_ccube)
            tbccube = fccube[0].data
            y_obs = np.array([ sum(sum(tbe)) for tbe in tbccube ])
        else:
            y_obs = self.like._Nobs()
        yerr_obs = np.sqrt(y_obs)
        ax_cspec_fit[0].errorbar(x_cspec_fit, y_obs, yerr=yerr_obs, fmt='o',label='Counts')
        ax_cspec_fit[0].legend(loc=0, fontsize=12, fancybox=True, framealpha=0.5)
        ax_cspec_fit[0].set_xlabel('Energy [MeV]')
        ax_cspec_fit[0].set_ylabel('[counts]')
        ax_cspec_fit[0].set_title('RoI of '+self.target.name)
        fig_cspec_fit.savefig("{0}/Count_spectrum_{1}{2}.png".format(self.dir_work, self.target.name, self.suffix))
        fig_cspec_fit.clf()
        return (fig_cspec_fit, ax_cspec_fit)
        #fig_cspec_fit.clf()


class CalOnlyData:
    def __init__(self, tr_onevt, htg_onexp, htg_bkgevt, ebins=None, on_energy=(0, sys.maxint), on_time=(0, sys.maxint), on_classes=(4096, 8192, 16384, 32768), on_zenith=(0., 100.), on_theta=(0., 65.)):
        logger.info('CalOnly data set is constructing...')
        logger.info('Energy limit: {0}'.format(on_energy))
        self.tr_onevt = tr_onevt
        onevt_times = []
        onevt_energies = []
        self.rebin([htg_onexp])
        #self.rebin([htg_onexp, htg_bkgevt])

        if htg_onexp.GetZaxis().GetNbins()!=htg_bkgevt.GetXaxis().GetNbins():
            logger.critical('Bin Numbers of the histograms {0} and {1} do NOT match!!!'.format(htg_onexp.GetName(), htg_bkgevt.GetName()))

        eedges = []
        ecenters = []
        nebins_included = []
        for ie in range(1, htg_onexp.GetZaxis().GetNbins()+1):
            ecenter = htg_onexp.GetZaxis().GetBinCenter(ie)
            if ecenter>=np.log10(on_energy[0]) and ecenter<=np.log10(on_energy[1]):
                eedges.append(htg_onexp.GetZaxis().GetBinLowEdge(ie))
                ecenters.append(ecenter)
                nebins_included.append(ie)
        logger.info('Included bins: {0}'.format(nebins_included))
        if len(nebins_included)<1:
            self.energies=None
            self.energies_center = None
            self.ne_bins = 0
            self.onexp = None
            self.onevt = None
            self.func_onexp = None
            self.onevt_unbinned = None
            self.bkgevt = None
        else:
            eedges.append(htg_onexp.GetZaxis().GetBinUpEdge(max(nebins_included)))
            self.energies = np.array(eedges)
            self.energies_center = np.array(ecenters) #(self.energies[:-1]+self.energies[1:])/2.
            self.ne_bins = len(self.energies)-1

            self._onexp = np.zeros(htg_onexp.GetZaxis().GetNbins()) #self.ne_bins)
            energies_center_original = np.zeros_like(self._onexp)
            for ie in range(len(self._onexp)):
                self._onexp[ie] = htg_onexp.Integral(1, htg_onexp.GetXaxis().GetNbins(), 1, htg_onexp.GetYaxis().GetNbins(), ie+1, ie+1) * 1e4
                energies_center_original[ie] = htg_onexp.GetZaxis().GetBinCenter(ie+1)
            logger.info('CalOnly exposure: {0}'.format(self._onexp))
            self.func_onexp = interpolate.interp1d(energies_center_original, self._onexp, fill_value=0, bounds_error=False)

            self.onevt = np.zeros(self.ne_bins) 
            logger.debug('On event tree: {0}'.format(tr_onevt))
            for evt in tr_onevt:
                if (evt.e>=np.log10(on_energy[0]) and evt.e<np.log10(on_energy[1])) and evt.s in on_classes and evt.FLAG_PSF68==True and (evt.t>=on_time[0] and evt.t<=on_time[1]) and (evt.z>=on_zenith[0] and evt.z<=on_zenith[1]) and (evt.th>=on_theta[0] and evt.th<=on_theta[1]):
                    onevt_times.append(evt.t)
                    onevt_energies.append(evt.e)
                    for iebin, (e0, e1) in enumerate(zip(self.energies[:-1], self.energies[1:])):
                        if e0<=evt.e and evt.e<e1:
                            self.onevt[iebin] += 1
                            break
            logger.info("""CalOnly ON events:
{0}""".format(self.onevt))

            self.onevt_unbinned = {'time':np.array(onevt_times),
                                   'energy': np.array(onevt_energies)}
            self.bkgevt = np.zeros_like(self.onevt) #(htg_bkgevt.GetXaxis().GetNbins())
            for ie,ne in enumerate(nebins_included): #for ie in range(htg_bkgevt.GetXaxis().GetNbins()):
                self.bkgevt[ie] = htg_bkgevt.GetBinContent(ne)
            logger.info("""CalOnly predicted backgrounds:
{0}""".format(self.bkgevt))

            if ebins is not None:
                self.get_rebinned_events(ebins)


    def rebin(self, histograms, energy=28, theta=10, zenith=18):
        dict_newbins = {'X':theta,
                        'Y':zenith,
                        'Z':energy}
        for htg in histograms:
            dict_axes = {'X':htg.GetXaxis(),
                         'Y':htg.GetYaxis(),
                         'Z':htg.GetZaxis()}
            dict_nrebin = {'X':1,
                           'Y':1,
                           'Z':1}
            for ax in ('X', 'Y', 'Z'):
                mbins = dict_axes[ax].GetNbins()
                nbins = dict_newbins[ax]
                if mbins%nbins==0:
                    dict_nrebin[ax] = mbins / nbins
                else:
                    logger.error('Bin number {old} of {htg} {ax} cannnot rebinned into {new}}'.format(htg=htg.GetName(), ax=ax, old=mbins, new=nbins))
                    sys.exit(1)
            if dict_nrebin['X']!=1 or dict_nrebin['Y']!=1 or dict_nrebin['Z']!=1:
                htg.Rebin3D(dict_nrebin['X'], dict_nrebin['Y'], dict_nrebin['Z'])


    def get_rebinned_events(self, ebins):

        onevt_new = np.zeros(len(ebins)-1)
        for e in self.onevt_unbinned['energy']:
            for ie, (e0, e1) in enumerate(zip(ebins[:-1], ebins[1:])):
                if e0<=e and e<e1:
                    onevt_new[ie] += 1
                    break
        #self.onevt = onevt_new

        # onexp_new = np.zeros(len(ebins)-1)
        # for jebin, nexp in enumerate(self.onexp):
        #     ebin_center = (self.energies[jebin]+self.energies[jebin+1]) / 2.
        #     for ie, (e0, e1) in enumerate(zip(ebins[:-1], ebins[1:])):
        #         if e0<=ebin_center and ebin_center<e1:
        #             onexp_new[ie] += nexp
        #             break
        #self.onexp = onexp_new

        bkgevt_new = np.zeros(len(ebins)-1)
        for jebin, nevt in enumerate(self.bkgevt):
            ebin_center = (self.energies[jebin]+self.energies[jebin+1]) / 2.
            for ie, (e0, e1) in enumerate(zip(ebins[:-1], ebins[1:])):
                if e0<=ebin_center and ebin_center<e1:
                    bkgevt_new[ie] += nevt
                    break
        #self.bkgevt = bkgevt_new

        #self.energies = ebins
        #self.energies_center = (self.energies[:-1]+self.energies[1:])/2.
        #self.ne_bins = len(self.energies)-1
        return (onevt_new, bkgevt_new)


##### Utilities #####
def check_required_files(dct_required):
    """Check whether the required files exist.
The input argument is a dictionary of nickname keys and path values.
"""
    sanity = 0
    for name, path in dct_required.iteritems():
        if path==None:
            logger.critical('{0} file is NOT set!!!'.format(name))
            sanity+=1
        elif os.path.exists(path)==False:
            logger.critical('{0} file does not exist!!!'.format(name))
            sanity+=1
    return sanity


def define_timephase(target, phase, tstop=100000.):
    """Returns tmin and tmax for your target and phase.
"""
    print 'Target T90:', target.t90
    print 'Target T95:', target.t95
    # Check phase argement
    if not phase in ("unified", "prompt", "primary", "intermittent", "afterglow", "earlyAG", "lateAG", "farAG", "T95to01ks", "01ksto10ks", "01ksto100ks", "T95to03ks", "03ksto10ks", "03ksto100ks"):
        logger.critical('Phase {0} is NOT avalable!!! Use "unified", "prompt", "afterglow", "earlyAG", "lateAG", "farAG", "T95to01ks", "01ksto10ks", "T95to03ks", "03ksto10ks", "01ksto100ks", "03ksto100ks".')
        sys.exit(1)
    logger.debug(phase)
    tmid_afterglow = target.t95+2.*target.t90
    if phase == 'prompt':
        return (0., target.t95)
    elif phase == 'primary':
        return (0., target.t25+1.5*target.t50)
    elif phase == 'intermittent':
        if target.t25+1.5*target.t50 < target.t95:
            return (target.t25+1.5*target.t50, target.t95)
        else:
            return (None, None)
    elif phase == 'afterglow':
        return (target.t95, tstop)
    elif phase == 'earlyAG':
        return (target.t95, tmid_afterglow)
    elif phase == 'lateAG':
        if tmid_afterglow < tstop:
            return (tmid_afterglow, tstop)
        else:
            logger.critical('tstop is earlier than start of the late afterglow!!!')
            sys.exit(1)
    elif phase == 'farAG':
        return (tstop, 100000.)
    elif phase == "T95to01ks":
        return (target.t95, target.t95+1000.)
    elif phase == "01ksto10ks":
        return (target.t95+1000., 10000.)
    elif phase == "T95to03ks":
        return (target.t95, target.t95+3000.)
    elif phase == "03ksto10ks":
        return (target.t95+3000., 10000.)
    elif phase == "01ksto100ks":
        return (target.t95+1000., 100000.)
    elif phase == "03ksto100ks":
        return (target.t95+3000., 100000.)
    elif phase == 'unified':
        return (0., tstop)


