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
import numpy as np
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
import gt_apps as my_apps
import pyLikelihood
from UnbinnedAnalysis import *
from BinnedAnalysis import *
import FluxDensity
from LikelihoodState import LikelihoodState
import SummedLikelihood
from astropy.io import fits
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
from FindGoodstatPeriods import find_goodstat_periods, get_entries, get_entries_roi
from DownloadFermiData import download_fermi_data_grb
import ReadLTFCatalogueInfo
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
from STLikelihoodAnalysis import get_module_logger
import GetGTI


##### Logger #####
logger = get_module_logger(__name__)


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


##### Target Classes #####
class _AstroTarget:
    def __init__(self, name, spatialmodel, spectraltype, spectralpars, spectralpars_scale, spectralpars_fixed, met0=0, redshift=np.nan):
        self.name = name
        self.spatialmodel = spatialmodel
        self.spectraltype = spectraltype
        if self.spectraltype =='PowerLaw':
            if not set(spectralpars.keys()) == set(['Prefactor', 'Index', 'Scale']):
                logger.error('Spectral parameters are NOT correct for PoweLaw!!!')
                logger.error(spectralpars)
                sys.exit(1)
        elif self.spectraltype =='PowerLaw2':
            if not set(spectralpars.keys()) == set(['Integral', 'Index', 'LowerLimit', 'UpperLimit']):
                logger.error('Spectral parameters are NOT correct for PoweLaw2!!!')
                logger.error(spectralpars)
                sys.exit(1)
        elif self.spectraltype =='ScaleFactor::PowerLaw2':
            if not set(spectralpars.keys()) == set(['Integral', 'Index', 'LowerLimit', 'UpperLimit', 'ScaleFactor']):
                logger.error('Spectral parameters are NOT correct for ScaleFactor::PowerLaw2!!!')
                logger.error(spectralpars)
                sys.exit(1)

        # Name of Normalization parameter
        if spectraltype in ('PowerLaw', 'BrokenPowerLaw', 'ExpCutoff'):
            self.norm_name = 'Prefactor' 
        elif spectraltype in ('PowerLaw2', 'BrokenPowerLaw2', 'ScaleFactor::PowerLaw2'):
            self.norm_name = 'Integral' 

        self.spectralpars = spectralpars
        self.spectralpars_scale = spectralpars_scale
        self.spectralpars_fixed = spectralpars_fixed
        #self.path_pickle = path_pickle
        self.met0 = met0
        self.redshift = redshift
        logger.debug('Target name: {0}'.format(self.name))


class PointTarget(_AstroTarget):
    def __init__(self, name, ra, dec, spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000.}, spectralpars_scale={'Prefactor':1, 'Index':1, 'Scale':1}, spectralpars_fixed=['Scale'], loc_err=0., met0=0, redshift=np.nan):
        _AstroTarget.__init__(self, name, 'PointSource', spectraltype, spectralpars, spectralpars_scale, spectralpars_fixed, met0=met0, redshift=redshift)
        logger.debug('MET0 = {0}'.format(self.met0))
        logger.debug('Target name: {0}'.format(self.name))
        self.ra = ra
        self.dec = dec
        self.loc_err = loc_err


class GRBTarget(PointTarget):
    def __init__(self, name, path_catalogue=GRB_CATALOGUE_LAT, spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000.}):

        tb_gbm = ReadGBMCatalogueInfo.open_table(1, GRB_CATALOGUE_GBM)
        self.path_catalogue = path_catalogue
        if path_catalogue[-5:] == '.fits':
            tb_lat = ReadLTFCatalogueInfo.open_table(1, self.path_catalogue)
            tb_one = ReadLTFCatalogueInfo.select_one_by_name(tb_lat, name)
        elif path_catalogue[-4:] == '.xml':
            tb_lat = ReadLATCatalogueInfo.open_table(GRB_CATALOGUE_LAT)
            tb_one = ReadLATCatalogueInfo.select_one_by_name(tb_lat, name, tb_gbm)
        else:
            logger.critical('Catalogue {0} is not in readable format!!!'.format(path_catalogue))
        if spectraltype=='PowerLaw':
            spectralpars_scale = {'Prefactor':1, 'Index':1, 'Scale':1}
            spectralpars_fixed = ['Scale']
        elif spectraltype=='PowerLaw2':
            spectralpars_scale = {'Integral':1, 'Index':1, 'LowerLimit':1, 'UpperLimit':1}
            spectralpars_fixed = ['LowerLimit', 'UpperLimit']
        elif spectraltype=='ScaleFactor::PowerLaw2':
            spectralpars_scale = {'Integral':1, 'Index':1, 'LowerLimit':1, 'UpperLimit':1, 'ScaleFactor':1}
            spectralpars_fixed = ['LowerLimit', 'UpperLimit', 'ScaleFactor']
        elif spectraltype=='ExpCutoff':
            spectralpars_scale = {'Prefactor':1, 'Index':1, 'Ebreak':1, 'P1':1, 'P2':1, 'P3':1}
            spectralpars_fixed = ['Scale', 'P2', 'P3']
        elif spectraltype=='BrokenPowerLaw':
            spectralpars_scale = {'Prefactor':1, 'Index1':1, 'Index2':1, 'BreakValue':1}
            spectralpars_fixed = ['BreakValue']

        PointTarget.__init__(self, name, float(tb_one['RA']), float(tb_one['DEC']), spectraltype=spectraltype, spectralpars=spectralpars, spectralpars_scale=spectralpars_scale, spectralpars_fixed=spectralpars_fixed, met0=pMETandMJD.ConvertMjdToMet(float(tb_one['TRIGGER_TIME'])), redshift=(float(tb_one['REDSHIFT']) if float(tb_one['REDSHIFT'])>0 else np.nan))

        if path_catalogue[-5:] == '.fits':
            tb_gbm_one = ReadGBMCatalogueInfo.select_one_by_name(tb_gbm, tb_one['GBM_assoc_key'])
            self.t05 = tb_one['T90_START']
            self.t90 = tb_one['T90']
            self.t25 = tb_gbm_one['T50_START']
            self.t50 = tb_gbm_one['T50']
        elif path_catalogue[-4:] == '.xml':
            tb_gbm_one = ReadGBMCatalogueInfo.select_one_by_name(tb_gbm, tb_one['GBM']['NAME'])
            self.t05 = tb_one['GBM']['T90_START']
            self.t90 = tb_one['GBM']['T90']
            self.t25 = tb_one['GBM']['T50_START']
            self.t50 = tb_one['GBM']['T50']
        self.met05 = self.met0 + self.t05
        self.t95 = self.t05 + self.t90
        self.met95 = self.met05 + self.t90
        self.met25 = self.met0 + self.t25
        logger.debug('Target name: {0}'.format(self.name))


##### Analysis Classes #####
PATH_BASEDIR = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone'

class AnalysisConfig:
    def __init__(self, target, emin, emax, tmin, tmax, evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., index_fixed=None, suffix='', emin_fit=None, emax_fit=None, emin_eval=None, emax_eval=None, binned=False, psForce=False, roi_checked=5.0, gti_external=None):
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
        self.index_fixed = index_fixed
        self.rad_margin = rad_margin
        self.suffix = '' if suffix=='' else '_'+suffix
        self.roi_checked = roi_checked # Radius within which number of events is checked.

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
        self.gti_filter = '(DATA_QUAL==1)&&(LAT_CONFIG==1)&&(ANGSEP(RA_ZENITH,DEC_ZENITH,{ra},{dec}) + {rad} < {zen})'.format(ra=self.target.ra, dec=self.target.dec, rad=12., zen=self.zmax)
        #self.gti_filter = ' (DATA_QUAL==1)&&(LAT_CONFIG==1)' #'(DATA_QUAL>0)&&(LAT_CONFIG==1)'
        self.gti_external = gti_external
        self.roicut = True #False

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

        my_apps.maketime['scfile'] = self.path_ft2
        my_apps.maketime['filter'] = self.gti_filter
        my_apps.maketime['roicut'] = 'yes' if self.roicut==True else 'no'
        my_apps.maketime['evfile'] = self.path_filtered
        my_apps.maketime['outfile'] = self.path_filtered_gti
        my_apps.maketime.run()
        logger.info("""Making GTI finished.
""")
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


    def diffuse_responses(self, bforce=False):
        logger.info("""Calculating diffuse source responses is starting...""")
        check_required_files({'FT2':self.path_ft2, 'Filtered GTI event':self.path_filtered_gti, 'Source model':self.path_model_xml})

        my_apps.diffResps['evfile'] = self.path_filtered_gti
        my_apps.diffResps['scfile'] = self.path_ft2
        my_apps.diffResps['srcmdl'] = self.path_model_xml
        my_apps.diffResps['irfs'] = self.irfs
        my_apps.diffResps.run()
        logger.info("""Calculating diffuse source responses finished.
""")


    def setup(self, force={'download':False, 'filter':False, 'maketime':False, 'evtbin':False, 'livetime':False, 'exposure':False, 'model_3FGL_sources':False, 'diffuse_responses':False, 'srcmaps':False}, skip_zero_data=False):
        """Perform all preparation before fitting. Namely, downloading, filtering, making GTI, calculating livetime and exposure, modeling 3FGL sources, making diffuse source responses.
"""
        self.set_directories()
        self.download(bforce=force['download'])
        self.filter(bforce=force['filter'])
        self.maketime(bforce=force['maketime'])
        logger.info('Roughly {0} events.'.format(self.nevt_rough))
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
        #self.like.setEnergyRange(self.emin_fit, self.emax_fit)
        #logger.debug('New energy bounds: {0}'.format(self.like.energies))
        if self.binned==False:
            self.like.reset_ebounds(self.energies)
            logger.info('Energy bound has changed to {0}'.format(self.like.energies))

        # Target source
        if not self.target.name in self.like.sourceNames():
            logger.info("""Adding {0}.""".format(self.target.name))
            target_src = pyLike.PointSource(0, 0, self.obs.observation) #self.like.observation.observation)
            pl = pyLike.SourceFactory_funcFactory().create(self.target.spectraltype)
            if self.target.spectraltype=='PowerLaw':
                pl.setParamValues((self.target.spectralpars['Prefactor'], self.target.spectralpars['Index'], self.target.spectralpars['Scale']))
               #Prefactor
                prefactor = pl.getParam("Prefactor")
                prefactor.setBounds(0.0, 1e-4)
                prefactor.setScale(self.target.spectralpars_scale['Prefactor'])
                pl.setParam(prefactor)
              # Index
                indexPar = pl.getParam("Index")
                indexPar.setBounds(-9.0, 3.0)
                indexPar.setScale(self.target.spectralpars_scale['Index'])
                pl.setParam(indexPar)
               #Scale
                escale = pl.getParam("Scale")
                escale.setBounds(100., 100000.)
                escale.setScale(self.target.spectralpars_scale['Scale'])
                pl.setParam(escale)
                pl.setParamAlwaysFixed('Scale')
                #target_src.setSpectrum(pl)

            elif self.target.spectraltype=='PowerLaw2':
                pl.setParamValues((self.target.spectralpars['Integral'], self.target.spectralpars['Index'], self.target.spectralpars['LowerLimit'], self.target.spectralpars['UpperLimit']))
               #Integral
                prefactor = pl.getParam("Integral")
                prefactor.setBounds(0.0, 1e-8)
                prefactor.setScale(self.target.spectralpars_scale['Integral'])
                pl.setParam(prefactor)
              # Index
                indexPar = pl.getParam("Index")
                indexPar.setBounds(-9.0, 3.0)
                indexPar.setScale(self.target.spectralpars_scale['Index'])
                pl.setParam(indexPar)
               #LowerLimit
                ell = pl.getParam("LowerLimit")
                ell.setBounds(100., 100000.)
                ell.setScale(self.target.spectralpars_scale['LowerLimit'])
                pl.setParam(ell)
                pl.setParamAlwaysFixed('LowerLimit')
               #UpperLimit
                eul = pl.getParam("UpperLimit")
                eul.setBounds(100., 100000.)
                eul.setScale(self.target.spectralpars_scale['UpperLimit'])
                pl.setParam(eul)
                pl.setParamAlwaysFixed('UpperLimit')
                #target_src.setSpectrum(pl)

            elif self.target.spectraltype=='ScaleFactor::PowerLaw2':
                #pl.setParamValues((self.target.spectralpars['Integral'], self.target.spectralpars['Index'], self.target.spectralpars['LowerLimit'], self.target.spectralpars['UpperLimit']))#, self.target.spectralpars['ScaleFactor']))
               #Integral
                prefactor = pl.getParam("Integral")
                prefactor.setScale(self.target.spectralpars_scale['Integral'])
                prefactor.setValue(self.target.spectralpars['Integral'])
                prefactor.setBounds(0, 1e6)
                pl.setParam(prefactor)
              # Index
                indexPar = pl.getParam("Index")
                indexPar.setValue(self.target.spectralpars['Index'])
                indexPar.setScale(self.target.spectralpars_scale['Index'])
                indexPar.setBounds(-9.0, 3.0)
                pl.setParam(indexPar)
               #LowerLimit
                ell = pl.getParam("LowerLimit")
                ell.setScale(self.target.spectralpars_scale['LowerLimit'])
                ell.setValue(self.target.spectralpars['LowerLimit'])
                ell.setBounds(100., 100000.)
                pl.setParam(ell)
                pl.setParamAlwaysFixed('LowerLimit')
               #UpperLimit
                eul = pl.getParam("UpperLimit")
                eul.setScale(self.target.spectralpars_scale['UpperLimit'])
                eul.setValue(self.target.spectralpars['UpperLimit'])
                eul.setBounds(100., 100000.)
                pl.setParam(eul)
                pl.setParamAlwaysFixed('UpperLimit')
               #ScaleFactor
                scf = pl.getParam("ScaleFactor")
                scf.setScale(self.target.spectralpars_scale['ScaleFactor'])
                scf.setValue(self.target.spectralpars['ScaleFactor'])
                scf.setBounds(0., 100.)
                pl.setParam(scf)
                pl.setParamAlwaysFixed('ScaleFactor')

            elif self.target.spectraltype=='ExpCutoff':
                pl.setParamValues((self.target.spectralpars['Prefactor'], self.target.spectralpars['Index'], self.target.spectralpars['Scale'], self.target.spectralpars['Ebreak'], self.target.spectralpars['P1'], self.target.spectralpars['P2'], self.target.spectralpars['P3'])) # Prefactor, Index, Scale, Ebreak, P1, P2, P3
                indexPar = pl.getParam("Index")
                indexPar.setBounds(-9.0, 3.0)
                pl.setParam(indexPar)
                prefactor = pl.getParam("Prefactor")
                prefactor.setBounds(0.0, 1e-4)
                prefactor.setScale(1)
                pl.setParam(prefactor)
                #escale = pl.getParam("Scale")
                #ebreak.setBounds(100., 10000.)
                #pl.setParam(escale)
                ebreak = pl.getParam("Ebreak")
                ebreak.setBounds(-1000., 1000.)
                pl.setParam(ebreak)
                #target_src.setSpectrum(pl)
                ecutoff = pl.getParam("P1")
                ecutoff.setBounds(100., 1000000.)
                pl.setParam(ecutoff)
                #target_src.setSpectrum(pl)

            elif self.target.spectraltype=='BrokenPowerLaw':
                pl.setParamValues((self.target.spectralpars['Prefactor'], self.target.spectralpars['Index1'], self.target.spectralpars['Index2'], self.target.spectralpars['BreakValue'])) # Prefactor, Index1, Index2, BreakValue
                indexPar1 = pl.getParam("Index1")
                indexPar1.setBounds(-9.0, 3.0)
                pl.setParam(indexPar1)
                indexPar2 = pl.getParam("Index2")
                indexPar2.setBounds(-9.0, 3.0)
                pl.setParam(indexPar2)
                prefactor = pl.getParam("Prefactor")
                prefactor.setBounds(0.0, 1e-4)
                prefactor.setScale(1)
                pl.setParam(prefactor)
                ebreak = pl.getParam("BreakValue")
                ebreak.setBounds(100., 100000.)
                pl.setParam(ebreak)

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

        self.obs = UnbinnedObs(self.path_filtered_gti, self.path_ft2, expMap=self.path_exposure, expCube=self.path_livetime, irfs=self.irfs)
        self.like = UnbinnedAnalysis(self.obs, path_model, optimizer='NewMinuit')
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
            if self.target.spectraltype=='ExpCutoff':
                eb_idx = self.like.par_index(self.target.name, 'Ebreak')
                self.like.freeze(eb_idx)
                logger.info(self.like.params()[eb_idx])
                logger.info('Free parameter of {src}: {free}'.format(src=self.target.name, free=self.like.freePars(self.target.name)))

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


    def summarize_fit_results(self):
        """Summarize the results after fitting and return a dictonary. The contents are the model parameters, photon flux, TS, return code.
"""
        # Model parameters
        if self.target.spectraltype=='PowerLaw':
            model_pars = ('Prefactor', 'Index', 'Scale')
        elif self.target.spectraltype=='PowerLaw2':
            model_pars = ('Integral', 'Index', 'LowerLimit', 'UpperLimit')
        elif self.target.spectraltype=='ScaleFactor::PowerLaw2':
            model_pars = ('Integral', 'Index', 'LowerLimit', 'UpperLimit', 'ScaleFactor')
        elif self.target.spectraltype=='ExpCutoff':
            model_pars = ('Prefactor', 'Index', 'Scale', 'Ebreak', 'P1', 'P2', 'P3')
        for name_param in model_pars:
            param = self.like.model[self.target.name].funcs['Spectrum'].getParam(name_param)
            self.dct_summary_results[name_param] = {'value':param.value(), 'error':param.error()}

        # Flux
        flux_and_err = self.eval_flux_and_error(self.target.name, self.emin_eval, self.emax_eval)

        # TS
        name = self.target.name
        logger.debug('TS of {0}:'.format(name))
        self.dct_summary_results['TS'] = self.like.Ts(str(name))

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
        e2 = eref if (eref is not None or self.target.spectraltype!='PowerLaw') else self.target.spectralpars['Scale']

        norm_value = self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).value()
        norm_error = self.like.model[self.target.name].funcs['Spectrum'].getParam(self.target.norm_name).error()
        norm_idx = self.like.par_index(self.target.name, self.target.norm_name)
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
        self.like.normPar(self.target.name).setBounds(xvals[0], xvals[-1])
        logger.info("""Profile normalization factor: 
{0}
{1} points.""".format(xvals, len(xvals)))

        # Index parameter
        index_name = 'Index'
        index_idx = self.like.par_index(self.target.name, index_name)
        index_value = self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).value()
        #if not index_value==index_value:
        #    logger.error('Index value is NOT valid!!! {0}'.format(index_value))
        #    index_value = -2.0
        index_error = self.like.model[self.target.name].funcs['Spectrum'].getParam(index_name).error()
        #if not index_error==index_error:
        #    logger.error('Index error is NOT valid!!! {0}'.format(index_value))
        #    index_error = 0.0
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
                
            if self.target.spectraltype=='PowerLaw':
                v0['dnde'] = norm_value * (e2 / self.target.spectralpars['Scale']) ** index_values[str_index_assumed]
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
                #index_fit = index_values[str_index_assumed]
                if str_index_assumed == 'free':
                    try:
                        loglike1 = -self.like.fit(verbosity=0,covar=False,optObject=self.likeobj)
                    except RuntimeError:
                        logger.error('RuntimeError!!!')
                        logger.error('Fitting with free index failed!!')
                        failed_freeIndex = True
                        sys.exit(1)
                else:
                    loglike1 = -self.like()
                logger.debug('No.{ip} normalization factor = {xval:.3E}, expected count = {nph}, loglikelihood = {ll:.3E}'.format(ip=i, xval=x, nph=self.like._npredValues(), ll=loglike1)) #NpredValue(self.target.name)))
                o[str_index_assumed]['dloglike'][i] = loglike1 - loglike0
                o[str_index_assumed]['loglike'][i] = loglike1

            # limits
            if str_index_assumed=='free' and failed_freeIndex==True:
                break
            limits[str_index_assumed]['norm'] = get_parameter_limits(xval=xvals, loglike=o[str_index_assumed]['dloglike'], cl_limit=cl_2sigma, cl_err=cl_1sigma)
            if self.target.spectraltype=='PowerLaw':
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
            print 'Target T90:', target.t90
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
            download_fermi_data_grb(self.target.name, lst_ft=[1], path_outdir=self.path_dir_data)
            #download_fermi_data_grb(self.target.name, lst_ft=[1], path_catalogue=self.path_catalogue, path_outdir=self.path_dir_data)
        if self.path_ft2==None:
            logger.info('Downloading FT2 data...')
            os.chdir(self.path_dir_data)
            download_fermi_data_grb(self.target.name, lst_ft=[2], ft2_interval=self.ft2interval, path_outdir=self.path_dir_data)
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
        return (fig_cspec_fit, ax_cspec_fit)
        #fig_cspec_fit.clf()


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
