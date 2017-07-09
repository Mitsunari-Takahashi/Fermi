#!/usr/bin/env python

import os
import os.path
import subprocess
#subprocess.call(['slacsetup', '11-05-02'])
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import matplotlib.pyplot as plt
from math import log10
from pLsList import ls_list
#import re
import click

def AnalyzeGRB_fermipy(name, ra, dec, t0, suffix, force, skipts, skiplc, ebinedge, tbinedge, path_outdir):
#def AnalyzeGRB_fermipy(name, ra, dec, negativedec, t0, suffix, force, skipts, skiplc):
    NAME_TGT = name
    RA = ra
    DEC = dec
    #if negativedec==True:
     #   DEC = - dec
    T0 = t0
    SUFFIX = ''
    if suffix!='':
        SUFFIX = '_' + suffix
    #LST_RAN_TIME = [[0, 15], [15, 30], [30, 60], [60, 120], [120,240], [240,480], [480,960], [960,1920], [1920, 3840], [3840, 7680]]
    LST_RAN_TIME = tbinedge #[0, 15, 30, 60, 120, 240, 480, 960, 1920, 3840, 7680]
    LST_TIME_MET = [x+T0 for x in LST_RAN_TIME]
    #LST_TIME_MET = [x+T0 for inner_list in LST_RAN_TIME for x in inner_list]
    NRAN_TIME = len(LST_RAN_TIME)-1

#LST_RAN_ENR = [[100,1000.]]#, [100, 1000]] # in MeV
    #LST_RAN_ENR_LOG = [2., 2.5, 3., 3.5, 4.0, 4.5, 5.0] # in MeV
    LST_RAN_ENR_LOG = ebinedge # in MeV
    NRAN_ENR = len(LST_RAN_ENR_LOG)

    LST_SED_ITEM = ['index', 'ts', 'e_min', 'e_max', 'e_ctr', 'e_ref', 'flux', 'flux_err', 'flux_err_lo', 'flux_err_hi', 'eflux', 'eflux_err', 'eflux_err_lo', 'eflux_err_hi', 'dnde', 'dnde_err_lo', 'dnde_err_hi', 'dnde_ul95', 'ref_flux', 'ref_dnde', 'ref_dnde_e_min', 'ref_dnde_e_max']
    LST_SED_ITEM_CSV = ['e_min', 'e_max', 'e_ref', 'index', 'ts', 'flux', 'flux_err_lo', 'flux_err_hi', 'eflux', 'eflux_err_lo', 'eflux_err_hi', 'ref_flux', 'ref_dnde', 'ref_dnde_e_min', 'ref_dnde_e_max']

    print 'Going to start standard analysis of', NAME_TGT
    if path_outdir==None:
        path_outdir = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/' + NAME_TGT
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)
        
    path_anlaysis = '{0}/fit_model{1}.npy'.format(path_outdir, SUFFIX)
    if os.path.exists(path_anlaysis) and force==False:
        #return 0
        print 'Loading previous analysis...'
        gta = GTAnalysis.create(path_anlaysis)
    else:
        ft1_candidates = ls_list('/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/DATA/FITS/LAT/{0}_P8_P302_ALL_ft1_*.fits'.format(NAME_TGT))
#         str_lst_ft1 = ''
#         for path_ft1 in ft1_candidates:
#             str_lst_ft1 = str_lst_ft1 + path_ft1 + """
# """
#         with open("{0}/ft1.lst".format(path_outdir), 'w') as ft1_lst:
#             ft1_lst.write(str_lst_ft1)

        #if len(ft1_candidates)>1:
        #    print 'Ambiguity of input ft1 data files!!!'
        #    return 1
        ft2_candidates = ls_list('/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/DATA/FITS/LAT/{0}_P8_P302_BASE_ft2_*.fits'.format(NAME_TGT))
        #if len(ft2_candidates)>1:
        #    print 'Ambiguity of input ft2 data files!!!'
        #    return 1
        str_path_lt = path_outdir+'/ltcube_00.fits'
        if not os.path.exists(str_path_lt):
            str_path_lt = 'Null'

        str_config = """fileio:
  outdir : {0}

data:
  evfile : {1}
  scfile : {2}
  ltcube : {3}
  #ltcube : ltcube_00.fits

binning:
  roiwidth   : 30.0
  binsz      : 0.1
  binsperdec : 5

selection :
  emin : {4}
  emax : {5}
  zmax    : 90
  evclass : 32 #128 32
  evtype  : 3
  tmin    : {6}
  tmax    : {7}
  filter  : null
  ra      : {8}
  dec     : {9}

gtlike:
  edisp : True
  irfs : 'P8R2_SOURCE_V6' #'P8R2_TRANSIENT010E_V6' #'P8R2_SOURCE_V6'
  edisp_disable : ['isodiff','galdiff']

model:
  src_radius  : 22.0
#  galdiff  : '/afs/slac/g/glast/ground/GLAST_EXT/redhat6-x86_64-64bit-gcc44/diffuseModels/v5r0/gll_iem_v06.fits'
  galdiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/gll_iem_v06.fits'
#  isodiff  : '/afs/slac/g/glast/ground/GLAST_EXT/redhat6-x86_64-64bit-gcc44/diffuseModels/v5r0/iso_P8R2_SOURCE_V6_v06.txt'
  isodiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/iso_P8R2_SOURCE_V6_v06.txt'
  catalogs : ['3FGL', '2FHL']
  sources :
#    - {{ 'name' : 'GRB{10}', 'ra' : {8}, 'dec' :{9}, 'SpectrumType' : 'EblAtten::BrokenPowerLaw2', 'Integral' : {{ value : 2.0, max : 100., min : 1e-2, scale : !!float 1e-7, free : '1' }}, 'Index1' : {{ value : -2.1, max : 0.0, min : -5.0, scale : 1, free : '1' }}, 'Index2' : {{ value : -1.6, max : 0.0, min : -5.0, scale : 1, free : '1' }}, 'LowerLimit' : {{ value : 0.1, max : 100., min : 0.1, scale : !!float 1e3, free : '0' }}, 'UpperLimit' : {{ value : 100.0, max : 1000., min : 10., scale : !!float 1e3, free : '0' }}, 'BreakValue' : {{ value : 5, max : 50, min : 0.5, scale : !!float 1e3, free : '1' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : 1.17, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}, 'SpatialModel': 'PointSource'}}
    - {{ 'name' : 'GRB{10}', 'ra' : {8}, 'dec' :{9}, 'SpectrumType' : 'EblAtten::PowerLaw2', 'Integral' : {{ value : 2.0, max : !!float 1e5, min : !!float 1e-5, scale : !!float 1e-8, free : '1' }}, 'Index' : {{ value : 1.8, min : 0.0, max : 5.0, scale : -1, free : '1' }}, 'LowerLimit' : {{ value : 0.1, max : 100., min : 0.1, scale : !!float 1e3, free : '0' }}, 'UpperLimit' : {{ value : 10.0, max : 1000., min : 10., scale : !!float 1e3, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : 1.17, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}, 'SpatialModel': 'PointSource'}}
#    - {{ 'name' : 'GRB{10}', 'ra' : {8}, 'dec' :{9}, 'SpectrumType' : 'BrokenPowerLaw', 'Prefactor' : !!float 1E-10, 'Index1' : 2.1, 'Index2' : 1.6, 'BreakValue' : 10000., 'SpatialModel': 'PointSource'}}
""".format(path_outdir, ft1_candidates[0], ft2_candidates[0], str_path_lt, pow(10, LST_RAN_ENR_LOG[0]), pow(10, LST_RAN_ENR_LOG[-1]), T0+LST_RAN_TIME[0], T0+LST_RAN_TIME[-1], RA, DEC, NAME_TGT)

        with open("{0}/config.yaml".format(path_outdir), 'w') as conf:
            conf.write(str_config)
        with open("{0}/config.yaml".format(path_outdir), 'r') as conf:
            for row in conf:
                print row

        print 'Setting up...'
        gta = GTAnalysis('{0}/config.yaml'.format(path_outdir),logging={'verbosity' : 3})

        gta.setup()
        gta.optimize(npred_threshold=0.01)
        gta.print_roi()
            
        gta.free_sources(free=False)
        #gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index1', 'Index2','BreakValue'])
        gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index'])
        #gta.free_source('GRB'+NAME_TGT, pars=['Prefactor', 'Index1', 'Index2', 'BreakValue', 'P1'], free=True)
        #gta.free_source('GRB'+NAME_TGT, pars=['Prefactor', 'Index1', 'Index2', 'BreakValue'], free=True)
        #gta.free_source('GRB'+NAME_TGT, pars='Eabs', free=False)
        #gta.free_source('galdiff', free=True)
        gta.print_roi()

        print 'Fitting...'
        gta.fit()
        gta.write_roi('fit_model'+SUFFIX)
        gta.print_roi()
        print ' Fitting finished.'
        return 0
    
    print '===== Fitting parameters ====='
    c = np.load('{0}/fit_model{1}.npy'.format(path_outdir, SUFFIX)).flat[0]
    np_src = c['sources']['GRB'+NAME_TGT]
    print 'Model counts:'
    print np_src['model_counts']
    print 'N_pred:'
    print np_src['npred']
    for (ipar, par) in enumerate(np_src['param_names']):
        print par, ':', np_src['param_values'][ipar], '+/-', np_src['param_errors'][ipar]

    #TS map
    if skipts==False:
        model_ts = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
        maps_ts = gta.tsmap('fit_ts',model=model_ts, make_plots=True)
        gta.plotter.make_tsmap_plots(maps_ts, roi=gta.roi)

    #SED
    if c['sources']['GRB'+NAME_TGT]['npred']>0:
        print 'Going to SED analysis...'
        sed  = gta.sed('GRB'+NAME_TGT, use_local_index=True, make_plots=True, outfile='sed{0}.fits'.format(SUFFIX), loge_bins=LST_RAN_ENR_LOG) #loge_bins=[LST_RAN_ENR_LOG[ie][0], LST_RAN_ENR_LOG[ie][1]], 
        print '----- MODEL FLUX -----'
        for (ie, enr) in enumerate(sed['model_flux']['energies']):
            print enr, ' ', sed['model_flux']['dnde'][ie]
        print '----------------------'
        str_outcsv = ''
        for (item, tem) in enumerate(LST_SED_ITEM_CSV):
            str_outcsv = str_outcsv + str(tem)
            if item < len(LST_SED_ITEM_CSV)-1:
                str_outcsv = str_outcsv + ','
            else:
                str_outcsv = str_outcsv + """
"""
        gta.plotter.make_sed_plots(sed, roi=gta.roi)
        for ienran in range(NRAN_ENR-1):
            for (item, tem) in enumerate(LST_SED_ITEM_CSV):
                str_outcsv = str_outcsv + str(sed[tem][ienran])
                if item < len(LST_SED_ITEM_CSV)-1:
                    str_outcsv = str_outcsv + ','
                else:
                    str_outcsv = str_outcsv + """
"""

        with open("{0}/GRB{1}.csv".format(path_outdir, NAME_TGT), 'w') as text:
            print str_outcsv
            text.write(str_outcsv)
    else:
        print 'Predicted count is', c['sources']['GRB'+NAME_TGT]['npred']
        print 'Skipping SED analysis...'

    # Residual map
    model_resid = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    maps_resid = gta.residmap('fit_resid',model=model_resid, make_plots=True)
    gta.plotter.make_residmap_plots(maps_resid, roi=gta.roi)

    # Light curve
    if skiplc==False:
        lc = gta.lightcurve('GRB'+NAME_TGT, make_plots=True, time_bins=LST_TIME_MET)


@click.command()
@click.argument('name', type=str)
@click.argument('ra', type=float)
@click.argument('dec', type=float)
@click.argument('t0', type=float)
@click.argument('tmin', type=float)
@click.argument('tmax', type=float)
@click.option('-s', '--suffix', type=str, default='')
@click.option('-f', '--force', is_flag=True)
@click.option('--skipts', is_flag=True)
@click.option('--skiplc', is_flag=True)
@click.option('--logemin', default=2.)
@click.option('--logemax', default=5.)
@click.option('--nebindecade', default=1.)
@click.option('--tbinedges', '-t', multiple=True, default=None)
@click.option('--outpath', '-o', default=None)
#@click.option('-n', '--negativedec', is_flag=True)
def main(name, ra, dec, t0, tmin, tmax, tbinedges, suffix, force, skipts, skiplc, logemin, logemax, nebindecade, outpath):
    nebin = int((logemax-logemin)*nebindecade+0.5)
    webin = (logemax-logemin)/float(nebin)
    lst_ebin = []
    for iebin in range(nebin+1):
        lst_ebin.append(logemin+iebin*webin)
    print 'Energy bin edges:', lst_ebin
    lst_tbin = [tmin]
    if tbinedges is not None:
        for tedge in tbinedges:
            lst_tbin.append(tedge)
    lst_tbin.append(tmax)
    print 'Time bin edges:', lst_tbin
    AnalyzeGRB_fermipy(name, ra, dec, t0, suffix, force, skipts, skiplc, lst_ebin, lst_tbin, outpath)


if __name__ == '__main__':
    main()
