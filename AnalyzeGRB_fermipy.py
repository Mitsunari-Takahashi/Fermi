#!/usr/bin/env python

import os
import os.path
import subprocess
#subprocess.call(['slacsetup', '11-05-02'])
import matplotlib as mpl
import matplotlib.pyplot as plt
from fermipy.utils import get_parameter_limits
from fermipy.gtanalysis import GTAnalysis
#import FluxDensity
#from LikelihoodState import LikelihoodState
#from fermipy.gtutils import BinnedAnalysis, SummedLikelihood
#import pyLikelihood as pyLike
import numpy as np
from math import log10
from pLsList import ls_list
from pReadGBMCatalogueInfo import ReadGBMCatalogueOneLine
#import re
import click
from FindCrossEarthlimb import find_cross_earthlimb
from FindGoodstatPeriods import find_goodstat_periods, get_entries


def AnalyzeGRB_fermipy(name, ft1_candidates, ft2_candidates, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, ebinedge, path_outdir, mode, catalogues, goodstat, shiftenergies): #, eranges=[[100, 316228]]):
    NAME_TGT = name
    dct_grb = ReadGBMCatalogueOneLine(NAME_TGT, '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/Highest-GBM-fluence-GRBs.csv')
    print dct_grb

    RA = dct_grb['ra']
    DEC = dct_grb['dec']
    T0 = dct_grb['trigger_time']
    if (dct_grb['z'] is not '') and (dct_grb['z'] is not '-') and (dct_grb['z'] is not '?'):
        REDSHIFT = dct_grb['z']
        #if thirtygev==True:
        #    eranges.append([31622.8/(1+REDSHIFT), 316228/(1+REDSHIFT)])
        if shiftenergies==True:
            for eset in ebinedge:
                eset[0] = eset[0]/(1+REDSHIFT)
                eset[1] = eset[1]/(1+REDSHIFT)
            print 'Shifted energy ranges:', ebinedge
    else:
        REDSHIFT = 1
    SUFFIX = ''
    if suffix!='':
        SUFFIX = '_' + suffix

    # Model catalogues
    if len(catalogues)>0:
        str_catalogues = 'catalogs : ['
        for cata in catalogues:
            str_catalogues = str_catalogues + cata + ','
        str_catalogues = str_catalogues[:-1] + ']'
    else:
        str_catalogues = ''

    ZCUT = 100.
    lst_tbin = []
    mstatag = 0
    if mode == 'unified':
        lst_tbin = [[tmin, tmax]]
    else:
        validtimes = find_cross_earthlimb(ft2_candidates[0], RA, DEC, T0+tmin, T0+tmax, ZCUT, T0)
        print validtimes
        for (ivt, vt) in enumerate(validtimes):
            if goodstat>0:
                if mode=='prompt':
                    periods_goodstat = find_goodstat_periods(ft1_candidates[0], T0+vt[0], T0+vt[1], goodstat)
                    print periods_goodstat
                    for igs in range(len(periods_goodstat)-1):
                        lst_tbin.append([periods_goodstat[igs]-T0, periods_goodstat[igs+1]-T0])
                elif mode=='afterglow':
                    mstatag += get_entries(ft1_candidates[0], vt[0]+T0, vt[1]+T0)
                    if ivt==0:
                        lst_tbin.append([vt[0]])
                    if mstatag>=goodstat:
                        lst_tbin[-1].append(vt[1])
                        mstatag = 0
                        if ivt<len(validtimes)-1:
                            lst_tbin.append([validtimes[ivt+1][0]])
                    if ivt==len(validtimes)-1:
                        if len(lst_tbin)>1:
                            lst_tbin = lst_tbin[:-1]
                        if len(lst_tbin[-1])<2:
                            lst_tbin[-1].append(vt[1])
                        else:
                            lst_tbin[-1][1] = vt[1]
            else:
                if ivt==0:
                    if tmin>=vt[0] and tmin<vt[1]:
                        lst_tbin = [[tmin]]
                    else:
                        lst_tbin = [[vt[0]]]
                    if tbinedges is not None and len(tbinedges)>0:
                        for tedge in tbinedges:
                            if vt[1]-tedge>tedge:
                                if tedge < vt[1]:
                                    lst_tbin[-1].append(tedge)
                                    lst_tbin.append([tedge])
                            else:
                                lst_tbin[-1].append(vt[1])
                    else:
                        lst_tbin[-1].append(vt[1])
                else: 
                    lst_tbin.append([vt[0]])
                    lst_tbin[-1].append(vt[1])
            
    print 'Time bin edges:', lst_tbin

    LST_RAN_TIME = lst_tbin
    NRAN_TIME = len(LST_RAN_TIME)

    #LST_RAN_ENR_LOG = ebinedge # in MeV
    #NRAN_ENR = len(LST_RAN_ENR_LOG)

    #LST_SED_ITEM = ['index', 'ts', 'e_min', 'e_max', 'e_ctr', 'e_ref', 'flux', 'flux_err', 'flux_err_lo', 'flux_err_hi', 'eflux', 'eflux_err', 'eflux_err_lo', 'eflux_err_hi', 'dnde', 'dnde_err_lo', 'dnde_err_hi', 'dnde_ul95', 'ref_flux', 'ref_dnde', 'ref_dnde_e_min', 'ref_dnde_e_max']
    LST_SED_ITEM_CSV = ['e_min', 'e_max', 'e_ref', 'index', 'ts', 'flux', 'flux_err_lo', 'flux_err_hi', 'eflux', 'eflux_err_lo', 'eflux_err_hi', 'ref_flux', 'ref_dnde', 'ref_dnde_e_min', 'ref_dnde_e_max']

    print 'Going to start standard analysis of', NAME_TGT
    if path_outdir=='.':
        path_outdir = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/' + NAME_TGT
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)

    for eedges in ebinedge:
        strenergies = 'E{0:0>6}-{1:0>6}MeV'.format(int(eedges[0]+0.5), int(eedges[1]+0.5))
        print '%%%%%%%%%%%%%%%%%%'
        print strenergies
        print '%%%%%%%%%%%%%%%%%%'
        for itime in range(NRAN_TIME):
            strtime = 'T{0:0>6}-{1:0>6}s'.format(int(0.5+LST_RAN_TIME[itime][0]), int(0.5+LST_RAN_TIME[itime][1]))
            print '=====', strtime, '====='
            path_subdir = '{0}/{1}/{2}'.format(path_outdir, strenergies, strtime)
            if not os.path.exists(path_subdir):
                os.makedirs(path_subdir)
            os.chdir(path_subdir)
            path_anlaysis = '{0}/fit_model{1}.npy'.format(path_subdir, SUFFIX)
            if os.path.exists(path_anlaysis) and force==False:
                print 'Loading previous analysis...'
                gta = GTAnalysis.create(path_anlaysis)
            else:
                str_path_lt = path_subdir+'/ltcube_00.fits'
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
  roiwidth   : 20.0
  binsz      : 0.1
  binsperdec : 4

selection :
  emin : {4}
  emax : {5}
  zmax    : {12}
  evclass : 128 # 8
  evtype  : 3
  tmin    : {6}
  tmax    : {7}
  filter  : null
  ra      : {8}
  dec     : {9}

gtlike:
  edisp : False
  irfs : 'P8R2_SOURCE_V6' #'P8R2_TRANSIENT020E_V6' #'P8R2_SOURCE_V6'
  edisp_disable : ['isodiff','galdiff']

model:
  src_radius  : 22.0
  galdiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/gll_iem_v06.fits'
  isodiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/iso_P8R2_SOURCE_V6_v06.txt'
#  isodiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/iso_P8R2_TRANSIENT020E_V6_v06.txt'
  {13}
  sources :
    - {{ 'name' : 'GRB{10}', 'ra' : {8}, 'dec' :{9}, 'SpectrumType' : 'EblAtten::PowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index' : {{ value : 2.0, min : -1.0, max : 8.0, scale : -1, free : '0' }}, 'LowerLimit' : {{ value : {4}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {5}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {11}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}, 'SpatialModel': 'PointSource'}}
""".format(path_subdir, ft1_candidates[0], ft2_candidates[0], str_path_lt, eedges[0], eedges[1], int(T0+LST_RAN_TIME[itime][0]), int(T0+LST_RAN_TIME[itime][1]), RA, DEC, NAME_TGT, REDSHIFT, ZCUT, str_catalogues)

                with open("{0}/config.yaml".format(path_subdir), 'w') as conf:
                    conf.write(str_config)
                with open("{0}/config.yaml".format(path_subdir), 'r') as conf:
                    print conf

                print 'Setting up...'
                gta = GTAnalysis('{0}/config.yaml'.format(path_subdir),logging={'verbosity' : 3})

                gta.setup()
                gta.optimize()
                gta.print_roi()
            
                gta.free_sources(free=False)
                if (dct_grb['z'] is not '') and (dct_grb['z'] is not '-') and (dct_grb['z'] is not '?'):
                    gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index'])
                else:
                    gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index', 'redshift'])

                print 'Fitting...'
                gta.fit()
                gta.write_roi('fit_model'+SUFFIX)
                gta.print_roi()
                print ' Fitting finished.'
                continue
    
            print '===== Fitting parameters ====='
            c = np.load('{0}/fit_model{1}.npy'.format(path_subdir, SUFFIX)).flat[0]
            np_src = c['sources']['GRB'+NAME_TGT]
            #print 'Model counts:'
            #print np_src['model_counts']
            #print 'N_pred:'
            #print np_src['npred']
            
            norm_lims95 = get_parameter_limits(np_src['norm_scan'], np_src['loglike_scan'], 0.95)
            norm_lims68 = get_parameter_limits(np_src['norm_scan'], np_src['loglike_scan'], 0.68)
            print '* Integral:', np_src['param_values'][0], '+/-', np_src['param_errors'][0]
            print '  95% limits:', np_src['param_values'][0]*norm_lims95['ll']/norm_lims95['x0'], '-', np_src['param_values'][0]*norm_lims95['ul']/norm_lims95['x0']
            print '  68% limits:', np_src['param_values'][0]*norm_lims68['ll']/norm_lims95['x0'], '-', np_src['param_values'][0]*norm_lims68['ul']/norm_lims95['x0']
            print '* Index:', np_src['param_values'][1], '+/-', np_src['param_errors'][1]
            print ''
            for (ipar, par) in enumerate(np_src['param_names']):
                print par, ':', np_src['param_values'][ipar], '+/-', np_src['param_errors'][ipar]

            str_lc = """#start:stop:ts:Integral:Integral_err:Integral_ul95:Integral_ll95:Integral_ul68:Integral_ll68:Index:Index_err:flux:flux_err
{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}
""".format(LST_RAN_TIME[itime][0], LST_RAN_TIME[itime][1], np_src['ts'], np_src['param_values'][0], np_src['param_errors'][0], np_src['param_values'][0]*norm_lims95['ul']/norm_lims95['x0'], np_src['param_values'][0]*norm_lims95['ll']/norm_lims95['x0'], np_src['param_values'][0]*norm_lims68['ul']/norm_lims68['x0'], np_src['param_values'][0]*norm_lims68['ll']/norm_lims68['x0'], np_src['param_values'][1], np_src['param_errors'][1], np_src['flux'], np_src['flux_err'])
            with open("{0}/GRB{1}_{2}_lc.csv".format(path_subdir, NAME_TGT, strtime), 'w') as text:
                print str_lc
                text.write(str_lc)

    #TS map
            if skipts==False:
                model_ts = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
                maps_ts = gta.tsmap('fit_ts',model=model_ts, make_plots=True)
                gta.plotter.make_tsmap_plots(maps_ts, roi=gta.roi)
            
            # logemin_extrapolated = 4.75
            # logemax_extrapolated = 5.5
            # loge_extrapolated = np.linspace(logemin_extrapolated, logemax_extrapolated, 50)
            # bowtie_extrapolate = gta.bowtie(NAME_TGT, fd=None, loge=loge_extrapolated)
            # print bowtie_extrapolate
                
    #SED
            if skipsed==False:
                if c['sources']['GRB'+NAME_TGT]['npred']>0:
                    print 'Going to SED analysis...'
                    sed  = gta.sed('GRB'+NAME_TGT, use_local_index=True, make_plots=True, outfile='sed{0}.fits'.format(SUFFIX))
            #print '----- MODEL FLUX -----'
            #for (ie, enr) in enumerate(sed['model_flux']['energies']):
            #    print enr, ' ', sed['model_flux']['dnde'][ie]
            #print '----------------------'
                    str_outcsv = ''
                    for (item, tem) in enumerate(LST_SED_ITEM_CSV):
                        str_outcsv = str_outcsv + str(tem)
                        if item < len(LST_SED_ITEM_CSV)-1:
                            str_outcsv = str_outcsv + ','
                        else:
                            str_outcsv = str_outcsv + """
"""
                    gta.plotter.make_sed_plots(sed, roi=gta.roi, prefix=strtime)
#                     for ienran in range(NRAN_ENR-1):
#                         for (item, tem) in enumerate(LST_SED_ITEM_CSV):
#                             str_outcsv = str_outcsv + str(sed[tem][ienran])
#                             if item < len(LST_SED_ITEM_CSV)-1:
#                                 str_outcsv = str_outcsv + ','
#                             else:
#                                 str_outcsv = str_outcsv + """
# """

#                     with open("{0}/GRB{1}_{2}_sed.csv".format(path_subdir, NAME_TGT, strtime), 'w') as text:
#                         print str_outcsv
#                         text.write(str_outcsv)
                else:
                    print 'Predicted count is', c['sources']['GRB'+NAME_TGT]['npred']
                    print 'Skipping SED analysis...'

    # Residual map
            if skipresid==False:
                model_resid = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
                maps_resid = gta.residmap('fit_resid',model=model_resid, make_plots=True)
                gta.plotter.make_residmap_plots(maps_resid, roi=gta.roi)

    # Light curve
    #if skiplc==False:
     #   gta.free_source('GRB'+NAME_TGT, free=False, pars=['Index', 'redshift'])
      #  lc = gta.lightcurve('GRB'+NAME_TGT, make_plots=True, time_bins=LST_TIME_MET, free_params=['Integral'], write_npy=True, write_fits=True)


@click.command()
@click.argument('name', type=str)
@click.argument('tmin', type=float)
@click.argument('tmax', type=float)
@click.option('-s', '--suffix', type=str, default='')
@click.option('-f', '--force', is_flag=True)
@click.option('--skipts', is_flag=True)
@click.option('--skipsed', is_flag=True)
@click.option('--skipresid', is_flag=True)
@click.option('--logemin', default=2.)
@click.option('--logemax', default=5.5)
@click.option('--nebindecade', default=0)
@click.option('--tbinedges', '-t', multiple=True, default=None, type=float)
@click.option('--outpath', '-o', default='.')
@click.option('--mode', '-m', type=click.Choice(['prompt', 'afterglow', 'unified']))
@click.option('--catalogues', '-c', multiple=True, default=None, type=str)
@click.option('--goodstat', '-g', type=int, default=0)
@click.option('--thirtygev', is_flag=True)
@click.option('--shiftenergies', is_flag=True)
def main(name, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, logemin, logemax, nebindecade, outpath, mode, catalogues, goodstat, thirtygev, shiftenergies):
    lst_ebin = []
    if nebindecade>0:
        nebin = int((logemax-logemin)*nebindecade+0.5)
        webin = (logemax-logemin)/float(nebin)
        for iebin in range(1, 1+nebin):
            lst_ebin.append([logemin, logemin+iebin*webin])
    else:
        lst_ebin.append([logemin, logemax])
    print 'Energy bin edges:', lst_ebin
    ft1_candidates = [None]
    ft2_candidates = [None]
    if mode=='prompt':
        ft1_candidates = ls_list(outpath+'/*-ft1*.fits')
        ft2_candidates = ls_list(outpath+'/*-ft2*.fits')
    elif mode=='afterglow' or mode=='unified':
        ft1_candidates = ls_list(outpath+'/*_PH??.fits')
        ft2_candidates = ls_list(outpath+'/*_SC??.fits')
    AnalyzeGRB_fermipy(name, ft1_candidates, ft2_candidates, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, lst_ebin, outpath, mode, catalogues, goodstat, shiftenergies) #, thirtygev, shiftenergies)


if __name__ == '__main__':
    main()
