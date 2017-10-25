#!/usr/bin/env python

import sys
import os
import os.path
import matplotlib as mpl
#mpl.use('tkagg')
mpl.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib import gridspec
import gt_apps as my_apps
from pyLikelihood import *
from UnbinnedAnalysis import *
from BinnedAnalysis import *
from Composite2 import Composite2
from CompositeLikelihood import CompositeLikelihood
#from bdlikeSED import *
import click
from astropy.io import fits
import numpy as np
import itertools
import gc
import ReadLTFCatalogueInfo
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
import GetGTI
from pLsList import ls_list
import pickle_utilities

mpl.rcParams['font.size'] = 25

lst_extertnal_parameters = ['gbm_fluence',
                            #'gbm_fluence_per_t50',
                            #'gbm_t90',
                            'gbm_flux1024',
                            'gbm_flux64',
                            'epeak_band',
                            'lat_flux_lowE'
                            ]
dct_col_gbm_catarogue = {'gbm_fluence':'FLUENCE',
                         #'gbm_t90':'T90',
                         'gbm_flux1024':'FLUX_1024',
                         'gbm_flux64':'FLUX_64',
                         'epeak_band':'FLNC_BAND_EPEAK'}


def get_category(chara, ncategory, path_dct='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/QuantiledGRBs_longonly3_LC.pickle'):
    print 'Loading {0}...'.format(path_dct)
    dct = pickle_utilities.load(path_dct)
    if ncategory==0:
        a = []
        for l in dct[chara]:
            a = a+l
        return a
    else:
        return dct[chara][ncategory-1]


def run_composite2(chara_cat, path_input, path_outdir, names_params_tied_universal=['Prefactor', 'Index'], ncat_analyzed=0, str_suffix='', eedges=None, binned=False, longonly=False, tolerr=180., path_dct_category='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_LC.pickle', fig_forms=('png', 'pdf')):

    # Open table
    #tb_ltf = ReadLTFCatalogueInfo.open_table()
    #if longonly==True:
    #    tb_ltf = ReadLTFCatalogueInfo.select_long(tb_ltf)
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    #if longonly==True:
    #    tb_gbm = ReadGBMCatalogueInfo.select_long(tb_gbm)
    tb_lat = ReadLATCatalogueInfo.open_table()
    lst_lat = ReadLATCatalogueInfo.read_all(tb_lat, tb_gbm)
    if longonly==True:
        lst_lat = ReadLATCatalogueInfo.select_long(lst_lat)
    lst_lat = ReadLATCatalogueInfo.select_small_error(lst_lat, tolerr)

    # Definition of GBM fluence categories
    FLUENCE_CUT = [1.09e-04, 3.06e-05] #[1.45E-4, 3.70E-5] # Top 10%, 35%
    NCATEGORIES_FLUENCE = len(FLUENCE_CUT)+1
    #dct_category_fluence = {}
#    rh_fluence_weightNobs = ROOT.TH1D('roohtg', 'GBM Fluence', 100, -7, -2)

    path_base = os.getcwd()
    #os.chdir(path_outdir)
    lst_name_subdir = ['plots', 'xml', 'fits']
    for name_subdir in lst_name_subdir:
        path_subdir = '{0}/{1}'.format(path_outdir, name_subdir)
        if not os.path.exists(path_subdir):
            os.makedirs(path_subdir)
    if str_suffix != '':
        str_suffix = '_' + str_suffix
    irfs = 'P8R2_SOURCE_V6' # To be used with P301 data
    optimizer='Minuit'
    level = 2.71
    CompositeLike = Composite2(optimizer=optimizer)
    like={}
    targets = []
    lst_fluence_gbm = []
    lst_fluence_gbm_err = []
    #lst_nobs_lat = []
    targets_analyzed = []
    lst_ncat_analyzed = get_category(chara_cat, ncat_analyzed, path_dct_category)
    dct_chara_title = {'GBM_FLUENCE':'Prompt fluence in GBM', 
                       #'gbm_fluence_per_t50':'Prompt fluence / T50 in GBM', 
                       #'gbm_t90':'GBM T90', 
                       'GBM_FLUX_1024':'Prompt peak flux in GBM for 1024ms', 
                       'GBM_FLUX_64':'Prompt peak flux in GBM for 64ms',
                       #'epeak_band':'Peak energy of prompt Band component', 
                       #'gbm_intermittent':'GBM intermittent duration',
                       #'spec_index':'Spectral index in 0.178 - 5.62 GeV', 
                       #'flux_gev':'Flux in 0.178 - 5.62 GeV', 
                       'LC_INDEX':'Light curve index in 0.1 - 100 GeV'
                       }
    norm_name = 'Integral'
    if isinstance(chara_cat, int):
        dct_chara_title[chara_cat] = 'Simulation No.{0:0>5}'.format(chara_cat)
    print dct_chara_title[chara_cat], 'Category No.', ncat_analyzed-1
    print lst_ncat_analyzed
    print len(lst_ncat_analyzed), 'GRBs.'

    nenergies = int(np.log10(eedges[1]/eedges[0])*4+0.5)
    energies = 10 ** np.linspace(np.log10(eedges[0]), np.log10(eedges[1]), nenergies+1)

    srcs_virtual = {} # For calculating sum of LAT fluence

    for (itarget, target_info) in enumerate(lst_lat):
#    for (itarget, path_target) in enumerate(lst_inputs):
        target = target_info['GRBNAME'] #name_base[3:12]
        path_target = '{0}/{1}/GRB{0}_P8_P302_BASE_T00-999-101000_r030'.format(target, path_input)
        path_base, name_base = os.path.split(path_target)
        targets.append(target)
        print '##### No.{0} {1} #####'.format(itarget, target)
        #dct_category_fluence[target] = judge_category_fluence(tb, target, FLUENCE_CUT) 

        if target not in lst_ncat_analyzed:
            print 'skipped.'
            continue

        targets_analyzed.append(target)
        ltcube = '/'.join((path_base, name_base+'_ft1_ltCube.fits'))
        expMap = '/'.join((path_base, name_base+'_ft1_expMap.fits'))
        srcMap = '/'.join((path_base, name_base+'_ft1_srcmap.fits'))
        ccube = '/'.join((path_base, name_base+'_ft1_ccube.fits'))
        #srcModel = '/'.join((path_base, name_base+'_ft1_model.xml'))
        srcModel = '/'.join((path_base, name_base+'_ft1_model_new.xml'))
        evt = '/'.join((path_base, name_base+'_ft1_filtered_gti.fits'))
        sc = '/'.join((path_base, '../../../../..', name_base.replace('_P8_P302_BASE_T00-999-101000_r030', '_T00-999-101000_ft2-30s.fits')))

        if itarget==0:
            print 'Files of the first target.'
            print '  Event:', evt
            print '  Spacecraft:', sc
            print '  Livetime cube:', ltcube
            print '  Exposure map:', expMap
            print '  Source model:', srcModel
        srcs_virtual[target] = {}
        srcs_virtual[target]['time'] = GetGTI.get_duration(evt)

        try:
            if binned==False:
                like[target] = unbinnedAnalysis(evfile=evt,
                                            scfile=sc,
                                            expmap=expMap,
                                            expcube=ltcube,
                                            irfs=irfs,
                                            srcmdl=srcModel,
                                            optimizer=optimizer)
                if energies is None:
                    energies = (100., 100000.)
                like[target].reset_ebounds(energies)
                print 'Energy bound has reset to {0}'.format(like[target].energies)
            else:
                like[target] = binnedAnalysis(irfs=irfs,
                                          expcube=ltcube,
                                          srcmdl=srcModel,
                                          optimizer=optimizer,
                                          cmap=srcMap,
                                          bexpmap=expMap)
        except RuntimeError:
        # Diffuse responses
            my_apps.diffResps['evfile'] = evt
            my_apps.diffResps['scfile'] = sc
            my_apps.diffResps['srcmdl'] = srcModel
            my_apps.diffResps['irfs'] = irfs
            my_apps.diffResps.run()

            if binned==False:
                like[target] = unbinnedAnalysis(evfile=evt,
                                            scfile=sc,
                                            expmap=expMap,
                                            expcube=ltcube,
                                            irfs=irfs,
                                            srcmdl=srcModel,
                                            optimizer=optimizer)
                if energies is None:
                    energies = (100., 100000.)
                like[target].reset_ebounds(energies)
                print 'Energy bound has reset to {0}'.format(like[target].energies)
            else:
                like[target] = binnedAnalysis(irfs=irfs,
                                          expcube=ltcube,
                                          srcmdl=srcModel,
                                          optimizer=optimizer,
                                          cmap=srcMap,
                                          bexpmap=expMap)

        for source in like[target].sourceNames():
            if source not in (target):
                like[target].normPar(source).setFree(False)
        sys.stdout.flush()

        CompositeLike.addComponent(like[target])
        sys.stdout.flush()
        #del like[target]
        gc.collect()
    print 'Analyzed GRBs:', targets_analyzed

    # Tying parameters universaly
    tiedParams_universal = {}
    for par in names_params_tied_universal:
        tiedParams_universal[par] = []
        for target in lst_ncat_analyzed:
            if target in targets:
                tiedParams_universal[par].append(tuple([like[target], target, par]))
        CompositeLike.tieParameters(tiedParams_universal[par])
    print '* Parameters tied universaly:'
    print tiedParams_universal

    #minuit = eval("pyLike.%s(CompLike.composite)"%optimizer)
    #minuit.setStrategy(2)
    #likeobj = pyLike.NewMinuit(like.logLike)
    fit_result = CompositeLike.fit(covar=True,tol=1.e-5,optimizer=optimizer)
    print '== Fitting result =='
    print fit_result

    dct_params = {}
    for tiedpar in names_params_tied_universal:
        dct_params[tiedpar] = {}
        dct_params[tiedpar]['value'] = like[targets_analyzed[0]].model[targets_analyzed[0]].funcs['Spectrum'].getParam(tiedpar).value()
        dct_params[tiedpar]['error'] = like[targets_analyzed[0]].model[targets_analyzed[0]].funcs['Spectrum'].getParam(tiedpar).error()
        print '{n}: {v} +/- {e}'.format(n=tiedpar, v=dct_params[tiedpar]['value'], e=dct_params[tiedpar]['error'])
    print ''

    x_stacked = (like[targets_analyzed[0]].energies[:-1] + like[targets_analyzed[0]].energies[1:])/2.
    print len(x_stacked), 'energy bins.'
    model_sum_stacked = np.zeros_like(like[targets_analyzed[0]]._srcCnts(like[targets_analyzed[0]].sourceNames()[0]))
    model_grb_stacked = np.zeros_like(model_sum_stacked)
    model_others_stacked = np.zeros_like(model_sum_stacked)
    nobs_sum_stacked = np.zeros_like(model_sum_stacked)
    #srcs_virtual = {} # Store characteristics under assuming resultant parameters of stacking analysis
    #npred = {}

    # Loop over GRBs
    for target in targets_analyzed:
        print target
        sys.stdout.flush()
        model_sum = np.zeros_like(model_sum_stacked)
        model_grb = np.zeros_like(model_sum_stacked)
        model_others = np.zeros_like(model_sum_stacked)
        if binned==False:
            nobs_sum = like[target]._Nobs()
        else:
            fccube = fits.open(ccube)
            tbccube = fccube[0].data
            nobs_sum = np.array([ sum(sum(tbe)) for tbe in tbccube ])
        nobs_sum_stacked = nobs_sum_stacked + nobs_sum

        for sourceName in like[target].sourceNames():
            model_sum = model_sum + like[target]._srcCnts(sourceName)
            if sourceName==target:
                model_grb = model_grb + like[target]._srcCnts(sourceName)
            else:
                model_others = model_others + like[target]._srcCnts(sourceName)
        model_sum_stacked = model_sum_stacked + model_sum
        model_grb_stacked = model_grb_stacked + model_grb
        model_others_stacked = model_others_stacked + model_others

        #npred[target] = like[target].NpredValue(target)
        #srcs_virtual[target] = {}
        srcs_virtual[target]['npred'] = like[target].NpredValue(target)
        srcs_virtual[target]['flux'] = {}
        srcs_virtual[target]['flux']['value'] = like[target].flux(target)
        srcs_virtual[target]['flux']['error'] = srcs_virtual[target]['flux']['value']*dct_params[norm_name]['error']/dct_params[norm_name]['value']#like[target].fluxError(target) }
        srcs_virtual[target]['eflux'] = {}
        srcs_virtual[target]['eflux']['value'] = like[target].energyFlux(target)
        srcs_virtual[target]['eflux']['error'] = srcs_virtual[target]['eflux']['value']*dct_params[norm_name]['error']/dct_params[norm_name]['value']#like[target].fluxError(target) }
        #srcs_virtual[target]['eflux'] = { 'value': like[target].energyFlux(target), 'error':like[target].energyFluxError(target) }
        srcs_virtual[target]['fluence'] = {}
        for t in ('nominal', 'conservative'):
            srcs_virtual[target]['fluence'][t] = {}
            for l in ('value', 'error'):
                srcs_virtual[target]['fluence'][t][l] = srcs_virtual[target]['time'][t]*srcs_virtual[target]['eflux'][l]

    # Count spectrum
    if fig_forms is not None:
        fig_stacked, ax_stacked = plt.subplots(2, 1, figsize=(16, 10))
        print 'Plotting count spectra...'
        print 'X-axis:', x_stacked
        print 'Y-axis:', model_sum_stacked
        print 'Y-axis:', model_grb_stacked
        print 'Y-axis:', model_others_stacked
        print 'Y-axis:', nobs_sum_stacked
        ax_stacked[0].set_ylim(0.1, 2E4)
        ax_stacked[0].loglog(x_stacked, model_sum_stacked, label='Sum of models')
        ax_stacked[0].loglog(x_stacked, model_grb_stacked, label='GRBs')
        ax_stacked[0].loglog(x_stacked, model_others_stacked, label='Others')
        ax_stacked[0].errorbar(x_stacked, nobs_sum_stacked, yerr=np.sqrt(nobs_sum_stacked), fmt='o',label='Counts')
        ax_stacked[0].legend(loc=1, fontsize=20)
        ax_stacked[0].set_title('{chara} Category No.{cat}'.format(chara=dct_chara_title[chara_cat], cat=ncat_analyzed))
        ax_stacked[0].set_ylabel('[counts]')
        ax_stacked[0].set_xticklabels([])
        ax_stacked[0].grid(ls='-', lw=0.5, alpha=0.2)
       #ax_stacked[0].set_xlabel(r'$\log_{10}Energy$ [MeV]')

        print 'Plotting residuals...'
        print 'X-axis:', x_stacked
        resid_stacked = (nobs_sum_stacked - model_sum_stacked) / model_sum_stacked
        print 'Y-axis:', resid_stacked
        resid_stacked_err = np.sqrt(nobs_sum_stacked) / model_sum_stacked
        ax_stacked[1].set_xscale('log')
        ax_stacked[1].set_yscale('linear')
        ax_stacked[1].errorbar(x_stacked, resid_stacked, yerr=resid_stacked_err, fmt='o')
        ax_stacked[1].axhline(0.0,ls=':')
        ax_stacked[1].grid(ls='-', lw=0.5, alpha=0.2)
        #ax_stacked].set_xlabel(r'$\log{10}Energy$ [MeV]')
        ax_stacked[1].set_xlabel(r'$\log_{10}Energy \rm{[MeV]}$')
        ax_stacked[1].set_ylabel('Fractional residual')

        fig_stacked.tight_layout()
        fig_stacked.subplots_adjust(hspace=0)
        ax_stacked[1].set_yticks([y for y in ax_stacked[1].get_yticks() if y<ax_stacked[1].get_ylim()[1]])

        for ff in fig_forms:
            fig_stacked.savefig('{0}/plots/StackedSpectrum{1}_category{2}.{3}'.format(path_outdir, str_suffix, ncat_analyzed, ff))

    dct_summary = {'dloglike_inv':fit_result, 'targets':targets_analyzed, 'parameters':dct_params, 'srcs_virtual':srcs_virtual }
    pickle_utilities.dump('{0}/Summary_StackedAnalysis{1}_category{2}.pickle'.format(path_outdir, str_suffix, ncat_analyzed), dct_summary)

    return dct_summary #(fit_result, dct_params, npred)


@click.command()
@click.argument('inputs', type=str)
@click.option('--pathout', '-o', type=str, default='.')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--tieuniv', multiple=True, type=str)
@click.option('--category', type=click.Choice(['0', '1', '2', '3']), help='0: all GRBs 1,2,3: only GRBs of each category')
@click.option('--chara', type=click.Choice(['GBM_FLUENCE', 'GBM_FLUX_1024', 'GBM_FLUX_64', 'LC_INDEX']), help='GRB characterristics which categorize them.')
#@click.option('--chara', type=click.Choice(['gbm_fluence', 'gbm_t90', 'spec_index', 'flux_gev', 'lc_index', 'gbm_intermittent', 'gbm_flux1024', 'gbm_flux64', 'epeak_band', 'gbm_fluence_per_t50']), help='GRB characterristics which categorize them.')
@click.option('--energies', '-e', type=(float, float), default=(100., 100000.))
#@click.option('--escale', type=float, default=1000.)
@click.option('--longonly', '-l', is_flag=True)
@click.option('--tolerr', '-t', type=float, default=180.)
@click.option('--binned', '-b', is_flag=True)
@click.option('--quantile', '-q', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_LC.pickle')
def main(chara, inputs, pathout, category, tieuniv, suffix, energies, binned, longonly, tolerr, quantile):
        if category==None:
            categories = [0,1,2,3]
        else:
            categories = [int(category)]

        for category in categories:
            print '##### Category {0} #####'.format(category)
            run_composite2(chara, inputs, '{0}/Category{1}'.format(pathout, category), tieuniv, category, suffix, energies, binned, longonly, tolerr, path_dct_category=quantile)
            print ''


if __name__ == '__main__':
    main()
