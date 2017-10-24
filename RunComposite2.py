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

def judge_category_fluence(tb, name, lst_cut):
    tb = ReadLTFCatalogueInfo.select_gbm_exist(tb)
    tb1 = ReadLTFCatalogueInfo.select_one_by_name(tb, name)
    ncategory = len(lst_cut)
    for ic in range(len(lst_cut)):
        ncategory -= int(tb1['FLUENCE']>=lst_cut[ic])
    print 'Fluence:', tb1['FLUENCE'], '-> Category:', ncategory
    return ncategory


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


def correlate_external_pars(extpars, like, norm_value, targets_analyzed, npred, tb_ltf, tb_gbm):
    print '$$$$$$$$$$$'
    print extpars
    print '$$$$$$$$$$$'
    npred_sum = sum(npred.values())
    print 'Sum of target models: {0}'.format(npred_sum)
    par_values = {}
    mpred_temp = {}
    mpred = {}
    loglike = {}
    print 'Reading external values...'
    for target in targets_analyzed:
        print target
        tb_ltf_one = ReadLTFCatalogueInfo.select_one_by_name(tb_ltf, target)
        name_gbm = tb_ltf_one['GBM_assoc_key']
        tb_gbm_one = ReadGBMCatalogueInfo.select_one_by_name(tb_gbm, name_gbm)
        tb_lat = ReadLATCatalogueIndo.open_table()
        tb_lat_one = ReadLATCatalogueInfo.select_one_by_name(tb_gbm, target, tb_gbm)

        path_lat_result = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/Summary_{name}_afterglow_ver4.pickle'.format(name=target)
        lat_result = pickle_utilities.load(path_lat_result)
        if len(extpars)==1:
            extpar = extpars[0]
            if extpar in ('gbm_fluence', 'gbm_t90', 'gbm_flux1024', 'gbm_flux64', 'epeak_band'):
                par_values[target] = tb_lat_one['GBM'][dct_col_gbm_catarogue[extpar]]
            elif extpar in ('lat_flux_lowE'):
                try:
                    par_values[target] = lat_result['lower_energies']['flux']['value']
                except KeyError:
                    print 'KeyError!!!'
                    par_values[target] = 0
        elif len(extpars)>1:
            par_values[target] = 1
            for extpar in extpars:
           #if extpars in itertools.product(('gbm_fluence', 'gbm_t90', 'gbm_flux1024', 'gbm_flux64', 'epeak_band'), repeat=2):
                if extpar in ('gbm_fluence', 'gbm_t90', 'gbm_flux1024', 'gbm_flux64', 'epeak_band'):
                    vpar = tb_lat_one['GBM'][dct_col_gbm_catarogue[extpar]] * tb_lat_one['GBM'][dct_col_gbm_catarogue[extpar]]
                elif extpar=='lat_flux_lowE':
                    try:
                        vpar = lat_result['lower_energies']['flux']['value']
                    except KeyError:
                        print 'KeyError!!!'
                        vpar = 0
                else:
                    print 'Wrong external parameter!!'
                par_values[target] = par_values[target]*vpar
        mpred_temp[target] = par_values[target] / norm_value * npred[target]
    factor_npred = npred_sum/sum(mpred_temp.values())
#    npa_mpred = np.array()
    for target in targets_analyzed:
        norm_idx = like[target].par_index(target, 'Prefactor')
        like[target][norm_idx] = par_values[target] * factor_npred
        loglike[target] = like[target].logLike.value()
        mpred[target] = like[target].NpredValue(target)
        print '{0}: Prefactor={1}, loglike={2}, npred={3}'.format(target, like[target].params()[norm_idx].parameter.getTrueValue(), loglike[target], mpred[target])

    mpred_sum = sum(mpred.values())
    print 'Sum of target models: {0}'.format(mpred_sum)
    loglike_sum = sum(loglike.values())
    print 'Sum of loglikelihood: {0}'.format(loglike_sum)
    return loglike_sum


def run_composite2(chara_cat, path_input, path_outdir, names_params_tied_universal=['Prefactor', 'Index'], ncat_analyzed=0, str_suffix='', eedges=None, escale=1000., binned=False, longonly=False, tolerr=180., path_dct_category='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_LC.pickle', fig_forms=('png', 'pdf')):

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
    if isinstance(chara_cat, int):
        dct_chara_title[chara_cat] = 'Simulation No.{0:0>5}'.format(chara_cat)
    print dct_chara_title[chara_cat], 'Category No.', ncat_analyzed-1
    print lst_ncat_analyzed
    print len(lst_ncat_analyzed), 'GRBs.'

    nenergies = int(np.log10(eedges[1]/eedges[0])*4+0.5)
    energies = 10 ** np.linspace(np.log10(eedges[0]), np.log10(eedges[1]), nenergies+1)

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

        # tb_one = ReadLTFCatalogueInfo.select_one_by_name(tb, target)
        # gbm_flux = tb_one['FLUX_1024']
        # print 'GBM flux:', gbm_flux
        # norm_idx = like[target].par_index(target, 'Prefactor')
        # like[target].params()[norm_idx].parameter.setBounds(0, 100)
        # like[target].params()[norm_idx].parameter.setScale(gbm_flux*1e-12)
        # print 'Normalization scale is set to', like[target].params()[norm_idx].parameter.getScale()

        #escale_idx = like[target].par_index(target, 'Scale')
        #like[target].params()[escale_idx].parameter.setBounds(100, 100000)
        #like[target][escale_idx] = escale
        #like[target].freeze(escale_idx)

        for source in like[target].sourceNames():
            if source not in (target):
                like[target].normPar(source).setFree(False)
        sys.stdout.flush()

        CompositeLike.addComponent(like[target])
        sys.stdout.flush()
        #del like[target]
        gc.collect()
    print 'Analyzed GRBs:', targets_analyzed

    # for icat in range(NCATEGORIES_FLUENCE):
    #     if ncat_analyzed-1 not in (icat, -1):
    #         print 'skipped.'
    #         continue
    #     print '======================'
    #     print '===== Category', icat, '====='
    #     print '======================'
    #     print 'Target:', len(targets_analyzed), 'GRBs.'
    #     #for target in targets:
    #     #    if dct_category_fluence[target]==icat or ncat_analyzed==0:
    #     print lst_ncat_analyzed
    #     print ''

       # # Tying parameters for each fluence category separately
       #  tiedParams_category = {}
       #  for par in names_params_tied_category:
       #      tiedParams_category[par] = []
       #      for target in lst_ncat_analyzed:
       #          if target in targets:
       #              tiedParams_category[par].append(tuple([like[target], target, par]))
       #      CompositeLike.tieParameters(tiedParams_category[par])
       #  print '* Parameters tied by each category:'
       #  print tiedParams_category

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
    #print minuit.getRetCode()

    #fig_stacked, ax_stacked = plt.subplots(2, 1, figsize=(16, 10))
    x_stacked = (like[targets_analyzed[0]].energies[:-1] + like[targets_analyzed[0]].energies[1:])/2.
    print len(x_stacked), 'energy bins.'
    model_sum_stacked = np.zeros_like(like[targets_analyzed[0]]._srcCnts(like[targets_analyzed[0]].sourceNames()[0]))
    model_grb_stacked = np.zeros_like(model_sum_stacked)
    model_others_stacked = np.zeros_like(model_sum_stacked)
    nobs_sum_stacked = np.zeros_like(model_sum_stacked)
    npred = {}

    # Loop over GRBs
    for target in targets_analyzed: #lst_ncat_analyzed:
        print target
        #ncategory = dct_category_fluence[target]
        #print '  Producing plots...'
        sys.stdout.flush()
        #path_xml = '{0}/xml/likelihood_status_{1}{2}.xml'.format(path_outdir, target, str_suffix)
        #like[target].writeXml(path_xml)
        #path_spectra = '{0}/fits/counts_spectra_{1}{2}.fits'.format(path_outdir, target, str_suffix)
        #like[target].writeCountsSpectra(path_spectra)
        #fspec = fits.open(path_spectra)
        #tb_counts = fspec[1].data
        #tb_fluxes = fspec[2].data
        #tb_ebounds = fspec[3].data
        #fig, ax = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
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

        npred[target] = like[target].NpredValue(target)

        # for col_src in tb_counts.columns[1:]:
        #     model_sum = model_sum + tb_counts[col_src.name]
        #     model_sum_stacked = model_sum_stacked + tb_counts[col_src.name]
        #     if col_src.name == target:
        #         model_grb = model_grb + tb_counts[col_src.name]
        #         model_grb_stacked = model_grb_stacked + tb_counts[col_src.name]
        #     else:
        #         model_others = model_others + tb_counts[col_src.name]
        #         model_others_stacked = model_others_stacked + tb_counts[col_src.name]

        # nobs_sum = like[target]._Nobs() #tb_counts['ObsCounts']
        #lst_nobs_lat.append(sum(nobs_sum))
        #tb1 = ReadLTFCatalogueInfo.select_one_by_name(tb, target)
        #lst_fluence_gbm.append(tb1["FLUENCE"])
        #lst_fluence_gbm_err.append(tb1["FLUENCE_ERROR"])
        #rh_fluence_weightNobs.Fill(np.log10(lst_fluence_gbm[-1]), lst_nobs_lat[-1])

        # # Top3 in all categories
        # if lst_nobs_lat[-1]>sum(lst_tops[0]['nobs']):
        #     lst_tops[2] = lst_tops[1]
        #     lst_tops[1] = lst_tops[0]
        #     lst_tops[0] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}
        # elif lst_nobs_lat[-1]>sum(lst_tops[1]['nobs']):
        #     lst_tops[2] = lst_tops[1]
        #     lst_tops[1] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}
        # elif lst_nobs_lat[-1]>sum(lst_tops[2]['nobs']):
        #     lst_tops[2] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}                

    #     for ihiest in (1, 2, 3, 4):
    #         if nobs_sum[-ihiest]>0:
    #             print nobs_sum[-ihiest], 'events in the', ihiest, '-th highest energy bin.'
        #nobs_sum_stacked = nobs_sum_stacked + nobs_sum #tb_counts['ObsCounts']
    #     try:
    #         ax[0].loglog(x_stacked, model_sum, label='Sum of models')
    #         ax[0].loglog(x_stacked, model_grb, label=target)
    #         ax[0].loglog(x_stacked, model_others, label='Others')
    #         ax[0].errorbar(x_stacked, nobs_sum, yerr=np.sqrt(nobs_sum), fmt='o',label='Counts')
    #         ax[0].legend(loc=1, fontsize=20)
    #         ax[0].set_ylabel('[counts]')
    #         ax[0].set_title(target)
    #         ax[0].set_xticklabels([])
    #         ax[0].grid(ls='-', lw=0.5, alpha=0.2)
    #         resid = (nobs_sum - model_sum) / model_sum
    #         resid_err = np.sqrt(nobs_sum) / model_sum
    #         ax[1].set_xscale('log')
    #         ax[1].errorbar(x_stacked, resid, yerr=resid_err, fmt='o')
    #         ax[1].axhline(0.0,ls=':')
    #         ax[1].grid(ls='-', lw=0.5, alpha=0.2)
    #         ax[1].set_xlabel(r'$\log_{10}Energy \rm{[MeV]}$')
    #         ax[1].set_ylabel('Fractional residual')
    #         fig.tight_layout()
    #         fig.subplots_adjust(hspace=0)
    #         ax[1].set_yticks([y for y in ax[1].get_yticks() if y<ax[1].get_ylim()[1]])

    #         fig.savefig('{0}/plots/Spectrum{1}{2}.png'.format(path_outdir, target, str_suffix))
    #         plt.close()
            
    #     except ValueError:
    #         continue
        
    # # Histogram of GBM fluence
    # fig2d = plt.figure()
    # ax2d = fig2d.add_axes((0.1, 0.1, 0.8, 0.8))
    # npa_fluence_gbm = np.array(lst_fluence_gbm)
    # npa_fluence_gbm_err = np.array(lst_fluence_gbm_err)
    # npa_nobs_lat = np.array(lst_nobs_lat)
    # ax2d.set_xscale('log')
    # ax2d.set_yscale('log')
    # #ax2d.set_ylim(0.5, 200)
    # ax2d.errorbar(x=npa_fluence_gbm, y=npa_nobs_lat, xerr=npa_fluence_gbm_err, fmt='o')
    # #ax2d.errorbar(x=npa_fluence_gbm, y=npa_nobs_lat, xerr=npa_fluence_gbm_err, yerr=np.sqrt(npa_nobs_lat), fmt='o')
    # ax2d.axvline(FLUENCE_CUT[0],ls=':')
    # ax2d.axvline(FLUENCE_CUT[1],ls=':')
    # ax2d.set_xlabel('Fluence in GBM [erg/cm^{2}]')
    # ax2d.set_ylabel('Photons in LAT [counts]')
    # fig2d.savefig('{0}/plots/nobs_vs_GBMfluence{1}.png'.format(path_outdir, str_suffix))

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

       # E^2 x count plot
       # fig_E2stacked, ax_E2stacked = plt.subplots(2, 1, figsize=(16, 10))

       # fig_E2stacked.tight_layout()
       # fig_E2stacked.subplots_adjust(hspace=0)
       # ax_E2stacked[1].set_yticks([y for y in ax_E2stacked[1].get_yticks() if y<ax_E2stacked[1].get_ylim()[1]])

        for ff in fig_forms:
            fig_stacked.savefig('{0}/plots/StackedSpectrum{1}_category{2}.{3}'.format(path_outdir, str_suffix, ncat_analyzed, ff))

    dct_summary = {'dloglike_inv':fit_result, 'targets':targets_analyzed, 'parameters':dct_params, 'npred':npred }
    pickle_utilities.dump('{0}/Summary_StackedAnalysis{1}_category{2}.pickle'.format(path_outdir, str_suffix, ncat_analyzed), dct_summary)

    # print 'Check Correlation of parameters...'
    # loglike_sum = {}
    # for extpar in lst_extertnal_parameters:
    #     loglike_sum[extpar] = correlate_external_pars((extpar,), like, norm_value, targets_analyzed, npred, tb_ltf, tb_gbm)
    # for extpars in list(itertools.combinations_with_replacement(lst_extertnal_parameters, 2)):
    #     loglike_sum[str(extpars)] = correlate_external_pars(extpars, like, norm_value, targets_analyzed, npred, tb_ltf, tb_gbm)
#    for extpar in lst_extertnal_parameters:
    #for k, v in loglike_sum.items():
    #    print '{0}: dloglike={1}'.format(k, max(loglike_sum.values())-v)

    return (fit_result, dct_params, npred)


@click.command()
@click.argument('inputs', type=str)
@click.option('--pathout', '-o', type=str, default='.')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--tieuniv', multiple=True, type=str)
#@click.option('--tiecat', multiple=True, type=str)
@click.option('--category', type=click.Choice(['0', '1', '2', '3']), help='0: all GRBs 1,2,3: only GRBs of each category')
@click.option('--chara', type=click.Choice(['GBM_FLUENCE', 'GBM_FLUX_1024', 'GBM_FLUX_64', 'LC_INDEX']), help='GRB characterristics which categorize them.')
#@click.option('--chara', type=click.Choice(['gbm_fluence', 'gbm_t90', 'spec_index', 'flux_gev', 'lc_index', 'gbm_intermittent', 'gbm_flux1024', 'gbm_flux64', 'epeak_band', 'gbm_fluence_per_t50']), help='GRB characterristics which categorize them.')
@click.option('--energies', '-e', type=(float, float), default=(100., 100000.))
@click.option('--escale', type=float, default=1000.)
@click.option('--longonly', '-l', is_flag=True)
@click.option('--tolerr', '-t', type=float, default=180.)
@click.option('--binned', '-b', is_flag=True)
@click.option('--quantile', '-q', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_LC.pickle')
def main(chara, inputs, pathout, category, tieuniv, suffix, energies, escale, binned, longonly, tolerr, quantile):
        if category==None:
            categories = [0,1,2,3]
        else:
            categories = [int(category)]

        for category in categories:
            print '##### Category {0} #####'.format(category)
            run_composite2(chara, inputs, '{0}/Category{1}'.format(pathout, category), tieuniv, category, suffix, energies, escale, binned, longonly, tolerr, path_dct_category=quantile)
            print ''


if __name__ == '__main__':
    main()
