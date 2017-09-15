#!/usr/bin/env python

import sys
import os
import os.path
import matplotlib as mpl
mpl.use('tkagg')
import matplotlib.pyplot as plt
import gt_apps as my_apps
from pyLikelihood import *
from UnbinnedAnalysis import *
from Composite2 import Composite2
from CompositeLikelihood import CompositeLikelihood
from fermipy.utils import get_parameter_limits
#from bdlikeSED import *
import click
from astropy.io import fits
import numpy as np
#import ROOT
#ROOT.gROOT.SetBatch()
import ReadLTFCatalogueInfo
from pLsList import ls_list


def judge_category_fluence(tb, name, lst_cut):
    tb = ReadLTFCatalogueInfo.select_gbm_exist(tb)
    tb1 = ReadLTFCatalogueInfo.select_one_by_name(tb, name)
    ncategory = len(lst_cut)
    for ic in range(len(lst_cut)):
        ncategory -= int(tb1['FLUENCE']>=lst_cut[ic])
    print 'Fluence:', tb1['FLUENCE'], '-> Category:', ncategory
    return ncategory


def run_composite2(lst_inputs, path_outdir, names_params_tied_universal=['Index'], names_params_tied_category=['Prefactor'], ncategory=0, str_suffix=''):
    # Open table
    tb = ReadLTFCatalogueInfo.open_table()
    # Definition of GBM fluence categories
    FLUENCE_CUT = [1.09e-04, 3.06e-05] #[1.45E-4, 3.70E-5] # Top 10%, 35%
    NCATEGORIES_FLUENCE = len(FLUENCE_CUT)+1
    dct_category_fluence = {}
#    rh_fluence_weightNobs = ROOT.TH1D('roohtg', 'GBM Fluence', 100, -7, -2)

    path_base = os.getcwd()
    os.chdir(path_outdir)
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

    like={}
    targets = []
    lst_fluence_gbm = []
    lst_fluence_gbm_err = []
    lst_nobs_lat = []
    for (itarget, path_target) in enumerate(lst_inputs):
        path_base, name_base = os.path.split(path_target)
        target = name_base[3:12]
        targets.append(target)
        print '##### No.{0} {1} #####'.format(itarget, target)
        dct_category_fluence[target] = judge_category_fluence(tb, target, FLUENCE_CUT) 
        if ncategory-1 not in (dct_category_fluence[target], -1):
            print 'skipped.'
            continue

        ltcube = '/'.join((path_base, name_base+'_ft1_ltCube.fits'))
        expMap = '/'.join((path_base, name_base+'_ft1_expMap.fits'))
        srcModel = '/'.join((path_base, name_base+'_ft1_model.xml'))
        evt = '/'.join((path_base, name_base+'_ft1_filtered.fits'))
        sc = '/'.join((path_base, '../../../../..', name_base.replace('_P8_P302_BASE_T00-999-101000_r030', '_T00-999-101000_ft2-30s.fits')))
        if itarget==0:
            print 'Files of the first target.'
            print '  Event:', evt
            print '  Spacecraft:', sc
            print '  Livetime cube:', ltcube
            print '  Exposure map:', expMap
            print '  Source model:', srcModel

        # Diffuse responses
        my_apps.diffResps['evfile'] = evt
        my_apps.diffResps['scfile'] = sc
        my_apps.diffResps['srcmdl'] = srcModel
        my_apps.diffResps['irfs'] = irfs
        my_apps.diffResps.run()

        like[target] = unbinnedAnalysis(evfile=evt,
                                        scfile=sc,
                                        expmap=expMap,
                                        expcube=ltcube,
                                        irfs=irfs,
                                        srcmdl=srcModel,
                                        optimizer=optimizer)
        for source in like[target].sourceNames():
            if source not in (target):
                like[target].normPar(source).setFree(False)
        sys.stdout.flush()

    CompositeLike = []
    # Scan ecutoff
    xvals = 10 ** np.linspace(2.0, 5.125, 26)
    ecutoff_lim95 = {}
    ebreak_fixed = 10.
    for icat in range(NCATEGORIES_FLUENCE):
        print '======================'
        print '===== Category', icat, '====='
        print '======================'
        if ncategory-1 not in (icat, -1):
            print 'skipped.'
            continue
        fit_results = []
        for i, x in enumerate(xvals):
            print '---------------'
            print 'Cutoff energy:', x, 'MeV'
            tiedParams_category = {}
            tiedParams_universal = {}
            for par in names_params_tied_category:
                tiedParams_category[par] = []
            for par in names_params_tied_universal:
                tiedParams_universal[par] = []
            CompositeLike = Composite2(optimizer=optimizer)
            for target in targets:
                if dct_category_fluence[target]==icat:
                    print target
                  # Fixing cutoff energy
                    ecutoff_index = like[target].par_index(target, 'P1')
                    print like[target][ecutoff_index]
                    like[target][ecutoff_index] = x
                    like[target].freeze(ecutoff_index)
                    print "  {0}'s P1 is freezed as {1}".format(target, like[target][ecutoff_index])
                  # Fixing break energy
                    ebreak_index = like[target].par_index(target, 'Ebreak')
                    print like[target][ebreak_index]
                    like[target][ebreak_index] = ebreak_fixed
                    like[target].freeze(ebreak_index)
                    print "  {0}'s Ebreak is freezed as {1}".format(target, like[target][ebreak_index])
                    CompositeLike.addComponent(like[target])
                   # Tying parameters for each fluence category separately
                    for par in names_params_tied_category:
                        tiedParams_category[par].append(tuple([like[target], target, par]))
                   # Tying parameters universaly
                    for par in names_params_tied_universal:
                        tiedParams_universal[par].append(tuple([like[target], target, par]))
                    sys.stdout.flush()
            print '* Parameters tied by each category:'
            print tiedParams_category
            print '* Parameters tied universaly:'
            print tiedParams_universal
            for par in names_params_tied_category:
                CompositeLike.tieParameters(tuple(tiedParams_category[par]))
            for par in names_params_tied_universal:
                CompositeLike.tieParameters(tuple(tiedParams_universal[par]))

            fit_results.append(CompositeLike.fit(covar=False,tol=1.e-2,optimizer=optimizer))

        # Limit of ecutoff
        loglike_inversed_scanned = np.array(fit_results)
        loglike_inversed_min = min(loglike_inversed_scanned)
        print '* Negative log-likelihood :'
        for (i,s) in enumerate(loglike_inversed_scanned):
            print '{0} ({1}) at {2} MeV'.format(s, s-loglike_inversed_min, xvals[i])
        ecutoff_lim95[icat] = get_parameter_limits(xvals, -1*loglike_inversed_scanned)
    print '* 95% limit of cutoff energy:'
    print ecutoff_lim95


    # fig_stacked, ax_stacked = plt.subplots(2, 2, figsize=(16, 10))
    # x_stacked = (like[targets[0]].energies[:-1] + like[targets[0]].energies[1:])/2.
    # print len(x_stacked), 'energy bins.'
    # model_sum_stacked = np.zeros_like(like[targets[0]]._srcCnts(like[targets[0]].sourceNames()[0]))
    # model_grb_stacked = np.zeros_like(model_sum_stacked)
    # model_others_stacked = np.zeros_like(model_sum_stacked)
    # nobs_sum_stacked = np.zeros_like(model_sum_stacked)

    # # Objects for plotting three GRB categories based on GBM fluence
    # lst_fig_stacked_subs = []
    # lst_ax_stacked_subs = []
    # lst_nobs_sum_stacked_sub = []
    # lst_model_sum_stacked_sub = []
    # lst_model_grb_stacked_sub = []
    # lst_model_others_stacked_sub = []
    # lst_tops = [{'name':'', 'fluence':0, 'nobs':np.zeros_like(model_sum_stacked)} for x in (1, 2, 3)]
    # lst_subtops = []
    # for isub in range(NCATEGORIES_FLUENCE):
    #     subpl = plt.subplots(2, 2, figsize=(16, 10))
    #     lst_fig_stacked_subs.append(subpl[0])
    #     lst_ax_stacked_subs.append(subpl[1])
    #     lst_nobs_sum_stacked_sub.append(np.zeros_like(model_sum_stacked))
    #     lst_model_sum_stacked_sub.append(np.zeros_like(model_sum_stacked))
    #     lst_model_grb_stacked_sub.append(np.zeros_like(model_sum_stacked))
    #     lst_model_others_stacked_sub.append(np.zeros_like(model_sum_stacked))
    #     # Top3 in each category
    #     lst_subtops.append([{'name':'', 'fluence':0, 'nobs':np.zeros_like(model_sum_stacked)} for x in (1, 2, 3)])

    # # Loop over GRBs
    # for target in targets:
    #     print target
    #     ncategory = dct_category_fluence[target]
    #     print '  Producing plots...'
    #     sys.stdout.flush()
    #     path_xml = '{0}/xml/likelihood_status_{1}{2}.xml'.format(path_outdir, target, str_suffix)
    #     like[target].writeXml(path_xml)
    #     path_spectra = '{0}/fits/counts_spectra_{1}{2}.fits'.format(path_outdir, target, str_suffix)
    #     like[target].writeCountsSpectra(path_spectra)
    #     fspec = fits.open(path_spectra)
    #     tb_counts = fspec[1].data
    #     #tb_fluxes = fspec[2].data
    #     #tb_ebounds = fspec[3].data
    #     fig, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    #     model_sum = np.zeros_like(model_sum_stacked)
    #     model_grb = np.zeros_like(model_sum_stacked)
    #     model_others = np.zeros_like(model_sum_stacked)
    #     nobs_sum = np.zeros_like(model_sum)
    #     for col_src in tb_counts.columns[1:]:
    #         model_sum = model_sum + tb_counts[col_src.name]
    #         model_sum_stacked = model_sum_stacked + tb_counts[col_src.name]
    #         lst_model_sum_stacked_sub[ncategory] = lst_model_sum_stacked_sub[ncategory] + tb_counts[col_src.name]
    #         if col_src.name == target:
    #             model_grb = model_grb + tb_counts[col_src.name]
    #             model_grb_stacked = model_grb_stacked + tb_counts[col_src.name]
    #             lst_model_grb_stacked_sub[ncategory] = lst_model_grb_stacked_sub[ncategory] + tb_counts[col_src.name]
    #         else:
    #             model_others = model_others + tb_counts[col_src.name]
    #             model_others_stacked = model_others_stacked + tb_counts[col_src.name]
    #             lst_model_others_stacked_sub[ncategory] = lst_model_others_stacked_sub[ncategory] + tb_counts[col_src.name]

    #     nobs_sum = tb_counts['ObsCounts']
    #     lst_nobs_lat.append(sum(nobs_sum))
    #     tb1 = ReadLTFCatalogueInfo.select_one_by_name(tb, target[3:])
    #     lst_fluence_gbm.append(tb1["FLUENCE"])
    #     lst_fluence_gbm_err.append(tb1["FLUENCE_ERROR"])
    #     #rh_fluence_weightNobs.Fill(np.log10(lst_fluence_gbm[-1]), lst_nobs_lat[-1])

    #     # Top3 in all categories
    #     if lst_nobs_lat[-1]>sum(lst_tops[0]['nobs']):
    #         lst_tops[2] = lst_tops[1]
    #         lst_tops[1] = lst_tops[0]
    #         lst_tops[0] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}
    #     elif lst_nobs_lat[-1]>sum(lst_tops[1]['nobs']):
    #         lst_tops[2] = lst_tops[1]
    #         lst_tops[1] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}
    #     elif lst_nobs_lat[-1]>sum(lst_tops[2]['nobs']):
    #         lst_tops[2] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}                

    #     # Top3 in each category
    #     if lst_nobs_lat[-1]>sum(lst_subtops[ncategory][0]['nobs']):
    #         lst_subtops[ncategory][2] = lst_subtops[ncategory][1]
    #         lst_subtops[ncategory][1] = lst_subtops[ncategory][0]
    #         lst_subtops[ncategory][0] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}
    #     elif lst_nobs_lat[-1]>sum(lst_subtops[ncategory][1]['nobs']):
    #         lst_subtops[ncategory][2] = lst_subtops[ncategory][1]
    #         lst_subtops[ncategory][1] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}
    #     elif lst_nobs_lat[-1]>sum(lst_subtops[ncategory][2]['nobs']):
    #         lst_subtops[ncategory][2] = {'name':target, 'fluence':tb1["FLUENCE"], 'nobs':nobs_sum}        

    #     for ihiest in (1, 2, 3, 4):
    #         if nobs_sum[-ihiest]>0:
    #             print nobs_sum[-ihiest], 'events in the', ihiest, '-th highest energy bin.'
    #     nobs_sum_stacked = nobs_sum_stacked + tb_counts['ObsCounts']
    #     lst_nobs_sum_stacked_sub[ncategory] = lst_nobs_sum_stacked_sub[ncategory] + tb_counts['ObsCounts']
    #     ax[0].loglog(x_stacked, model_sum, label='Sum of models')
    #     ax[0].loglog(x_stacked, model_grb, label=target)
    #     ax[0].loglog(x_stacked, model_others, label='Others')
    #     ax[0].errorbar(x_stacked, nobs_sum, yerr=np.sqrt(nobs_sum), fmt='o',label='Counts')
    #     ax[0].legend(loc=1, fontsize=12)
    #     resid = (nobs_sum - model_sum) / model_sum
    #     resid_err = np.sqrt(nobs_sum) / model_sum
    #     ax[1].set_xscale('log')
    #     ax[1].errorbar(x_stacked, resid, yerr=resid_err, fmt='o')
    #     ax[1].axhline(0.0,ls=':')
    #     fig.savefig('{0}/plots/Spectrum{1}{2}.png'.format(path_outdir, target, str_suffix))
    #     plt.close()
        
    # #rh_fluence_weightNobs.SaveAs('htg_GBM-fluence_weightNobs.root')
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

    # # Count spectrum
    # ax_stacked[0,0].loglog(x_stacked, model_sum_stacked, label='Sum of models')
    # ax_stacked[0,0].loglog(x_stacked, model_grb_stacked, label='GRBs')
    # ax_stacked[0,0].loglog(x_stacked, model_others_stacked, label='Others')
    # ax_stacked[0,0].errorbar(x_stacked, nobs_sum_stacked, yerr=np.sqrt(nobs_sum_stacked), fmt='o',label='Counts')
    # ax_stacked[0,0].legend(loc=1, fontsize=12)
    # ax_stacked[0,0].set_ylabel('[counts]')
    # ax_stacked[0,0].set_xlabel(r'$\log_{10}Energy$ [MeV]')

    # resid_stacked = (nobs_sum_stacked - model_sum_stacked) / model_sum_stacked
    # resid_stacked_err = np.sqrt(nobs_sum_stacked) / model_sum_stacked
    # ax_stacked[0,1].set_xscale('log')
    # ax_stacked[0,1].errorbar(x_stacked, resid_stacked, yerr=resid_stacked_err, fmt='o')
    # ax_stacked[0,1].axhline(0.0,ls=':')
    # #ax_stacked[0].set_xlabel(r'$\log{10}Energy$ [MeV]')
    # ax_stacked[0,1].set_xlabel(r'$\log_{10}Energy$ [MeV]')
    # ax_stacked[0,1].set_ylabel('Fractional residual')

    # nobs_denominator = np.zeros_like(nobs_sum_stacked)
    # for (inobs_den, nobs_den) in enumerate(nobs_sum_stacked):
    #     nobs_denominator[inobs_den] = max(1, nobs_den)

    # ax_stacked[1,0].stackplot(x_stacked, lst_tops[0]['nobs']/nobs_denominator, lst_tops[1]['nobs']/nobs_denominator, lst_tops[2]['nobs']/nobs_denominator, labels=['No.1 '+lst_tops[0]['name'], 'No.2 '+lst_tops[1]['name'], 'No.3 '+lst_tops[2]['name']])
    # ax_stacked[1,0].legend(loc=2, fontsize=12, fancybox=True, framealpha=0.5)
    # ax_stacked[1,0].set_xscale('log')
    # ax_stacked[1,0].set_xlabel(r'$\log_{10}Energy$ [MeV]')
    # ax_stacked[1,0].set_ylabel('Occupation rate')
    # fig_stacked.savefig('{0}/plots/StackedSpectrum{1}.png'.format(path_outdir, str_suffix))

    # for isub in range(NCATEGORIES_FLUENCE):
    #     lst_ax_stacked_subs[isub][0,0].loglog(x_stacked, lst_model_sum_stacked_sub[isub], label='Sum of models')
    #     lst_ax_stacked_subs[isub][0,0].loglog(x_stacked, lst_model_grb_stacked_sub[isub], label='GRBs')
    #     lst_ax_stacked_subs[isub][0,0].loglog(x_stacked, lst_model_others_stacked_sub[isub], label='Others')
    #     lst_ax_stacked_subs[isub][0,0].errorbar(x_stacked, lst_nobs_sum_stacked_sub[isub], yerr=np.sqrt(lst_nobs_sum_stacked_sub[isub]), fmt='o',label='Counts')
    #     lst_ax_stacked_subs[isub][0,0].legend(loc=0, fontsize=12)
    #     lst_ax_stacked_subs[isub][0,0].set_xlabel(r'$\log_{10}Energy$ [MeV]')
    #     lst_ax_stacked_subs[isub][0,0].set_ylabel('[counts]')

    #     resid_stacked_sub = (lst_nobs_sum_stacked_sub[isub] - lst_model_sum_stacked_sub[isub]) / lst_model_sum_stacked_sub[isub]
    #     resid_stacked_sub_err = np.sqrt(lst_nobs_sum_stacked_sub[isub]) / lst_model_sum_stacked_sub[isub]
    #     lst_ax_stacked_subs[isub][0,1].set_xscale('log')
    #     lst_ax_stacked_subs[isub][0,1].errorbar(x_stacked, resid_stacked_sub, yerr=resid_stacked_sub_err, fmt='o')
    #     lst_ax_stacked_subs[isub][0,1].axhline(0.0,ls=':')
    #     #lst_ax_stacked_subs[isub][0,0].set_xlabel(r'$\log{10}Energy$ [MeV]')
    #     lst_ax_stacked_subs[isub][0,1].set_xlabel(r'$\log_{10}Energy$ [MeV]')
    #     lst_ax_stacked_subs[isub][0,1].set_ylabel('Fractional residual')

    #     nobs_denominator = np.zeros_like(lst_nobs_sum_stacked_sub[isub])
    #     for (inobs_den, nobs_den) in enumerate(lst_nobs_sum_stacked_sub[isub]):
    #         nobs_denominator[inobs_den] = max(1, nobs_den)

    #     lst_ax_stacked_subs[isub][1,0].set_xscale('log')
    #     lst_ax_stacked_subs[isub][1,0].set_ylim(0, 1)
    #     lst_ax_stacked_subs[isub][1,0].stackplot(x_stacked, lst_subtops[isub][0]['nobs']/nobs_denominator, lst_subtops[isub][1]['nobs']/nobs_denominator, lst_subtops[isub][2]['nobs']/nobs_denominator, labels=['No.1 '+lst_subtops[isub][0]['name'], 'No.2 '+lst_subtops[isub][1]['name'], 'No.3 '+lst_subtops[isub][2]['name']])
    #     lst_ax_stacked_subs[isub][1,0].legend(loc=2, fontsize=12, fancybox=True, framealpha=0.5)
    #     #lst_ax_stacked_subs[isub][1,0].set_xlabel(r'$\log{10}Energy$ [MeV]')
    #     lst_ax_stacked_subs[isub][1,0].set_xlabel(r'$\log_{10}Energy$ [MeV]')
    #     lst_ax_stacked_subs[isub][1,0].set_ylabel('Occupation rate')

    #     lst_fig_stacked_subs[isub].savefig('{0}/plots/StackedSpectrum_category{1}{2}.png'.format(path_outdir, isub, str_suffix))


@click.command()
@click.argument('inputs', type=str)
@click.option('--pathout', '-o', type=str, default='.')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--tieuniv', multiple=True, type=str)
@click.option('--tiecat', multiple=True, type=str)
@click.option('--category', type=click.Choice(['0', '1', '2', '3']), help='0: all GRBs 1,2,3: only GRBs of each category')
def main(inputs, pathout, tieuniv, tiecat, category, suffix):
    with open(inputs, "r") as filein:
        str_paths = filein.read()
        input_paths = str_paths.split('\n')[:-1]
        print input_paths
        run_composite2(input_paths, pathout, tieuniv, tiecat, int(category), suffix)


if __name__ == '__main__':
    main()
