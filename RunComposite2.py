#!/usr/bin/env python

import sys
import os
import os.path
import subprocess
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


DCT_EDGE_ENERGIES = {'wholeE':(100, 100000), 'lowE':(100, 1000), 'highE':(10000, 100000), 'midE':(1000, 10000)}
DCT_DIR_ENERGIES = {'wholeE':'E0000100-0100000MeV', 'lowE':'E0000100-0001000MeV', 'highE':'E0010000-0100000MeV', 'midE':'E0001000-0010000MeV'}
DCT_DIR_ROI = {'wholeE':'r12deg', 'lowE':'r12deg', 'highE':'r01deg', 'midE':'r03deg'}
DCT_NORM_NAME = {'PowerLaw':'Prefactor', 'PowerLaw2':'Integral', 'ScaleFactor::PowerLaw2':'Integral'}
#DCT_EDGE_ENERGIES = {'wholeE':(100, 100000), 'lowE':(100, 10**3.75), 'highE':(10000, 100000)}
#DCT_DIR_ENERGIES = {'wholeE':'E0000100-0100000MeV', 'lowE':'E0000100-0005623MeV', 'highE':'E0010000-0100000MeV'}
#DCT_DIR_ROI = {'wholeE':'r12deg', 'lowE':'r12deg', 'highE':'r01deg'}
#DCT_NORM_NAME = {'PowerLaw':'Prefactor', 'PowerLaw2':'Integral', 'ScaleFactor::PowerLaw2':'Integral'}


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


def run_composite2(path_input, path_outdir, names_params_tied_universal=['Prefactor', 'Index'], str_suffix='', eedges=None, binned=False, longonly=False, tolerr=180., path_dct_category='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_LC.pickle', fig_forms=('png', 'pdf'), norm_name = 'Integral', exclude=None):

    # Open table
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    tb_lat = ReadLATCatalogueInfo.open_table()
    lst_lat = ReadLATCatalogueInfo.read_all(tb_lat, tb_gbm)
    lst_lat = ReadLATCatalogueInfo.remove_by_name(lst_lat, exclude)
    #lst_lat = ReadLATCatalogueInfo.select_by_name(tb_lat, '080000000', '100000000', tb_gbm)
    if longonly==True:
        lst_lat = ReadLATCatalogueInfo.select_long(lst_lat)
    lst_lat = ReadLATCatalogueInfo.select_small_error(lst_lat, tolerr)
    lst_lat_name = []
    for g in lst_lat:
        lst_lat_name.append(g['GRBNAME'])

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
    like={}
    targets = []
    grbs_skipped = ('130427324', '150314205') #, '141022087'

    dct_chara_title = {'GBM_FLUENCE':'Fluence in GBM (prompt)', 
                       'GBM_FLUX_1024':'Peak flux in GBM for 1024ms (prompt)', 
                       'GBM_FLUX_64':'Peak flux in GBM for 64ms (prompt)',
                       #'epeak_band':'Peak energy of prompt Band component', 
                       #'gbm_intermittent':'GBM intermittent duration',
                       #'spec_index':'Spectral index in 0.178 - 5.62 GeV', 
                       #'flux_gev':'Flux in 0.178 - 5.62 GeV', 
                       'LC_INDEX':'Light curve index in 0.1 - 100 GeV'
                       }
    dct_categories = {'ALL': lst_lat_name,
                      'GBM_FLUENCE_1':get_category('GBM_FLUENCE', 1, path_dct_category),
                      'GBM_FLUENCE_2':get_category('GBM_FLUENCE', 2, path_dct_category),
                      'GBM_FLUENCE_3':get_category('GBM_FLUENCE', 3, path_dct_category),
                      'GBM_FLUX_1024_1':get_category('GBM_FLUX_1024', 1, path_dct_category),
                      'GBM_FLUX_1024_2':get_category('GBM_FLUX_1024', 2, path_dct_category),
                      'GBM_FLUX_1024_3':get_category('GBM_FLUX_1024', 3, path_dct_category),
                      'LC_INDEX_1':get_category('LC_INDEX', 1, path_dct_category),
                      'LC_INDEX_2':get_category('LC_INDEX', 2, path_dct_category)                      
                      }

    #lst_ncat_analyzed = get_category(chara_cat, ncat_analyzed, path_dct_category)
    #if isinstance(chara_cat, int):
    #    dct_chara_title[chara_cat] = 'Simulation No.{0:0>5}'.format(chara_cat)
    #print dct_chara_title[chara_cat], 'Category No.', ncat_analyzed-1
    #print lst_ncat_analyzed
    #print len(lst_ncat_analyzed), 'GRBs.'

    nenergies = int(np.log10(eedges[1]/eedges[0])*4+0.5)
    energies = 10 ** np.linspace(np.log10(eedges[0]), np.log10(eedges[1]), nenergies+1)

    srcs_virtual = {} # For calculating sum of LAT fluence

    for (itarget, target_info) in enumerate(lst_lat):
        target = target_info['GRBNAME'] 
        path_target = '{0}/{1}/GRB{0}_P8_P302_BASE_T00-999-101000_r030'.format(target, path_input)
        path_base, name_base = os.path.split(path_target)
        print '##### No.{0} {1} #####'.format(itarget, target)
        #if target not in lst_ncat_analyzed:
        #    print 'skipped.'
        #    continue

        #targets_analyzed.append(target)
        evt = '/'.join((path_base, name_base+'_ft1_filtered_gti.fits'))
        sc = '/'.join((path_base, '../../../../..', name_base.replace('_P8_P302_BASE_T00-999-101000_r030', '_T00-999-101000_ft2-30s.fits')))
        ltcube = '/'.join((path_base, name_base+'_ft1_ltCube.fits'))
        expMap = '/'.join((path_base, name_base+'_ft1_expMap.fits'))
        srcMap = '/'.join((path_base, name_base+'_ft1_srcmap.fits'))
        ccube = '/'.join((path_base, name_base+'_ft1_ccube.fits'))
        #srcModel = '/'.join((path_base, name_base+'_ft1_model.xml'))
        srcModel = '/'.join((path_base, name_base+'_ft1_model_new.xml'))

        if itarget==0:
            print 'Files of the first target.'
            print '  Event:', evt
            print '  Spacecraft:', sc
            print '  Livetime cube:', ltcube
            print '  Exposure map:', expMap
            print '  Source model:', srcModel
        srcs_virtual[target] = {}
        srcs_virtual[target]['time'] = GetGTI.get_duration(evt)
        if target in grbs_skipped:
            print 'Skipping {0}...'.format(target)
            continue
        if srcs_virtual[target]['time']['conservative']<=0: 
            print 'Valid observation time: {0}'.format(srcs_virtual[target]['time']['conservative'])
            print 'Skipping {0}...'.format(target)
            continue
        targets.append(target)
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
                #if energies is None:
                #    energies = (100., 100000.)
                #like[target].reset_ebounds(energies)
                #print 'Energy bound has reset to {0}'.format(like[target].energies)
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

    dct_summary = {} #'dloglike_inv':fit_result, 'targets':lst_cat, 'parameters':dct_params, 'srcs_virtual':srcs_virtual }

    for key_cat, lst_cat in dct_categories.items():
        print '=====', key_cat, '====='
        dct_summary[key_cat] = {}
        print 'Analyzed GRBs:', lst_cat
        dct_summary[key_cat]['targets'] = lst_cat
        CompositeLike = Composite2(optimizer=optimizer)
        for tgt_cat in lst_cat:
            if tgt_cat in grbs_skipped:
                print 'Skipping {0}...'.format(tgt_cat)
                continue
            if srcs_virtual[tgt_cat]['time']['conservative']<=0:
                print 'Valid observation time: {0}'.format(srcs_virtual[tgt_cat]['time']['conservative'])
                print 'Skipping {0}...'.format(tgt_cat)
                continue
            CompositeLike.addComponent(like[tgt_cat])
        #sys.stdout.flush()
        #gc.collect()

        # Tying parameters universaly
        print '* Parameters tied universaly:'
        print names_params_tied_universal #tiedParams_universal
        tiedParams_universal = {}
        for par in names_params_tied_universal:
            tiedParams_universal[par] = []
            for target in lst_cat: #lst_ncat_analyzed:
                if target in grbs_skipped:
                    continue
                if srcs_virtual[target]['time']['conservative']<=0:
                    continue
                tiedParams_universal[par].append(tuple([like[target], target, par]))
            CompositeLike.tieParameters(tiedParams_universal[par])

       #minuit = eval("pyLike.%s(CompLike.composite)"%optimizer)
       #minuit.setStrategy(2)
       #likeobj = pyLike.NewMinuit(like.logLike)
        fit_result = CompositeLike.fit(covar=True,tol=1.e-5,optimizer=optimizer)
        print '== Fitting result =='
        print fit_result
        dct_summary[key_cat]['dloglike_inv'] = fit_result

        dct_params = {}
        for tiedpar in names_params_tied_universal:
            dct_params[tiedpar] = {}
            dct_params[tiedpar]['value'] = like[lst_cat[0]].model[lst_cat[0]].funcs['Spectrum'].getParam(tiedpar).value()
            dct_params[tiedpar]['error'] = like[lst_cat[0]].model[lst_cat[0]].funcs['Spectrum'].getParam(tiedpar).error()
            print '{n}: {v} +/- {e}'.format(n=tiedpar, v=dct_params[tiedpar]['value'], e=dct_params[tiedpar]['error'])
        dct_summary[key_cat]['parameters'] = dct_params
        print ''

        x_stacked = (like[lst_cat[0]].energies[:-1] + like[lst_cat[0]].energies[1:])/2.
        print len(x_stacked), 'energy bins.'
        model_sum_stacked = np.zeros_like(like[lst_cat[0]]._srcCnts(like[lst_cat[0]].sourceNames()[0]))
        model_grb_stacked = np.zeros_like(model_sum_stacked)
        model_others_stacked = np.zeros_like(model_sum_stacked)
        nobs_sum_stacked = np.zeros_like(model_sum_stacked)

       # Loop over GRBs
        for target in lst_cat:
            print target
            if target in grbs_skipped:
                print 'Skipping {0}...'.format(target)
                continue
            if srcs_virtual[target]['time']['conservative']<=0:
                print 'Skipping {0}...'.format(target)
                continue
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

            srcs_virtual[target]['npred'] = like[target].NpredValue(target)
            srcs_virtual[target]['flux'] = {}
            srcs_virtual[target]['flux']['value'] = like[target].flux(target)
            srcs_virtual[target]['flux']['error'] = srcs_virtual[target]['flux']['value']*dct_params[norm_name]['error']/dct_params[norm_name]['value'] if dct_params[norm_name]['value']>0 else 0 #like[target].fluxError(target) }
            srcs_virtual[target]['eflux'] = {}
            srcs_virtual[target]['eflux']['value'] = like[target].energyFlux(target)
            srcs_virtual[target]['eflux']['error'] = srcs_virtual[target]['eflux']['value']*dct_params[norm_name]['error']/dct_params[norm_name]['value'] if dct_params[norm_name]['value']>0 else 0 #like[target].fluxError(target) }

            srcs_virtual[target]['fluence'] = {}
            for t in ('nominal', 'conservative'):
                srcs_virtual[target]['fluence'][t] = {}
                for l in ('value', 'error'):
                    srcs_virtual[target]['fluence'][t][l] = srcs_virtual[target]['time'][t]*srcs_virtual[target]['eflux'][l]
        dct_summary[key_cat]['srcs_virtual'] = srcs_virtual

        # Count spectrum 
        if fig_forms is not None:
            fig_stacked, ax_stacked = plt.subplots(2, 1, figsize=(16, 10))
            print 'Plotting count spectra...'
            print 'X-axis:', x_stacked
            print 'Y-axis:', model_sum_stacked
            print 'Y-axis:', model_grb_stacked
            print 'Y-axis:', model_others_stacked
            print 'Y-axis:', nobs_sum_stacked
            ax_stacked[0].set_ylim(0.02, 2E4)
            ax_stacked[0].loglog(x_stacked, model_sum_stacked, label='Sum of models')
            ax_stacked[0].loglog(x_stacked, model_grb_stacked, label='GRBs')
            ax_stacked[0].loglog(x_stacked, model_others_stacked, label='Others')
            ax_stacked[0].errorbar(x_stacked, nobs_sum_stacked, yerr=np.sqrt(nobs_sum_stacked), fmt='o',label='Counts')
            ax_stacked[0].legend(loc=1, fontsize=20)
            ax_stacked[0].set_title(key_cat) #'{chara} Category No.{cat}'.format(chara=dct_chara_title[chara_cat], cat=ncat_analyzed))
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
                fig_stacked.savefig('{0}/plots/StackedSpectrum_{1}{2}.{3}'.format(path_outdir, key_cat, str_suffix, ff))
        

        pickle_utilities.dump('{0}/Summary_StackedAnalysis{1}{2}.pickle'.format(path_outdir, key_cat, str_suffix), dct_summary)

    return dct_summary #(fit_result, dct_params, npred)


@click.command()
#@click.argument('inputs', type=str)
@click.option('--pathout', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/Stacking')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--tieuniv', multiple=True, type=str)
@click.option('--energies', '-e', multiple=True, type=str)#type=(float, float), default=(100., 100000.))
@click.option('--phases', '-p', multiple=True, type=str)
@click.option('--longonly', '-l', is_flag=True)
@click.option('--tolerr', '-t', type=float, default=180.)
@click.option('--binned', is_flag=True)
@click.option('--quantile', '-q', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_LC.pickle')
@click.option('--func', '-f', type=str, default='PowerLaw2')
@click.option('--exclude', multiple=True, type=str, default=('130427324','150314205'))
@click.option('--bsub', '-b', is_flag=True)
def main(pathout, tieuniv, suffix, energies, phases, binned, longonly, tolerr, quantile, bsub, func, exclude):
    for e, p in itertools.product(energies, phases):
        path_input = '/'.join([DCT_DIR_ENERGIES[e], DCT_DIR_ROI[e], p, func, 'Binned' if binned==True else 'Unbinned'])
        path_output = '/'.join([pathout, p, e])
        if not os.path.exists(path_output):
            os.mkdir(path_output)
        if bsub==False:
            run_composite2(path_input, path_output, tieuniv, suffix, DCT_EDGE_ENERGIES[e], binned, longonly, tolerr, path_dct_category=quantile, norm_name=DCT_NORM_NAME[func], exclude=exclude)
        else:
            acmd = ['bsub', '-o','{0}/StackingGRBs_{1}_{2}{3}.log'.format(path_output, e, p, suffix if suffix=='' else '_'+suffix), '-J','RC{0}{1}'.format(p,e), '-W','900', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/RunComposite2.py', '-o', pathout, '-s', '{0}'.format(suffix), '--tieuniv', DCT_NORM_NAME[func], '--tieuniv', 'Index', '-e', e, '-p', p, '-t', str(tolerr), '-q', quantile, '--func', func]
            if longonly==True:
                acmd.append('-l')
            if binned==True:
                acmd.append('--binned')
            if exclude is not None and len(exclude)>0:
                for exc in exclude:
                    acmd.append('--exclude')
                    acmd.append(exc)
            print acmd
            subprocess.call(acmd)


if __name__ == '__main__':
    main()
