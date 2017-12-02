#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd
import click
import ReadLATCatalogueInfo
#import ReadLTFCatalogueInfo
import ReadGBMCatalogueInfo
import FindGoodstatPeriods
import pickle_utilities
import pIntegratePowerlawLightcurve

##### PATH of Catalogue #####
GRB_CATALOGUE_LTF = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'


def quantile_grbs(output, longonly, tolerr, suffix):

    str_output = '{lo}{te}{suf}'.format(lo='_longonly' if longonly==True else '', te='tolerr'.format(int(tolerr*10)) if tolerr<180 else '', suf='_'+suffix if suffix!='' else suffix)
    
    tb_lat = ReadLATCatalogueInfo.open_table()
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    lst_lat = ReadLATCatalogueInfo.read_all(tb_lat, tb_gbm)
    lst_lat = ReadLATCatalogueInfo.select_gbm_exist(lst_lat)
    if longonly==True:
        lst_lat = ReadLATCatalogueInfo.select_long(lst_lat)
    lst_lat = ReadLATCatalogueInfo.select_small_error(lst_lat, tolerr)
    grbs_skipped = ('130427324', '150314205') #, '141022087'
    ngrb_lat = len(lst_lat)-len(grbs_skipped)
    print 'Number of GRBs:', ngrb_lat

    charas = ['GRBNAME', 'GBM_NAME', 'GBM_T90', 'GBM_T90_START', 'GBM_FLUENCE', 'GBM_FLUX_1024', 'GBM_FLUX_64', 'RA', 'DEC', 'NCOUNTS', 'LC_INDEX']
    arrays = {}
    for chara in charas:
        if chara in ('GRBNAME', 'GBM_NAME'):
            arrays[chara] = np.chararray(ngrb_lat, itemsize=9 if chara=='GRBNAME' else 12)
        else:
            arrays[chara] = np.zeros(ngrb_lat)
    #spec_index = np.zeros(ngrb_lat)
    #flux_gev = np.zeros(ngrb_lat)
    #lc_index = np.zeros(ngrb_lat)

    nevent = 0
    nevent_validlc = 0
    lc_phases = ('prompt', 'T95to03ks', '03ksto100ks')
    integrated_lc = {}
    scaled_lc = {}
    lc_indices = (-1.0, -1.3)
    jgrb = 0
    for igrb, grb_lat in enumerate(lst_lat):
        print 'No.{0} {1}'.format(igrb, grb_lat['GRBNAME'])
        if grb_lat['GRBNAME'] in grbs_skipped:
            print 'Skipped...'
            continue
        path_lc_fit = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/E0000100-0100000MeV/r12deg/lightcurve/LightCurve_{name}_indexfreeLAT_N10in2deg_fit.pickle'.format(name=grb_lat['GRBNAME'])
        lc_fit = pickle_utilities.load(path_lc_fit)
        path_gtifile = '{name}/E0000100-0100000MeV/r12deg/unified/PowerLaw2/Unbinned/GRB{name}_P8_P302_BASE_T00-999-101000_r030_ft1_filtered_gti.fits'.format(name=grb_lat['GRBNAME'])
        integrated_lc[grb_lat['GRBNAME']] = pIntegratePowerlawLightcurve.integrate_time_phases(grb_lat['GRBNAME'], path_gtifile, indices=lc_indices, tpl_phases=lc_phases)

        for chara in charas:
            if chara=='NCOUNTS':
                arrays[chara][jgrb] = FindGoodstatPeriods.get_entries_roi(path_gtifile, None, None, rlim=5.0, ra=grb_lat['RA'], dec=grb_lat['DEC'])
            elif chara[:4]=='GBM_':
                arrays[chara][jgrb] = grb_lat['GBM'][chara[4:]]
            elif chara=='LC_INDEX':
                if len(lc_fit['fit']['flux']['lightcurve'].keys())>0:
                    arrays[chara][jgrb] = lc_fit['fit']['flux']['lightcurve']['index']['value']
                    #if not grb_lat['GRBNAME'] in grbs_skipped:
                    nevent_validlc += arrays['NCOUNTS'][jgrb]
                else:
                    arrays[chara][jgrb] = 0
            else:
                arrays[chara][jgrb] = grb_lat[chara]        
        #if not grb_lat['GRBNAME'] in grbs_skipped:
        nevent += arrays['NCOUNTS'][jgrb]
        #    print '  Skipping this GRB...'
        #    continue

    #ncounts_validlc = np.zeros_like(gbm_fluence)
    #ncounts_validintermittent = np.zeros_like(gbm_fluence)
    #ncounts_validaf = np.zeros_like(gbm_fluence)

        #afterglow = pickle_utilities.load('/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/Summary_{name}_afterglow_Eth100MeV_r12deg.pickle'.format(name=grb))
        #if afterglow['lower_energies']['TS']>0:
        #    ncounts_validaf[jgrb] = ncounts[jgrb]
        #    spec_index[jgrb] = afterglow['lower_energies']['Index']['value']
        #    flux_gev[jgrb] = afterglow['lower_energies']['flux']['value']
        #else:
        #    ncounts_validaf[jgrb] = 0
        #    spec_index[jgrb] = 0
        #    flux_gev[jgrb] = 0

        #if (gbm_t90_start[jgrb]+gbm_t90[jgrb])-(gbm_t50_start[jgrb]+1.5*gbm_t50[jgrb])>=10:
        #    ncounts_validintermittent[jgrb] = ncounts[jgrb]
        #else:
        #    ncounts_validintermittent[jgrb] = 0

        #lc = pickle_utilities.load('/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/E0000100-0100000MeV/r12deg/lightcurve/LightCurve_{name}_indexfree_fit.pickle'.format(name=grb))
        #if 'index' in lc['fit']['flux']['lightcurve'] and lc['fit']['flux']['lightcurve']['index']['value']==lc['fit']['flux']['lightcurve']['index']['value']:
        #    lc_index[jgrb] = lc['fit']['flux']['lightcurve']['index']['value']
        #    ncounts_validlc[jgrb] = ncounts[jgrb]
        #else:
        #    lc_index[jgrb] = 0
        #    ncounts_validlc[jgrb] = 0
        #print '  {0} events'.format(ncounts[jgrb])
        jgrb += 1
    
    #nevent = sum(arrays['NCOUNTS'])
    nevent_onethird = nevent/3.
    nevent_validlc_onehalf = nevent_validlc/2.

    #nevent_validlc = sum(ncounts_validlc)
    #nevent_validlc_onehalf = nevent_validlc/2.

    #nevent_validintermittent = sum(ncounts_validintermittent)
    #nevent_validintermittent_onehalf = nevent_validintermittent/2.

    #nevent_validaf = sum(ncounts_validaf)
    #nevent_validaf_onethird = nevent_validaf/3.

    categories = {}
    gbm_fluence_category_sum = {}
    fluence_scaled_sum = {}
    #gbm_fluence_category_sum_err = {}
    charas_quant = ('GBM_T90', 'GBM_FLUENCE', 'GBM_FLUX_1024', 'GBM_FLUX_64', 'LC_INDEX')
    for chara in charas_quant:
        categories[chara] = [[], [], []]
        gbm_fluence_category_sum[chara] = [{'value':0, 'error':0}, {'value':0, 'error':0}, {'value':0, 'error':0}]
        fluence_scaled_sum[chara] = [np.zeros((len(lc_indices), len(lc_phases))),np.zeros((len(lc_indices), len(lc_phases))),np.zeros((len(lc_indices), len(lc_phases)))]

    # categories = {'gbm_t90':[[], [], []], 
    #               'lat_count':[[], [], []], 
    #               'gbm_fluence':[[], [], []], 
    #               'gbm_flux1024':[[], [], []], 
    #               'gbm_flux64':[[], [], []]
    #               #'epeak_band':[[], [], []], 
    #               #'gbm_intermittent':[[], []],
    #               #'gbm_fluence_per_t50':[[], [], []],
    #               #'spec_index':[[], [], []], 
    #               #'flux_gev':[[], [], []], 
    #               #'lc_index':[[], [], []]
    #               }

    dct_df = {}
    for chara in charas:
        dct_df[chara] = arrays[chara]
    df = pd.DataFrame(dct_df)
    # df = pd.DataFrame({ 'name' : names,
    #                     #'redshift' : redshift,
    #                     'gbm_t90' : gbm_t90,
    #                     'gbm_t90_start' : gbm_t90_start,
    #                     'gbm_t50' : gbm_t50,
    #                     'gbm_t50_start' : gbm_t50_start,
    #                     'gbm_flux1024': gbm_flux1024, 
    #                     'gbm_flux64': gbm_flux64,
    #                     #'epeak_band': epeak_band,
    #                     #'gbm_intermittent' : (gbm_t90_start+gbm_t90)-(gbm_t50_start+1.5*gbm_t50),
    #                     #'gbm_fluence_per_t50': gbm_fluence/gbm_t50,
    #                     'lat_count' : ncounts,
    #                     'gbm_fluence' : gbm_fluence
    #                     #'spec_index' : spec_index,
    #                     #'flux_gev' : flux_gev,
    #                     #'lc_index': lc_index
    #                     })

    # GBM fluence, T90
    for col in charas_quant: #('gbm_fluence', 'gbm_t90', 'gbm_flux1024', 'gbm_flux64'): #, 'spec_index', 'flux_gev', 'epeak_band', 'gbm_fluence_per_t50'):
        print '#####', col, '#####'
        df_sorted = df.sort_values(by=col, ascending=True if col in ('LC_INDEX') else False)
        df_sorted.reset_index( drop = True, inplace=True )
        print df_sorted
        mcount = 0
        print '----- Category 1 -----'
        for irow in range(len(df_sorted.index)):
            name = df_sorted.loc[irow,['GRBNAME']][0]
            if name in grbs_skipped:
                print '  Skipping this GRB...'
                continue
            lst_lat_one = ReadLATCatalogueInfo.select_one_by_name(tb_lat, name, tb_gbm)
            gbm_fluence_value = lst_lat_one['GBM']['FLUENCE']
            gbm_fluence_error = lst_lat_one['GBM']['FLUENCE_ERROR']
            scaled_lc[name] = integrated_lc[name]*gbm_fluence_value
            lat_count = df_sorted.loc[irow,['NCOUNTS']][0]
            if col in ('spec_index', 'flux_gev'): # Requires LAT low-E results
                if df_sorted.loc[irow,['spec_index']][0]!=0:
                    mcount += lat_count
                    if mcount < nevent_validaf_onethird:
                        ncategory = 1
                    elif mcount < 2*nevent_validaf_onethird:
                        if len(categories[col][1])<1:
                            print '----- Category 2 -----'
                        ncategory = 2
                    else:
                        if len(categories[col][2])<1:
                            print '----- Category 3 -----'
                        ncategory = 3
                    categories[col][ncategory-1].append(name)
                    gbm_fluence_category_sum[col][ncategory-1]['value'] += gbm_fluence_value
                    gbm_fluence_category_sum[col][ncategory-1]['error'] += pow(gbm_fluence_error,2)
                    fluence_scaled_sum[col][ncategory-1] += scaled_lc[name]
            elif col in ('LC_INDEX'):
                if df_sorted.loc[irow,[col]][0]!=0:
                    mcount += lat_count
                    if mcount < nevent_validlc_onehalf:
                        ncategory = 1
                    else:
                        if len(categories[col][1])<1:
                            print '----- Category 2 -----'
                        ncategory = 2
                else:
                    if len(categories[col][2])<1:
                        print '----- Category 3 -----'
                    ncategory = 3
                categories[col][ncategory-1].append(name)
                gbm_fluence_category_sum[col][ncategory-1]['value'] += gbm_fluence_value
                gbm_fluence_category_sum[col][ncategory-1]['error'] += pow(gbm_fluence_error,2)
                fluence_scaled_sum[col][ncategory-1] += scaled_lc[name]
            else:
                mcount += lat_count
                if mcount < nevent_onethird:
                    ncategory = 1
                elif mcount < 2*nevent_onethird:
                    if len(categories[col][1])<1:
                        print '----- Category 2 -----'
                    ncategory = 2
                else:
                    if len(categories[col][2])<1:
                        print '----- Category 3 -----'
                    ncategory = 3
                categories[col][ncategory-1].append(name)
                gbm_fluence_category_sum[col][ncategory-1]['value'] += gbm_fluence_value
                gbm_fluence_category_sum[col][ncategory-1]['error'] += pow(gbm_fluence_error,2)
                fluence_scaled_sum[col][ncategory-1] += scaled_lc[name]

            print '{name}: {col} (LAT count: {cnt})'.format(name=name, cnt=lat_count, col=df_sorted.loc[irow,[col]][0])
        for icat in (1, 2, 3):
            print 'Category {0}: {1}'.format(icat, categories[col][icat-1])
            print '{0} GRBs'.format(len(categories[col][icat-1]))
            gbm_fluence_category_sum[col][icat-1]['error'] = np.sqrt(gbm_fluence_category_sum[col][icat-1]['error'])
            print 'Summed GBM fluence: {0} +/- {1}'.format(gbm_fluence_category_sum[col][icat-1]['value'], gbm_fluence_category_sum[col][icat-1]['error'])
            print ''


#     col = 'gbm_intermittent'
#     print '#####', 'GBM intermittent time', '#####'
#     df_sorted = df.sort_values(by=col, ascending=False)
#     df_sorted.reset_index( drop = True, inplace=True )
#     print df_sorted
#     mcount = 0
#     print '----- Category 1 -----'
#     for irow in range(len(df_sorted.index)):
#         name = df_sorted.loc[irow,['name']][0]
#         if name in ('130427324', '141022087'):
#             print '  Skipping this GRB...'
#             continue
#         lat_count = df_sorted.loc[irow,['lat_count']][0]
#         time_intermittent = df_sorted.loc[irow,['gbm_intermittent']][0]
#         if time_intermittent>=10:
#             #mcount += lat_count
# #            if mcount < nevent_validintermittent_onehalf:
#             categories[col][0].append(name)
#             #else:
#              #   if len(categories[col][1])<1:
#               #      print '----- Category 2 -----'
#                # categories[col][1].append(name)
#         else:
#             if len(categories[col][1])<1:
#                 print '----- Category 2 -----'
#             categories[col][1].append(name)
#         print '{name}: {col} (LAT count: {cnt})'.format(name=name, cnt=lat_count, col=df_sorted.loc[irow,[col]][0])
#     for icat in (1, 2):
#         print 'Category {0}: {1}'.format(icat, categories[col][icat-1])
#         print '{0} GRBs'.format(len(categories[col][icat-1]))

#     col = 'lc_index'
#     print '#####', 'lightcurve index', '#####'
#     df_sorted = df.sort_values(by=col, ascending=False)
#     df_sorted.reset_index( drop = True, inplace=True )
#     print df_sorted
#     mcount = 0
#     print '----- Category 1 -----'
#     for irow in range(len(df_sorted.index)):
#         name = df_sorted.loc[irow,['name']][0]
#         if name in ('130427324', '141022087'):
#             print '  Skipping this GRB...'
#             continue
#         lat_count = df_sorted.loc[irow,['lat_count']][0]
#         lc_index = df_sorted.loc[irow,['lc_index']][0]
#         if lc_index!=0:
#             mcount += lat_count
#             if mcount < nevent_validlc_onehalf:
#                 categories[col][0].append(name)
#             else:
#                 if len(categories[col][1])<1:
#                     print '----- Category 2 -----'
#                 categories[col][1].append(name)
#         else:
#             if len(categories[col][2])<1:
#                 print '----- Category 3 -----'
#             categories[col][2].append(name)
#         print '{name}: {col} (LAT count: {cnt})'.format(name=name, cnt=lat_count, col=df_sorted.loc[irow,[col]][0])
#     for icat in (1, 2, 3):
#         print 'Category {0}: {1}'.format(icat, categories[col][icat-1])
#         print '{0} GRBs'.format(len(categories[col][icat-1]))
            
    pickle_utilities.dump('{dire}/QuantiledGRBs{suf}.pickle'.format(dire=output, suf=str_output ), 
                          {'categories':categories, 
                           'fluence_summed':gbm_fluence_category_sum, 
                           'fluence_scaled_gbm':{'sum':fluence_scaled_sum, 
                                                 'normalized':integrated_lc, 
                                                 'scaled': scaled_lc,
                                                 'indices':lc_indices, 
                                                 'phases':lc_phases}
                           }
                          )


@click.command()
@click.option('--output', '-o', type=str, default='.')
@click.option('--tolerr', '-t', type=float, default=180.)
@click.option('--suffix', '-s', type=str, default='')
@click.option('--longonly', '-l', is_flag=True)
def main(output, longonly, tolerr, suffix): #name, sed, ra, dec, king, acceptance, livetime, suffix, nside):
    quantile_grbs(output, longonly, tolerr, suffix)


if __name__ == '__main__':
    main()
