#!/usr/bin/env python

import sys
import os
import numpy as np
import pandas as pd
import click
import ReadLTFCatalogueInfo
import ReadGBMCatalogueInfo
import FindGoodstatPeriods
import pickle_utilities

##### PATH of Catalogue #####
GRB_CATALOGUE_LTF = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'


def quantile_grbs(output, suffix):
    tb_ltf = ReadLTFCatalogueInfo.open_table(1, GRB_CATALOGUE_LTF)
    tb_ltf = ReadLTFCatalogueInfo.select_long(tb_ltf)
    tb_gbm = ReadGBMCatalogueInfo.open_table(1)
    gbm_t90 = tb_ltf['T90'].byteswap().newbyteorder()
    gbm_t90_start = tb_ltf['T90_START'].byteswap().newbyteorder()
    gbm_fluence = tb_ltf['FLUENCE'].byteswap().newbyteorder()
    gbm_flux1024 = tb_ltf['FLUX_1024'].byteswap().newbyteorder()
    gbm_flux64 = tb_ltf['FLUX_64'].byteswap().newbyteorder()
    ra = tb_ltf['RA'].byteswap().newbyteorder()
    dec = tb_ltf['DEC'].byteswap().newbyteorder()
    redshift = tb_ltf['REDSHIFT'].byteswap().newbyteorder()
    names = tb_ltf['GRBNAME'].byteswap().newbyteorder()
    gbm_names = tb_ltf['GBM_assoc_key'].byteswap().newbyteorder()
    spec_index = np.zeros_like(gbm_fluence)
    flux_gev = np.zeros_like(gbm_fluence)
    lc_index = np.zeros_like(gbm_fluence)

    # Associate GBM catalogue
    gbm_t50 = np.zeros_like(gbm_t90) 
    gbm_t50_start = np.zeros_like(gbm_t90) 
    epeak_band = np.zeros_like(gbm_t90) 
    for igrb, grb in enumerate(names):
        print grb
        tb_gbm1 = ReadGBMCatalogueInfo.select_one_by_name(tb_gbm, gbm_names[igrb])
        gbm_t50[igrb] = tb_gbm1['T50']
        gbm_t50_start[igrb] = tb_gbm1['T50_START']
        epeak_band[igrb] = tb_gbm1['FLNC_BAND_EPEAK']

    ncounts = np.zeros_like(gbm_fluence)
    ncounts_validlc = np.zeros_like(gbm_fluence)
    ncounts_validintermittent = np.zeros_like(gbm_fluence)
    ncounts_validaf = np.zeros_like(gbm_fluence)

    for igrb, grb in enumerate(names):
        print 'No.{0} {1}'.format(igrb, grb)
        if grb in ('130427324', '141022087'):
            print '  Skipping this GRB...'
            continue

        ncounts[igrb] = FindGoodstatPeriods.get_entries_roi('{name}/E0000100-0100000MeV/r12deg/unified/PowerLaw/IndexFree/GRB{name}_P8_P302_BASE_T00-999-101000_r030_ft1_filtered_gti.fits'.format(name=grb), None, None, rlim=5.0, ra=ra[igrb], dec=dec[igrb])

        afterglow = pickle_utilities.load('/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/Summary_{name}_afterglow_Eth100MeV_r12deg.pickle'.format(name=grb))
        if afterglow['lower_energies']['TS']>0:
            ncounts_validaf[igrb] = ncounts[igrb]
            spec_index[igrb] = afterglow['lower_energies']['Index']['value']
            flux_gev[igrb] = afterglow['lower_energies']['flux']['value']
        else:
            ncounts_validaf[igrb] = 0
            spec_index[igrb] = 0
            flux_gev[igrb] = 0

        if (gbm_t90_start[igrb]+gbm_t90[igrb])-(gbm_t50_start[igrb]+1.5*gbm_t50[igrb])>=10:
            ncounts_validintermittent[igrb] = ncounts[igrb]
        else:
            ncounts_validintermittent[igrb] = 0

        lc = pickle_utilities.load('/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/E0000100-0100000MeV/r12deg/lightcurve/LightCurve_{name}_indexfree_fit.pickle'.format(name=grb))
        if 'index' in lc['fit']['flux']['lightcurve'] and lc['fit']['flux']['lightcurve']['index']['value']==lc['fit']['flux']['lightcurve']['index']['value']:
            lc_index[igrb] = lc['fit']['flux']['lightcurve']['index']['value']
            ncounts_validlc[igrb] = ncounts[igrb]
        else:
            lc_index[igrb] = 0
            ncounts_validlc[igrb] = 0
        print '  {0} events'.format(ncounts[igrb])
    
    nevent = sum(ncounts)
    nevent_onethird = nevent/3.

    nevent_validlc = sum(ncounts_validlc)
    nevent_validlc_onehalf = nevent_validlc/2.

    nevent_validintermittent = sum(ncounts_validintermittent)
    nevent_validintermittent_onehalf = nevent_validintermittent/2.

    nevent_validaf = sum(ncounts_validaf)
    nevent_validaf_onethird = nevent_validaf/3.

    categories = {'gbm_t90':[[], [], []], 
                  'lat_count':[[], [], []], 
                  'gbm_fluence':[[], [], []], 
                  'gbm_flux1024':[[], [], []], 
                  'gbm_flux64':[[], [], []], 
                  'epeak_band':[[], [], []], 
                  'gbm_intermittent':[[], []],
                  'gbm_fluence_per_t50':[[], [], []],
                  'spec_index':[[], [], []], 
                  'flux_gev':[[], [], []], 
                  'lc_index':[[], [], []]}

    df = pd.DataFrame({ 'name' : names,
                        'redshift' : redshift,
                        'gbm_t90' : gbm_t90,
                        'gbm_t90_start' : gbm_t90_start,
                        'gbm_t50' : gbm_t50,
                        'gbm_t50_start' : gbm_t50_start,
                        'gbm_flux1024': gbm_flux1024, 
                        'gbm_flux64': gbm_flux64,
                        'epeak_band': epeak_band,
                        'gbm_intermittent' : (gbm_t90_start+gbm_t90)-(gbm_t50_start+1.5*gbm_t50),
                        'gbm_fluence_per_t50': gbm_fluence/gbm_t50,
                        'lat_count' : ncounts,
                        'gbm_fluence' : gbm_fluence,
                        'spec_index' : spec_index,
                        'flux_gev' : flux_gev,
                        'lc_index': lc_index})

    # GBM fluence, T90
    for col in ('gbm_fluence', 'gbm_t90', 'spec_index', 'flux_gev', 'gbm_flux1024', 'gbm_flux64', 'epeak_band', 'gbm_fluence_per_t50'):
        print '#####', col, '#####'
        df_sorted = df.sort_values(by=col, ascending=False)
        df_sorted.reset_index( drop = True, inplace=True )
        print df_sorted
        mcount = 0
        print '----- Category 1 -----'
        for irow in range(len(df_sorted.index)):
            name = df_sorted.loc[irow,['name']][0]
            if name in ('130427324', '141022087'):
                print '  Skipping this GRB...'
                continue
            lat_count = df_sorted.loc[irow,['lat_count']][0]
            if col in ('spec_index', 'flux_gev'): # Requires LAT low-E results
                if df_sorted.loc[irow,['spec_index']][0]!=0:
                    mcount += lat_count
                    if mcount < nevent_validaf_onethird:
                        categories[col][0].append(name)
                    elif mcount < 2*nevent_validaf_onethird:
                        if len(categories[col][1])<1:
                            print '----- Category 2 -----'
                        categories[col][1].append(name)
                    else:
                        if len(categories[col][2])<1:
                            print '----- Category 3 -----'
                        categories[col][2].append(name)
            else:
                mcount += lat_count
                if mcount < nevent_onethird:
                    categories[col][0].append(name)
                elif mcount < 2*nevent_onethird:
                    if len(categories[col][1])<1:
                        print '----- Category 2 -----'
                    categories[col][1].append(name)
                else:
                    if len(categories[col][2])<1:
                        print '----- Category 3 -----'
                    categories[col][2].append(name)
            print '{name}: {col} (LAT count: {cnt})'.format(name=name, cnt=lat_count, col=df_sorted.loc[irow,[col]][0])
        for icat in (1, 2, 3):
            print 'Category {0}: {1}'.format(icat, categories[col][icat-1])
            print '{0} GRBs'.format(len(categories[col][icat-1]))


    col = 'gbm_intermittent'
    print '#####', 'GBM intermittent time', '#####'
    df_sorted = df.sort_values(by=col, ascending=False)
    df_sorted.reset_index( drop = True, inplace=True )
    print df_sorted
    mcount = 0
    print '----- Category 1 -----'
    for irow in range(len(df_sorted.index)):
        name = df_sorted.loc[irow,['name']][0]
        if name in ('130427324', '141022087'):
            print '  Skipping this GRB...'
            continue
        lat_count = df_sorted.loc[irow,['lat_count']][0]
        time_intermittent = df_sorted.loc[irow,['gbm_intermittent']][0]
        if time_intermittent>=10:
            #mcount += lat_count
#            if mcount < nevent_validintermittent_onehalf:
            categories[col][0].append(name)
            #else:
             #   if len(categories[col][1])<1:
              #      print '----- Category 2 -----'
               # categories[col][1].append(name)
        else:
            if len(categories[col][1])<1:
                print '----- Category 2 -----'
            categories[col][1].append(name)
        print '{name}: {col} (LAT count: {cnt})'.format(name=name, cnt=lat_count, col=df_sorted.loc[irow,[col]][0])
    for icat in (1, 2):
        print 'Category {0}: {1}'.format(icat, categories[col][icat-1])
        print '{0} GRBs'.format(len(categories[col][icat-1]))

    col = 'lc_index'
    print '#####', 'lightcurve index', '#####'
    df_sorted = df.sort_values(by=col, ascending=False)
    df_sorted.reset_index( drop = True, inplace=True )
    print df_sorted
    mcount = 0
    print '----- Category 1 -----'
    for irow in range(len(df_sorted.index)):
        name = df_sorted.loc[irow,['name']][0]
        if name in ('130427324', '141022087'):
            print '  Skipping this GRB...'
            continue
        lat_count = df_sorted.loc[irow,['lat_count']][0]
        lc_index = df_sorted.loc[irow,['lc_index']][0]
        if lc_index!=0:
            mcount += lat_count
            if mcount < nevent_validlc_onehalf:
                categories[col][0].append(name)
            else:
                if len(categories[col][1])<1:
                    print '----- Category 2 -----'
                categories[col][1].append(name)
        else:
            if len(categories[col][2])<1:
                print '----- Category 3 -----'
            categories[col][2].append(name)
        print '{name}: {col} (LAT count: {cnt})'.format(name=name, cnt=lat_count, col=df_sorted.loc[irow,[col]][0])
    for icat in (1, 2, 3):
        print 'Category {0}: {1}'.format(icat, categories[col][icat-1])
        print '{0} GRBs'.format(len(categories[col][icat-1]))
            
    pickle_utilities.dump('{dire}/QuantiledGRBs{suf}.pickle'.format(dire=output, suf='_'+suffix if suffix!='' else suffix), categories)


@click.command()
#@click.argument('name', type=str)
#@click.argument('dst', nargs=-1)
@click.option('--output', '-o', type=str, default='.')
@click.option('--suffix', '-s', type=str, default='')
#@click.option('--values', multiple=True)
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
#@click.option('--shout', is_flag=True)
def main(output, suffix): #name, sed, ra, dec, king, acceptance, livetime, suffix, nside):
    quantile_grbs(output, suffix)


if __name__ == '__main__':
    main()
