#!/usr/bin/env python

import numpy as np
import pandas as pd
import logging
from logging import getLogger, StreamHandler
import click


##### Logging #####
logger = getLogger(__name__)
handler = StreamHandler()
handler.setLevel('WARNING')
logger.setLevel('WARNING')
logger.addHandler(handler)


DCT_TABLES = {'burst_info':'table_1_burst_info.csv',
              'durations':'table_2_durations.csv',
              'spectral_fits':'table_3_spectral_fits.csv',
              'energetics':'table_4_energetics.csv',
              'collimation':'table_5_collimation.csv'}


def open_table(table=None, path_catalogue='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/KW_GRBs_z_trig_tables'):
    """Open the ASCII catalogue file you assigned and return the data frame.
"""
    if table is None:
        return [pd.read_csv('/'.join([path_catalogue, t])) for t in DCT_TABLES.values()]
    else:
        path_file = '/'.join([path_catalogue, DCT_TABLES[table]])
        return pd.read_csv(path_file)


def select_one_by_name(tb, name, specmodel='Band', latname=False):
    """Return a masked table consists of one paticular GRB. The name must be string.
"""
    if latname==True:
        grbname = 'GRB'+DCT_LAT_NAME_ASSOC[name]
    else:
        grbname = 'GRB'+name
    if 'Model' in tb.columns:
        tb_masked = tb[(tb.Burst_name == grbname)&(tb.Model == specmodel)&(tb.Spec_type == 'i')]
        #print 'Spec table:',tb_masked
    else:
        tb_masked = tb[tb.Burst_name == grbname]
    logger.debug(tb_masked)
    if len(tb_masked)==1:
        return tb_masked
    elif len(tb_masked)<1:
        logger.info('No data of {0} in the Konus-Wind redshift-known GRB catalogue.'.format(grbname))
        return 1
    elif len(tb_masked)>1:
        logger.info('{0} ambiguous candidates of {0} in the Konus-Wind redshift-known GRB catalogue.'.format(grbname))
        return 1


def select_by_name(tb, name_min='0', name_max='200000Z'):
    """Return a masked table consists of GRBs from name_min to name_max. The names should be string.
"""
    grbname_min = 'GRB'+name_min
    grbname_max = 'GRB'+name_max
    return tb[tb.Burst_name>=grbname_min and tb.Burst_name<=grbname_max] #tb[(tb['NAME'] >= name_min) * (tb['NAME'] < name_max)]


@click.command()
@click.argument('grb', type=str)
@click.option('--table', '-t',default="energetics", help="Name of the table.")
@click.option('--items', '-i', multiple=True, default=None)
def main(grb, table, items):
    ##### Logger #####
    loglevel='INFO'
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    tb = open_table(table)
    logger.debug(tb.columns)
    tb_masded = select_one_by_name(tb, grb)
    if items is None:
        logger.info(tb_masded)
    elif len(items)>1:
        logger.info(tb_masded[list(items)])
    elif len(items)==1:
        logger.info(tb_masded[items[0]])
    elif len(items):
        logger.info(tb_masded[items])


if __name__ == '__main__':
    main()


# LAT name association
DCT_LAT_NAME_ASSOC = {'171212434':'171212B',
                      '171210493':'171210A',
                      '171124235':'171124A',
                      '171120556':'171120A',
                      '171102107':'171102A',
                      '171022885':'171022A',
                      '171010792':'171010A',
                      '170906030':'170906A',
                      '170810918':'170810A',
                      '170808936':'170808B',
                      '170522657':'170522A',
                      '170510217':'170510A',
                      '170409112':'170409A',
                      '170405777':'170405A',
                      '170329387':'170329A',
                      '170306588':'170306B',
                      '170228794':'170228A',
                      '170214649':'170214A',
                      '170115743':'170115B',
                      '161202970':'161202A',
                      '161109263':'161109A',
                      '160910722':'160910A',
                      '160905471':'160905A',
                      '160829334':'160829A',
                      '160821859':'160821A',
                      '160816730':'160816A',
                      '160709826':'160709A',
                      '160625945':'160625B',
                      '160623209':'160623A',
                      '160521385':'160521B',
                      '160509374':'160509A',
                      '160503567':'160503A',
                      '160422499':'160422A',
                      '160325291':'160325A',
                      '160314929':'160314B',
                      '160310016':'160310A',
                      '160101215':'160101B',
                      '151006413':'151006A',
                      '150902733':'150902A',
                      '150724782':'150724B',
                      '150702998':'150702A',
                      '150627183':'150627A',
                      '150523396':'150523A',
                      '150514774':'150514A',
                      '150513855':'150513A',
                      '150510139':'150510A',
                      '150416773':'150416A',
                      '150403913':'150403A',
                      '150314205':'150314A',
                      '150210935':'150210A',
                      '150202999':'150202B',
                      '150127398':'150127A',
                      '150118409':'150118B',
                      '141222298':'141222A',
                      '141207800':'141207A',
                      '141102536':'141102A',
                      '141028455':'141028A',
                      '140928437':'140928A',
                      '140810782':'140810A',
                      '140729026':'140729A',
                      '140723067':'140723A',
                      '140619475':'140619B',
                      '140523129':'140523A',
                      '140402007':'140402A',
                      '140329295':'140329A',
                      '140323433':'140323A',
                      '140219824':'140219A',
                      '140206275':'140206B',
                      '140110263':'140110A',
                      '140104731':'140104B',
                      '140102887':'140102A',
                      '131231198':'131231A',
                      '131216081':'131216A',
                      '131209547':'131209A',
                      '131108862':'131108A',
                      '131029973':'131029A',
                      '131018673':'131018B',
                      '131014215':'131014A',
                      '130907904':'130907A',
                      '130828306':'130828A',
                      '130821674':'130821A',
                      '130702004':'130702A',
                      '130606497':'130606B',
                      '130518580':'130518A',
                      '130504979':'130504C',
                      '130502327':'130502B',
                      '130427324':'130427A',
                      '130327350':'130327B',
                      '130325203':'130325A',
                      '130310840':'130310A',
                      '130305486':'130305A',
                      '130228111':'130228A',
                      '130206817':'130206A',
                      '121225417':'121225B',
                      '121011469':'121011A',
                      '120916173':'120916A',
                      '120911268':'120911B',
                      '120830297':'120830A',
                      '120729456':'120729A',
                      '120711115':'120711A',
                      '120709883':'120709A',
                      '120624933':'120624B',
                      '120328268':'120328B',
                      '120316008':'120316A',
                      '120226871':'120226A',
                      '120107384':'120107A',
                      '110731465':'110731A',
                      '110721200':'110721A',
                      '110709642':'110709A',
                      '110625881':'110625A',
                      '110529034':'110529A',
                      '110428388':'110428A',
                      '110328520':'110328B',
                      '110120666':'110120A',
                      '101123952':'101123A',
                      '101014175':'101014A',
                      '100826957':'100826A',
                      '100728095':'100728A',
                      '100724029':'100724B',
                      '100620119':'100620A',
                      '100414097':'100414A',
                      '100325275':'100325A',
                      '100225115':'100225A',
                      '100116897':'100116A',
                      '091208410':'091208B',
                      '091031500':'091031',
                      '091003191':'091003',
                      '090926181':'090926A',
                      '090902462':'090902B',
                      '090720710':'090720B',
                      '090626189':'090626',
                      '090531775':'090531B',
                      '090510016':'090510',
                      '090328401':'090328',
                      '090323002':'090323',
                      '090227772':'090227B',
                      '090217206':'090217',
                      '081024891':'081024B',
                      '081006604':'081006',
                      '080916009':'080916C',
                      '080825593':'080825C'}

DCT_FLUENCE_GCN = {'110625A': {'value': 6.1E-5, 'err_hi':np.nan, 'err_lo':np.nan}, #(6.1 =B1 0.6)x10-5 erg/cm2
                   '110709A': {'value': 3.7E-5, 'err_hi':np.nan, 'err_lo':np.nan}, #(3.7 =B1 0.3)x10-5 erg/cm2
                   '130305A': {'value': 8.5E-5 , 'err_hi':np.nan, 'err_lo':np.nan},#(8.5 =B1 0.6)x10-5 erg/cm2 +  (3.3 =B1 1.1)x10-5 erg/cm2 http://www.mpe.mpg.de/~jcg/grb130305A.html
                   '130504C': {'value': 2.0E-4, 'err_hi':np.nan, 'err_lo':np.nan}, #(2.0 =B1 0.1)x10-4 erg/cm2
                   '130606B': {'value': 4.3E-4 , 'err_hi':np.nan, 'err_lo':np.nan}, #(4.3 =B1 0.1)x10-4 erg/cm2
                   '130702A': {'value':6.70E-6, 'err_hi':+0.82E-6 , 'err_lo':0.80E-6}, #6.70(-0.80,+0.82)10^-6 erg/cm2
                   '131014A': {'value':2.05E-4 , 'err_hi':np.nan, 'err_lo':np.nan}, #(2.05 =B1 0.03)x10-4 erg/cm2
                   '140102A': {'value': 2.0E-5, 'err_hi':np.nan, 'err_lo':np.nan}, #(2.0 =B1 0.2)x10-5 erg/cm2
                   '140323A': {'value': 3.1E-5, 'err_hi':0.7E-5, 'err_lo':0.4E-5}, #3.1(-0.4,+0.7)x10^-5 erg/cm2
                   '140928A': {'value': 8.6E-5, 'err_hi':1.3E-5, 'err_lo':1.4E-5}, #8.6(-1.4,+1.3)x10^-5 erg/cm2
                   '150523A': {'value': 4.44E-5, 'err_hi':+0.71E-5, 'err_lo':0.69E-5}, #4.44(-0.69,+0.71)x10^-5 erg/cm2
                   '150724B': {'value': 3.63E-5, 'err_hi':+0.69E-5, 'err_lo':0.29E-5}, #3.63(-0.29,+0.69)x10^-5 erg/cm2
                   '150902A': {'value': 1.21E-4, 'err_hi':0.07E-4, 'err_lo':0.07E-4}, #1.21(-0.07,+0.07)x10^-4 erg/cm2
                   '151006A': {'value': 5.30E-5, 'err_hi':0.92E-5, 'err_lo':0.87E-5}, #5.30(-0.87,+0.92)x10^-5 erg/cm2
                   '160325A': {'value': 1.73E-5, 'err_hi':0.14E-5, 'err_lo':0.12E-5}, #1.73(-0.12,+0.14)x10^-5 erg/cm2
                   '160521B': {'value': 1.32E-5, 'err_hi':0.17E-5, 'err_lo':0.15E-5},#1.32(-0.15,+0.17)x10^-5 erg/cm2
                   '160816A': {'value': 3.16E-5, 'err_hi':0.11E-5, 'err_lo':0.10E-5},#3.16(-0.10,+0.11)x10^-5 erg/cm2
                   '160905A': {'value': 1.50E-4, 'err_hi':0.11E-4, 'err_lo':0.11E-4}, #1.50(-0.11,+0.11)x10^-4 erg/cm2
                   '160910A': {'value': 1.36E-4, 'err_hi':0.10E-4, 'err_lo':0.09E-4}, # 1.36(-0.09,+0.10)x10^-4 erg/cm2
                   '170214A': {'value': 2.41E-4, 'err_hi':np.nan, 'err_lo':np.nan}, #(2.41 =B1 0.17)x10^-4 erg/cm2
                   '170906A': {'value': 1.19E-4, 'err_hi':0.22E-4, 'err_lo':0.21E-4}, #1.19(-0.21,+0.22)x10^-4 erg/cm2}
                   }
