#!/usr/bin/env python

import os
import sys
from astropy.io import fits
import xml.etree.ElementTree as ET
import numpy as np
import click
import ReadGBMCatalogueInfo

# Catalogue Path
PATH_CATALOGUE = "/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LATBurstCatalogue.xml"
PATH_CATALOGUE_GBM = "/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/GBMcatalogue20171005.fits"
# GBM name association
DCT_GBM_NAME_ASSOC = {"080818945" : "GRB080818945",
                      "080825593" : "GRB080825593",
                      "080916009" : "GRB080916009",
                      "081006604" : "GRB081006604",
                      "081009140" : "GRB081009140",
                      "081024891" : "GRB081024891",
                      "081102365" : "GRB081102365",
                      "081122520" : "GRB081122520",
                      "081203581" : None        ,
                      "081224887" : "GRB081224887",
                      "090217206" : "GRB090217206",
                      "090227310" : "GRB090227310",
                      "090227772" : "GRB090227772",
                      "090228204" : "GRB090228204",
                      "090323002" : "GRB090323002",
                      "090328401" : "GRB090328401",
                      "090424592" : "GRB090424592",
                      "090427977" : None        ,
                      "090510016" : "GRB090510016",
                      "090531775" : "GRB090531775",
                      "090626189" : "GRB090626189",
                      "090720710" : "GRB090720710",
                      "090902462" : "GRB090902462",
                      "090926181" : "GRB090926181",
                      "091003191" : "GRB091003191",
                      "091003192" : "GRB091003191",
                      "091031500" : "GRB091031500",
                      "091120191" : "GRB091120191",
                      "091127976" : "GRB091127976",
                      "091208410" : "GRB091208410",
                      "100116897" : "GRB100116897",
                      "100225115" : "GRB100225115",
                      "100325275" : "GRB100325275",
                      "100414097" : "GRB100414097",
                      "100423244" : "GRB100423244",
                      "100511035" : "GRB100511035",
                      "100620119" : "GRB100620119",
                      "100724029" : "GRB100724029",
                      "100728095" : "GRB100728095",
                      "100826957" : "GRB100826957",
                      "101014175" : "GRB101014175",
                      "101107011" : "GRB101107011",
                      "101123952" : "GRB101123952",
                      "101227406" : "GRB101227406",
                      "110120666" : "GRB110120666",
                      "110123804" : "GRB110123804",
                      "110227229" : "GRB110227229",
                      "110328520" : "GRB110328520",
                      "110428388" : "GRB110428388",
                      "110518860" : None,
                      "110529034" : "GRB110529034",
                      "110625881" : "GRB110625881",
                      "110709642" : "GRB110709642",
                      "110721200" : "GRB110721200",
                      "110728056" : "GRB110728056",
                      "110731465" : "GRB110731465",
                      "110903111" : "GRB110903111",
                      "110921912" : "GRB110921912",
                      "120107384" : "GRB120107384",
                      "120226871" : "GRB120226871",
                      "120316008" : "GRB120316008",
                      "120328268" : "GRB120328268",
                      "120420858" : "GRB120420858",
                      "120526303" : "GRB120526303",
                      "120624933" : "GRB120624933",
                      "120709883" : "GRB120709883",
                      "120711115" : "GRB120711115",
                      "120729456" : "GRB120729456",
                      "120830297" : "GRB120830297",
                      "120911268" : None        ,
                      "120915000" : "GRB120915000",
                      "120916173" : "GRB120916173",
                      "121011469" : "GRB121011469",
                      "121029350" : "GRB121029350",
                      "121123442" : "GRB121123442",
                      "121216419" : "GRB121216419",
                      "121225417" : "GRB121225417",
                      "130206817" : "GRB130206817",
                      "130228111" : "GRB130228111",
                      "130305486" : "GRB130305486",
                      "130310840" : "GRB130310840",
                      "130325203" : "GRB130325203",
                      "130327350" : "GRB130327350",
                      "130427324" : "GRB130427324",
                      "130502327" : "GRB130502327",
                      "130504314" : "GRB130504314",
                      "130504979" : "GRB130504978",
                      "130504978" : "GRB130504978",
                      "130518580" : "GRB130518580",
                      "130606497" : "GRB130606497",
                      "130702004" : "GRB130702004",
                      "130804023" : "GRB130804023",
                      "130821674" : "GRB130821674",
                      "130828306" : "GRB130828306",
                      "130907904" : None,
                      "131014215" : "GRB131014215",
                      "131018673" : "GRB131018673",
                      "131029973" : "GRB131029973",
                      "131108862" : "GRB131108862",
                      "131209547" : "GRB131209547",
                      "131216081" : "GRB131216081",
                      "131231198" : "GRB131231198",
                      "140102887" : "GRB140102887",
                      "140104731" : "GRB140104731",
                      "140110263" : "GRB140110263",
                      "140124527" : "GRB140124527",
                      "140204547" : "GRB140204547",
                      "140206275" : "GRB140206275",
                      "140219824" : "GRB140219824",
                      "140323433" : "GRB140323433",
                      "140329295" : "GRB140329295",
                      "140402007" : "GRB140402007",
                      "140416060" : "GRB140416060",
                      "140523129" : "GRB140523129",
                      "140528837" : "GRB140528837",
                      "140619475" : "GRB140619475",
                      "140723067" : "GRB140723067",
                      "140724533" : "GRB140724533",
                      "140729026" : "GRB140729026",
                      "140810782" : "GRB140810782",
                      "140928437" : "GRB140928437",
                      "141012773" : "GRB141012773",
                      "141022087" : "GRB141022087",
                      "141028455" : "GRB141028455",
                      "141102536" : "GRB141102536",
                      "141113346" : "GRB141113346",
                      "141207800" : "GRB141207800",
                      "141221897" : "GRB141221897",
                      "141222298" : "GRB141222298",
                      "150118409" : "GRB150118409",
                      "150127398" : "GRB150127398",
                      "150202999" : "GRB150202999",
                      "150210935" : "GRB150210935",
                      "150314205" : "GRB150314205",
                      "150403913" : "GRB150403913",
                      "150416773" : "GRB150416773",
                      "150422703" : "GRB150422703",
                      "150510139" : "GRB150510139",
                      "150513855" : "GRB150513856",
                      "150514774" : "GRB150514774",
                      "150523396" : "GRB150523396",
                      "150627183" : "GRB150627183",
                      "150702998" : "GRB150702998",
                      "150724782" : "GRB150724782",
                      "150902733" : "GRB150902733",
                      "150913161" : "GRB150913161",
                      "151006413" : "GRB151006413",
                      "160310016" : "GRB160310016",
                      "160314929" : "GRB160314929",
                      "160325291" : "GRB160325291",
                      "160422499" : "GRB160422499",
                      "160101215" : "GRB160101215",
                      #"171010792" : "GRB171010792",
                      "170906030" : "GRB170906030",
                      "170810918" : "GRB170810918",
                      "170808936" : "GRB170808936",
                      "170522657" : "GRB170522657",
                      "170510217" : "GRB170510217",
                      "170409112" : "GRB170409112",
                      "170405777" : "GRB170405777",
                      "170329387" : "GRB170329387",
                      "170306588" : "GRB170306588",
                      "170228794" : "GRB170228794",
                      "170214649" : "GRB170214649",
                      "170115743" : "GRB170115743",
                      "161202970" : None, #"GRB161202970",
                      "161109263" : "GRB161109263",
                      "160910722" : "GRB160910722",
                      "160905471" : "GRB160905471",
                      "160829334" : "GRB160829334",
                      "160821859" : "GRB160821857", #"GRB160821859",
                      "160816730" : "GRB160816730",
                      "160709826" : "GRB160709826",
                      "160625945" : "GRB160625945",
                      "160623209" : "GRB160623209",
                      "160521385" : "GRB160521385",
                      "160509374" : "GRB160509374",
                      "160503567" : "GRB160503567"
                      }

def open_table(path=PATH_CATALOGUE):
    """Open your XML file and return its root.
"""
    if path in (None, ""):
        path=PATH_CATALOGUE
    f = ET.parse(path)
    rtXml = f.getroot()
    return rtXml


def read_one_row(grb, tb_gbm=None):
    #print 'Reading {0}...'.format(grb)
    dct_info = {}
    dct_info['GRBNAME'] = grb.findtext("./GRBNAME")
    #print dct_info['GRBNAME']
    dct_info['GCNNAME'] = grb.findtext("./GCNNAME")
    dct_info['LAT_TRIGGER_TIME'] = float(grb.findtext("./MET"))
    dct_info['TRIGGER_TIME'] = dct_info['LAT_TRIGGER_TIME']
    if not grb.findtext("./TS") in ("--", "", "NA"):
        dct_info['LAT_TS'] = float(grb.findtext("./TS"))
    else:
        dct_info['LAT_TS'] = 0

    if not grb.findtext("./ERROR") in ("--", "","NA") and not grb.findtext("./LATERROR") in ("--", "", "NA"):
        if float(grb.findtext("./ERROR"))<=float(grb.findtext("./LATERROR")):
            dct_info['ERROR'] = float(grb.findtext("./ERROR"))
            dct_info['RA'] = float(grb.findtext("./RA"))
            dct_info['DEC'] = float(grb.findtext("./DEC"))
        else:
            dct_info['ERROR'] = float(grb.findtext("./LATERROR")) 
            dct_info['RA'] = float(grb.findtext("./LATRA"))
            dct_info['DEC'] = float(grb.findtext("./LATDEC"))
    elif not grb.findtext("./ERROR") in ("--", "", "NA"):
        dct_info['ERROR'] = float(sys.maxint)
        dct_info['RA'] = float(grb.findtext("./RA"))
        dct_info['DEC'] = float(grb.findtext("./DEC"))
    elif not grb.findtext("./LATERROR") in ("--", "", "NA"):
        dct_info['ERROR'] = float(grb.findtext("./LATERROR")) 
        dct_info['RA'] = float(grb.findtext("./LATRA"))
        dct_info['DEC'] = float(grb.findtext("./LATDEC"))
    else:
        dct_info['ERROR'] = sys.maxint
        dct_info['RA'] = float(grb.findtext("./RA"))
        dct_info['DEC'] = float(grb.findtext("./DEC"))
        
    if tb_gbm is not None:
        if not dct_info['GRBNAME'] in DCT_GBM_NAME_ASSOC:
            print 'No GBM association in the dictionary!'
        elif DCT_GBM_NAME_ASSOC[dct_info['GRBNAME']] is not None:
            dct_info['GBM'] = ReadGBMCatalogueInfo.select_one_by_name(tb_gbm, DCT_GBM_NAME_ASSOC[dct_info['GRBNAME']])
            dct_info['TRIGGER_TIME'] = dct_info['GBM']['TRIGGER_TIME']
        else:
            print dct_info['GRBNAME'], ': NO GBM observation!!'
    # Redshift is not implemented yet.
    dct_info['REDSHIFT'] = 0
    return dct_info


def read_all(root_xml, tb_gbm=None):
    lst_info = []
    dct_info = {}
    for grb in root_xml:
        dct_info[grb.findtext("./GRBNAME")] = read_one_row(grb, tb_gbm)
    for k, v in sorted(dct_info.items()):
        lst_info.append(v)
    return lst_info


def select_one_by_name(root_xml, grbname, tb_gbm=None):
    """Return a masked table consists of one paticular GRB.
"""
    dct_info = {}
    for grb in root_xml:
        if grb.findtext("./GRBNAME")==grbname: 
            one_grb = read_one_row(grb, tb_gbm)
            return one_grb
    print 'No data of {0} in the LAT list!!'.format(grbname)
    return 1


def select_by_name(root_xml, name_min, name_max='999999999', tb_gbm=None):
    """Return a masked table consists of several GRBs.
"""
    lst_info = []
    dct_info = {}
    for grb in root_xml:
        if float(grb.findtext("./GRBNAME")) >= float(name_min) and float(grb.findtext("./GRBNAME")) <= float(name_max):
            dct_info[grb.findtext("./GRBNAME")] = read_one_row(grb, tb_gbm)
    for k, v in sorted(dct_info.items()):
        lst_info.append(v)
    print 'Selected', len(lst_info), 'GRBs by name.' 
    return lst_info


def remove_by_name(lst_grbs, excludes):
    if excludes is None:
        return lst_grbs
    lst_new = []
    for grb in lst_grbs:
        if not grb['GRBNAME'] in excludes:
            lst_new.append(grb)
    print 'Updated list after removing {0}:'.format(excludes)
    print lst_new
    return lst_new


def select_gbm_exist(lst_grbs):
    """Return a masked table consists of GRBs with GBM catalogue.
"""
    lst_info = []
    dct_info = {}
    for grb in lst_grbs:
        if 'GBM' in grb:
            dct_info[grb["GRBNAME"]] = grb
    for k, v in sorted(dct_info.items()):
        lst_info.append(v)
    print 'Selected', len(lst_info), 'GRBs by existance of GBM data.' 
    return lst_info


def select_small_error(lst_grbs, rad_tol=0.3):
    """Return a masked table consists of long-GRBs based on GBM catalogue.
"""
    lst_info = []
    dct_info = {}
    for grb in lst_grbs:
        if grb['ERROR']<rad_tol:
            dct_info[grb["GRBNAME"]] = grb
    for k, v in sorted(dct_info.items()):
        lst_info.append(v)
    print 'Selected', len(lst_info), 'GRBs with localization error smaller than', rad_tol, 'deg.' 
    return lst_info


def select_long(lst_grbs):
    """Return a masked table consists of long-GRBs based on GBM catalogue.
"""
    lst_info = []
    dct_info = {}
    for grb in lst_grbs:
        if 'GBM' in grb and grb['GBM']['T90']>=2.0:
            #lst_info.append(o)
            dct_info[grb["GRBNAME"]] = grb
    for k, v in sorted(dct_info.items()):
        lst_info.append(v)
    print 'Selected', len(lst_info), 'long-GRBs.' 
    return lst_info


def select_short(lst_grbs):
    """Return a masked table consists of short-GRBs based on GBM catalogue.
"""
    lst_info = []
    dct_info = {}
    for grb in lst_grbs:
        if 'GBM' in grb and grb['GBM']['T90']<2.0:
            dct_info[grb["GRBNAME"]] = grb
    for k, v in sorted(dct_info.items()):
        lst_info.append(v)
    print 'Selected', len(lst_info), 'short-GRBs.' 
    return lst_info


@click.command()
@click.argument('grb', type=str)
@click.option('--xml', '-x',default=PATH_CATALOGUE, help="Path of XML catalogue.")
@click.option('--gbm', '-t',default=PATH_CATALOGUE_GBM, help="Path of FITS table.")
@click.option('--items', '-i', multiple=True, default=None)
@click.option('--gbmitems', multiple=True, default=None)
def main(grb, xml, gbm, items, gbmitems):
    tb = open_table(xml)
    tb_gbm = ReadGBMCatalogueInfo.open_table(1, gbm)
    tb_masded = select_one_by_name(tb, grb, tb_gbm)
    if items is None:
        print tb_masded
    elif len(items)>0:
        for item in items:
            print item, tb_masded[item]
    if len(gbmitems)>0:
        print 'GBM'
        for item in gbmitems:
            print item, tb_masded['GBM'][item]


if __name__ == '__main__':
    main()
