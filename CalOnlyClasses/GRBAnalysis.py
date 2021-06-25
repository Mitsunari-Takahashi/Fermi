#!/usr/bin/env python

import sys
import os
import os.path as path
from astropy.table import Table
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import matplotlib as mpl
import matplotlib.pyplot as plt
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
ROOT.gROOT.SetBatch()
from pColor import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import EventData
import Catalogue


##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


FORMAT_DICT = {'fits': 'fits',
               'tsv': 'ascii.tab',
               'csv': 'ascii.csv'}
    

def analyze_grb(name, event_paths, ra, dec, t0, time_intervals=[{'>=':-100, '<':10000}], dist_max='PSF95', event_class=4, zenith_max=90.):
    center_coordinates = {'ra':ra, 'dec':dec}
    event_data = EventData.EventDataAnnuli(center_coordinates=center_coordinates, rinterval=[{'<':dist_max}], data_path=None, intervals={'EVENT_CLASS':[{'>=':event_class}], 'TIME':time_intervals, 'ENERGY':[{'>=':10**4.35, '<':10**5.75}],'ZENITH_ANGLE':[{'<':zenith_max}]}, offset_time=t0)
    for event_path in event_paths:
        event_data.add_events(event_path)
    event_data.add_angular_distance_friend(origin_name='center', origin_coordinates=center_coordinates)
    evttree = event_data.cut_events()
    for evt in evttree:
        logger.info('--------------------')
        str_evt = '''Class: {cla}
Arrival time: {atime:.2f} s after trigger time
Energy: {egev:.1f} GeV
Angular distance: {angdist:.1f} deg
PSF68%: {psf68:.1f} deg
PSF95%: {psf95:.1f} deg
RA: {ra:.1f} deg
DEC: {dec:.1f} deg
'''.format(cla=evt.s, atime=evt.t-t0, egev=evt.e/1e3, aangdist=evt.ANG_DIST, psf68=evt.PSF68, psf95=evt.PSF95, ra=evt.ra, dec=evt.dec)
        logger.warning(str_evt)
    
    
@click.command()
@click.argument('bursts', type=str)
@click.option('--eventfile', type=str, default='/nfs/farm/g/glast/u/mtakahas/data/EVENTS/AllSky/trAllSkyMap_CalOnlyDownlike_2008-2020.root')
@click.option('--suffix', type=str, default='')
@click.option('--tmin', type=float, default=-100.)
@click.option('--tmax', type=float, default=10000.)
@click.option('--emin', type=float, default=10**4.35)
@click.option('--emax', type=float, default=10**5.75)
@click.option('--catalogue', type=click.Choice(['GBM', 'LAT', 'Swift']), default='GBM')
@click.option('--bkg', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(bursts, eventfile, catalogue, suffix, tmin, tmax, emin, emax, bkg, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    if catalogue=='Swift':
        catalogue = Catalogue.SwiftCatalogue('Swift catalogue', '/nfs/farm/g/glast/u/mtakahas/data/catalogue/grb_table_xrt.csv')
    elif catalogue=='LAT':
        catalogue = Catalogue.XMLCatalogue('LAT catalogue', '/nfs/farm/g/glast/u/mtakahas/data/catalogue/PublicTableGRBs.xml')
    elif catalogue=='GBM':
        catalogue = Catalogue.GBMCatalogue('GBM catalogue', '/nfs/farm/g/glast/u/mtakahas/data/catalogue/GBM_catalog_20210625.fits')
        
    grb_list = []
    if path.isfile(bursts):
        with open(bursts, 'r') as grbfile:
            lines = grbfile.readlines()
        for line in lines:
            grb_list.append(line.strip('\n'))
            logger.debug(grb_list[-1])
    else:
        grb_list.append(bursts)

    for grb in grb_list:
        logger.info('\n===== {grb} ====='.format(grb))
        grb_info = catalogue.find_grb(grb)
        analyze_grb(name=grb, event_paths=eventfile, ra=grb_info['RA'], dec=grb_info['DEC'], t0=grb_info['MET'], time_intervals=[{'>=':-100, '<':10000}], dist_max='PSF95', event_class=4, zenith_max=90.)

        
if __name__ == '__main__':
    main()
