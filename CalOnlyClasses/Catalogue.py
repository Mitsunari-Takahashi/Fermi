#!/usr/bin/env python

import sys
import os
import os.path as path
from astropy.table import Table
import xml.etree.ElementTree as ET
from astropy.time import Time
import numpy as np
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
from pMETandMJD import ConvertMetToMjd, ConvertMjdToMet


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


class Catalogue:
    def __init__(self, name, file_path):
        self.name = name
        self.file_path = file_path
        self.file_ext = path.splitext(path.basename(file_path))[1]


class AstropyTableCatalogue(Catalogue):
    def __init__(self, name, file_path):
        Catalogue.__init__(self, name=name, file_path=file_path)
        if not self.file_ext.lower() in FORMAT_DICT.keys():
            logger.error('File extenstion must be one of {0}!!'.format(FORMAT_DICT.keys()))
        self.catalogue_table = Table.read(self.file_path, format=FORMAT_DICT[self.file_ext])


class SwiftCatalogue(AstropyTableCatalogue):
    def __init__(self, name, file_path):
        AstropyTableCatalogue.__init__(self, name=name, file_path=file_path)

        
    def find_grb(self, grbname):
        grb_info = {}
        grb = self.catalogue_table[self.catalogue_table['NAME']==grbname][0]
        grb_info['NAME'] = grb['GRB']
        astrotime = Time('20{y}-{m}-{d}T{t}'.format(y=grb_info['NAME'][:2], m=grb_info['NAME'][2:4], d==grb_info['NAME'][4:6], t=grb_info['Time [UT]']),
                              format='isot', scale='ut')
        grb_info['MJD'] = astrotime.mjd
        grb_info['MET'] = ConvertMjdToMet(grb_info['MJD'])
        grb_info['RA'] = grb['BAT RA (J2000)']
        grb_info['DEC'] = grb['BAT Dec (J2000)']
        return grb_info
        

class GBMCatalogue(AstropyTableCatalogue):
    def __init__(self, name, file_path):
        AstropyTableCatalogue.__init__(self, name=name, file_path=file_path)

        
    def find_grb(self, grbname):
        grb_info = {}
        grb = self.catalogue_table[self.catalogue_table['NAME']==grbname][0]
        grb_info['NAME'] = grb['NAME']
        astrotime = Time('20{y}-{m}-{d}T{t}'.format(y=grb_info['NAME'][:2], m=grb_info['NAME'][2:4], d==grb_info['NAME'][4:6], t=grb_info['Time [UT]']),
                              format='isot', scale='ut')
        grb_info['MJD'] = grb['TRIGGER_TIME']
        grb_info['MET'] = ConvertMjdToMet(grb_info['MJD'])
        grb_info['RA'] = grb['RA']
        grb_info['DEC'] = grb['DEC']
        grb_info['ERROR'] = grb['ERROR_RADIUS']
        return grb_info

    
class XMLCatalogue(Catalogue):
    def __init__(self, name, file_path):
        Catalogue.__init__(self, name=name, file_path=file_path)
        if self.file_ext.lower()!='xml':
            logger.error('File extenstion must be xml!!')
        fileList = ET.parse(self.file_path)
        self.catalogue_table = fileList.getroot()


    def find_grb(self, grbname):
        grb_info = {}
        for grb in self.catalogue_table:
            if grb.findtext("./GRBNAME")==grbnme:
                grb_info['NAME'] = grb.findtext("./GRBNAME")
                grb_info['MET'] = float(grb.findtext("./MET"))
                grb_info['MJD'] = ConvertMetToMjd(grb_info['MET'])
                if grb.findtext("./ERROR") == "--" or grb.findtext("./ERROR") == "":
                    if grb.findtext("./LATERROR") == "--" or grb.findtext("./LATERROR") == "":
                        grb_info['ERROR'] = np.nan
                    else:
                        grb_info['ERROR'] = float(grb.findtext("./LATERROR"))
                else:
                    if grb.findtext("./LATERROR") == "--" or grb.findtext("./LATERROR") == "" or float(grb.findtext("./ERROR"))<=float(grb.findtext("./LATERROR")):
                        grb_info['ERROR'] = float(grb.findtext("./ERROR"))                    
                        grb_info['RA'] = float(grb.findtext("./RA"))
                        grb_info['DEC'] = float(grb.findtext("./DEC"))
                    else :
                        grb_info['ERROR'] = float(grb.findtext("./LATERROR"))                    
                        grb_info['RA'] = float(grb.findtext("./LATRA"))
                        grb_info['DEC'] = float(grb.findtext("./LATDEC"))
        return grb_info

    
            
