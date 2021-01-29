#!/usr/bin/env python

import sys
import os
import os.path
import ROOT
from ROOT import TTree
from ROOT import TChain
import numpy as np
from astropy.io import fits
from bitarray import bitarray
ROOT.gROOT.SetBatch()


par = sys.argv[1:]

DICT_PATH_LUT_S18ZDIR020catTwoZDIR060_E28binx_Cth40bins = {4096: '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R100_perf.root',
							   8192: '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R030_perf.root',
							   16384: '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R010_perf.root'}


def table_srcon_fits(name_src='PSRJ0218', path_src='/nfs/farm/g/glast/u/mtakahas/data/AllSky/Pulsar/PlotsAllSky_TwoPSRs_2008.root', path_evt='/nfs/farm/g/glast/u/mtakahas/data/AllSky/events_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/trAllSkyMap_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060_2008.root', path_outdir='/nfs/farm/g/glast/u/mtakahas/data/AllSky/Pulsar/fits', dict_path_lut=DICT_PATH_LUT_S18ZDIR020catTwoZDIR060_E28binx_Cth40bins):

    if not os.path.isdir(path_outdir):
        os.path.makedirs(path_outdir)

    file_src = ROOT.TFile(path_src)
    print 'Searching for Directory {0}...'.format(name_src)
    dir_src = file_src.GetDirectory(name_src)
    print dir_src
    print 'Directory {0} found.'.format(dir_src.GetName())    
    name_tree_src = 'trFriend{0}'.format(name_src)
    tree_src = dir_src.Get(name_tree_src)
    print 'Object {0} found.'.format(tree_src.GetName())

    file_evt = ROOT.TFile(path_evt)
    tree_evt = file_evt.Get("EVENTS")
    print '{0} found.'.format(tree_evt.GetName())

    if tree_src.GetEntries()!=tree_evt.GetEntries():
        print "Event number mismatches!!"
        return 1
    else:
        print "Event number: {0}".format(tree_src.GetEntries())
    tree_evt.AddFriend(tree_src)

    print "ON CalOnly events: {0}".format(tree_evt.GetEntries("FLAG_ON==1 && Category==2 && CosTHETA>=0.2"))

    dict_file_lut = {}
    dict_hist_lut_psfq68 = {}
    dict_hist_lut_psfq95 = {}
    for nclass, path_lut in dict_path_lut.items():
        dict_file_lut[nclass] = ROOT.TFile(path_lut, "READ")
        print dict_file_lut[nclass].GetName()
        dict_hist_lut_psfq68[nclass] = dict_file_lut[nclass].Get('psf_cth_q68_hist')
        print '{0} found.'.format(dict_hist_lut_psfq68[nclass].GetName())
        dict_hist_lut_psfq95[nclass] = dict_file_lut[nclass].Get('psf_cth_q95_hist')
        print '{0} found.'.format(dict_hist_lut_psfq95[nclass].GetName())

    dict_list_q = {'ENERGY': [],
                   'THETA': [],
                   'EVENT_CLASS': [],
                   'EVENT_CLASS_LABEL': [],
                   'EVENT_TYPE': [],
                   'EVENT_ID': [],
                   'RUN_ID': [],
                   'RA': [],
                   'DEC': [],
                   'L': [],
                   'B': [],
                   'ZENITH_ANGLE': [],
                   'TIME': [],
                   'PSF_Q68': [],
                   'PSF_Q95': []
               }

    DICT_EVENT_CLASS_32X = {4096: np.array(19*[False]+[True]+12*[False]),
                            8192: np.array(18*[False]+2*[True]+12*[False]),
                            16384: np.array(17*[False]+3*[True]+12*[False])
}

    DICT_EVENT_CLASS_LABEL = {4096: 'CalOnly_R100',
                              8192: 'CalOnly_R030',
                              16384: 'CalOnly_R010'
}

    tbhdr = fits.Header()
    tbhdr.set("TELESCOP", "GLAST", "name of telescope generating data")
    tbhdr.set("INSTRUME", "LAT", "name of instrument generating data")
    tbhdr.set("OBSERVER", "Peter Michelson", "GLAST/LAT PI")
    tbhdr.set("ORIGIN", "Mitsunari Takahashi", "CalOnly analysis developer")
    tbhdr.set("HDUCLAS1", "EVENTS", "extension contains events")
    tbhdr.set("PASS_VER", "P8R2", "IRF pass version corresponding to a specific se")
    tbhdr.set("TIMEUNIT", "s", "units for the time related keywords")
    tbhdr.set("TSTART", "239557417.0", "mission time of the start of the observation")
    tbhdr.set("TSTOP", "504919671.0", "mission time of the end of the observation")
#    tbhdr.set("", "", "")
    tbhdr.set("TUNIT1", "MeV", "physical unit of field")
    tbhdr.set("TLIMIN1", "0.", "minimum value")
    tbhdr.set("TLIMAX1", "10000000.", "maximum value")
    tbhdr.set("TUNIT2", "deg", "physical unit of field")
    tbhdr.set("TLIMIN2", "0.", "minimum value")
    tbhdr.set("TLIMAX2", "360.", "maximum value")
    tbhdr.set("TUNIT3", "deg", "physical unit of field")
    tbhdr.set("TLIMIN3", "-90.", "minimum value")
    tbhdr.set("TLIMAX3", "90.", "maximum value")
    tbhdr.set("TUNIT4", "deg", "physical unit of field")
    tbhdr.set("TLIMIN4", "0.", "minimum value")
    tbhdr.set("TLIMAX4", "360.", "maximum value")
    tbhdr.set("TUNIT5", "deg", "physical unit of field")
    tbhdr.set("TLIMIN5", "-90.", "minimum value")
    tbhdr.set("TLIMAX5", "90.", "maximum value")
    tbhdr.set("TUNIT6", "deg", "physical unit of field")
    tbhdr.set("TLIMIN6", "0.", "minimum value")
    tbhdr.set("TLIMAX6", "180.", "maximum value")
    tbhdr.set("TUNIT7", "deg", "physical unit of field")
    tbhdr.set("TLIMIN7", "0.", "minimum value")
    tbhdr.set("TLIMAX7", "180.", "maximum value")
    tbhdr.set("TUNIT8", "s", "physical unit of field")
    tbhdr.set("TLIMIN8", "0.", "minimum value")
    tbhdr.set("TLIMAX8", "10000000000.", "maximum value")
    tbhdr.set("NDSKEYS", "1", "Number of data subspace keywords in header")#    tbhdr.set("NDSKEYS", "5", "Number of data subspace keywords in header")
    # tbhdr.set("DSTYP1", "TIME")
    # tbhdr.set("DSUNI1", "s")
    # tbhdr.set("DSVAL1", "TABLE")
    # tbhdr.set("DSREF1", ":GTI")
    # tbhdr.set("DSTYP2", "BIT_MASK(EVENT_CLASS,4096,P8R2)")
    # tbhdr.set("DSUNI2", "DIMENSIONLESS")
    # tbhdr.set("DSVAL2", "1:1")
    # tbhdr.set("DSTYP3", "POS(RA,DEC)")
    # tbhdr.set("DSUNI3", "deg")
    # tbhdr.set("DSVAL3", "circle(0.0,0.0,180.000000)")
    # tbhdr.set("DSTYP4", "TIME")
    # tbhdr.set("DSUNI4", "s")
    # tbhdr.set("DSVAL4", "239557417.0:504919671.0")
    tbhdr.set("DSTYP1", "ENERGY")
    tbhdr.set("DSUNI1", "MeV")
    tbhdr.set("DSVAL1", "22387.211385:562341.325190")
    tbhdr['COMMENT'] = "List of CalOnly gamma-like events around {0} with PSF_Q68 and PSF_Q95. The event cut is S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15.".format(name_src)

    tbhdr2 = fits.Header()
    tbhdr2.set("TIMEUNIT", "s", "units for the time related keywords")
    tbhdr2.set("TLMIN1", "0.0", "minimum value")
    tbhdr2.set("TLMAX1", "1.0E10", "maximum value")
    tbhdr2.set("TLMIN2", "0.0", "minimum value")
    tbhdr2.set("TLMAX2", "1.0E10", "maximum value")

    for evt in tree_evt:
        if evt.FLAG_ON==1 and evt.c==2 and evt.cth>=0.2:
            dict_list_q['ENERGY'].append(10**evt.e)
            dict_list_q['THETA'].append(evt.th)
            dict_list_q['EVENT_CLASS'].append(DICT_EVENT_CLASS_32X[evt.s])
            dict_list_q['EVENT_CLASS_LABEL'].append(DICT_EVENT_CLASS_LABEL[evt.s])
#            dict_list_q['EVENT_CLASS'].append(np.fromstring(np.binary_repr(evt.s), dtype='S1').astype(int)) ##evt.s)
#            dict_list_q['EVENT_CLASS'].append(int(np.binary_repr(evt.s)))
            dict_list_q['EVENT_TYPE'].append(evt.ty)
            dict_list_q['EVENT_ID'].append(evt.evid)
            dict_list_q['RUN_ID'].append(evt.run)
            dict_list_q['RA'].append(evt.ra)
            dict_list_q['DEC'].append(evt.dec)
            dict_list_q['L'].append(evt.l)
            dict_list_q['B'].append(evt.b)
            dict_list_q['ZENITH_ANGLE'].append(evt.z)
            dict_list_q['TIME'].append(evt.t)
            dict_list_q['PSF_Q68'].append(dict_hist_lut_psfq68[evt.s].Interpolate(evt.e, evt.cth))
            dict_list_q['PSF_Q95'].append(dict_hist_lut_psfq95[evt.s].Interpolate(evt.e, evt.cth))

    #print dict_list_q['EVENT_CLASS']
    #print np.array(dict_list_q['EVENT_CLASS'])
    coldefs = fits.ColDefs([fits.Column(name='ENERGY', format='E', array=np.array(dict_list_q['ENERGY'])),
                            fits.Column(name='RA', format='E', array=np.array(dict_list_q['RA']), unit='degree'),
                            fits.Column(name='DEC', format='E', array=np.array(dict_list_q['DEC']), unit='degree'),
                            fits.Column(name='L', format='E', array=np.array(dict_list_q['L']), unit='degree'),
                            fits.Column(name='B', format='E', array=np.array(dict_list_q['B']), unit='degree'),
                            fits.Column(name='THETA', format='E', array=np.array(dict_list_q['THETA']), unit='degree'),
                            fits.Column(name='ZENITH_ANGLE', format='E', array=np.array(dict_list_q['ZENITH_ANGLE']), unit='degree'),
                            fits.Column(name='TIME', format='D', array=np.array(dict_list_q['TIME']), unit='s'),
                            fits.Column(name='EVENT_CLASS', format='32X', array=np.array(dict_list_q['EVENT_CLASS'])),
                            fits.Column(name='EVENT_CLASS_LABEL', format='12A', array=np.array(dict_list_q['EVENT_CLASS_LABEL'])),
                            fits.Column(name='EVENT_TYPE', format='32X', array=np.array(dict_list_q['EVENT_TYPE'])),
                            fits.Column(name='EVENT_ID', format='J', array=np.array(dict_list_q['EVENT_ID'])),
                            fits.Column(name='RUN_ID', format='J', array=np.array(dict_list_q['RUN_ID'])),
                            fits.Column(name='PSF_Q68', format='E', array=np.array(dict_list_q['PSF_Q68']), unit='degree'),
                            fits.Column(name='PSF_Q95', format='E', array=np.array(dict_list_q['PSF_Q95']), unit='degree')])

    tbhdu = fits.BinTableHDU.from_columns(coldefs, tbhdr)
    #hdu = fits.BinTableHDU.from_columns(coldefs)
    print tbhdu.data['EVENT_CLASS'][0]
    print tbhdu.data['EVENT_CLASS_LABEL'][0]
    tbhdu.name = "EVENTS"

    coldefs2 = fits.ColDefs([fits.Column(name='START', format='D', array=np.array([239557417])),
                             fits.Column(name='STOP', format='D', array=np.array([504919671]))])
    tbhdu2 = fits.BinTableHDU.from_columns(coldefs2, tbhdr2)
    tbhdu2.name = "GTI"

    empty_primary = fits.PrimaryHDU() #header=hdr)
    hdulist = fits.HDUList([empty_primary, tbhdu, tbhdu2])

    name_fits = os.path.basename(path_evt).replace('.root', '.fits').replace('trAllSkyMap', name_src)
    path_fits = '/'.join([path_outdir, name_fits])
    #tbhdu.writeto(path_fits, overwrite=True)
    hdulist.writeto(path_fits)


    # hdulist = fits.open(path_fitsfile)
    # tbdata = hdulist[1].data
    #     tbdata.field('ENERGY')
    #     tbdata.field('THETA')
    #     tbdata.field('EVENT_CLASS')
    #     tbdata.field('EVENT_TYPE')


if __name__=='__main__':
    for year in [2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]:
        table_srcon_fits(name_src=par[0], path_src='/nfs/farm/g/glast/u/mtakahas/data/AllSky/Pulsar/PlotsAllSky_TwoPSRs_{0:.0f}.root'.format(year), path_evt='/nfs/farm/g/glast/u/mtakahas/data/AllSky/events_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/trAllSkyMap_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060_{0:.0f}.root'.format(year), path_outdir='/nfs/farm/g/glast/u/mtakahas/data/AllSky/Pulsar/fits/{0}'.format(par[0]))
        
