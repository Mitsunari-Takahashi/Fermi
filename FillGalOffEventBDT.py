#!/usr/bin/env python

import sys
import os.path
import ROOT
from ROOT import TFile, TTree, TChain, TH1, TH2, TH3, TH1F, TH2F, TH3F
import numpy as np
import commands
import click
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, log, log10
import subprocess
import healpy as hp
from healpy import pixelfunc as hppf
from pConvertSkyCoords import *
from pCalcAngDist_Catalogue import *
from pFindHEALPix import find_galoff_healpxs
from pLsList import ls_list

ROOT.gROOT.SetBatch()

class Variable:
    def __init__(self, name, path_file=None):
        self.NAME = name
        self.path_file = path_file
        self.array = array('f',[0.0])


    def get_value(self):
        return self.array[0]


dict_var_bdt_default = {'ZDirDown020100': Variable('G01_S050B050GO08to09_catFourRWZDirRWCALE_16_BDTG5e2NEvMi250'), #'S20_S050B050_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIRcatFourRWZDirRWCALE_16_UpDownZDir_NTree1000_BDTG1000D06NEvtMin250'),
                        'ZDirDown060100': Variable('S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06'),
                        'ZDirDown020060': Variable('S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06'),
                        'ZDirUp020050': Variable('T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR020050_BDTG500D06'),
                        'ZDirUp050100': Variable('T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR050100_BDTG500D06')}

var_updown_default = Variable('U01_V200909_62IRFBK_020RAWE20woRW_09ZDir_MinNode2_D06_MinNode00005')
        

def fill_flighteventBDT(hist3D, tree, list_healpxs, nhpside, zmax, dict_var_bdt=dict_var_bdt_default, var_updown=var_updown_default, met_start=247017601., met_stop=sys.maxint):
    """hist3D: X) log10(Energy) Y) BDT cut Z) Inclination
"""
    # Setting the branch addresses
    for i,v in enumerate(dict_var_bdt.values()):
        tree.SetBranchAddress(v.NAME, v.array)
    if var_updown!=None:
        tree.SetBranchAddress(var_updown.NAME, var_updown.array)

    # Event loop
    for iev, ev in enumerate(tree):
        if ev.EvtElapsedTime<met_start or ev.EvtElapsedTime>=met_stop:
            continue
        bdt = None
        sign_arrival = 1 if var_updown==None or var_updown.get_value()>=0 else -1

        ra_rad = radians(ev.FT1CalRa)
        dec_rad = radians(ev.FT1CalDec)
        phi_unsign = ra_rad
        theta_unsign = pi/2.-dec_rad
        phi_sign = phi_unsign if sign_arrival>0 else phi_unsign - pi
        theta_sign = theta_unsign if sign_arrival>0 else pi - theta_unsign
        kpix = hppf.ang2pix(nside=nhpside, theta=theta_sign, phi=phi_sign)

        if kpix in list_healpxs:
            if sign_arrival>0 and ev.Cal1MomZCrossSide840>=0 and ev.FT1CalZenithTheta<=zmax:
                bdt = dict_var_bdt_default['ZDirDown020100']
                # if ev.Cal1MomZDir>=0.6:
                #     bdt = dict_var_bdt_default['ZDirDown060100']
                # elif ev.Cal1MomZDir>=0.2:
                #     bdt = dict_var_bdt_default['ZDirDown020060']
            elif sign_arrival<0 and (180.-ev.FT1CalZenithTheta)<=zmax:
                    if ev.Cal1MomZDir>=0.5:
                        bdt = dict_var_bdt_default['ZDirUp050100']
                    elif ev.Cal1MomZDir>=0.2:
                        bdt = dict_var_bdt_default['ZDirUp020050']
            if bdt is not None:
                if not (ev.WP8CalOnlyEnergy>=10**4.35 and ev.WP8CalOnlyEnergy<=10**5.75):
                    #print "WP8CalOnlyEnergy = {0}".format(ev.WP8CalOnlyEnergy)
                    continue
                antilogbdt = 1.-(1.+bdt.get_value())/2.
                if antilogbdt<=0:
                    print "BDT = {0}".format(antilogbdt)
                    if antilogbdt==0:
                        logbdt = 100
                    continue
                else:
                    logbdt = -log10(antilogbdt)
                hist3D.Fill(log10(ev.WP8CalOnlyEnergy), logbdt, ev.Cal1MomZDir*sign_arrival)
                if iev%5000==0:
                    print "(l, b) = {s}({l}, {b}); Z = ({s}){z}; log10(BDT) = {t}; {p}% done.".format(s=sign_arrival, l=ev.FT1CalL, b=ev.FT1CalB, z=ev.FT1CalZenithTheta, t=logbdt, p=int(iev*100/tree.GetEntries()))


@click.command()
@click.argument('meritpath', type=str)
@click.option('--suffix', '-s', default='')
@click.option('--outdir', '-o', default='.')
@click.option('--zmax', '-z', default=90, type=float)
@click.option('--bdt', '-o', type=str, default=None)
@click.option('--bdtcth', '-c', default=None, type=click.Choice([None, 'ZDirDown020100', 'ZDirDown060100', 'ZDirDown020060', 'ZDirUp020050', 'ZDirUp050100']))
@click.option('--nobk', is_flag=True, default=False, help='Use when ignoring input events come from the backside.')
@click.option('--bsub', '-b', is_flag=True, default=False)
def main(meritpath, suffix, outdir, zmax, bdt, bdtcth, nobk, bsub):
    list_path_filein = ls_list(meritpath)

    for path_filein in list_path_filein:
        name_log = os.path.basename(path_filein).replace("LPA_", "Hist3D_").replace(".root",".log")
        path_dir_log = '{d}/logs'.format(d=outdir)
        if not os.path.isdir(path_dir_log):
            os.makedirs(path_dir_log)
        path_log = '{d}/{n}'.format(d=path_dir_log, n=name_log)
        if bsub is True:
            cmd = ['bsub', '-o'+path_log, '-JFill'+name_log[-7:-4], '-W300', 'python', '/nfs/farm/g/glast/u/mtakahas/PythonModuleMine/Fermi/FillGalOffEventBDT.py', path_filein, '--suffix', suffix, '--outdir', outdir, '--zmax', str(zmax)]
            if bdt!=None:
                cmd += ['--bdt', bdt]
            if bdtcth!=None:
                cmd += ['--bdtcth', bdtcth]
            if nobk==True:
                cmd += ['--nobk']
            print cmd
            subprocess.call(cmd)
        else:
            print 'Input file: {}'.format(path_filein)
            # HEALPix
            pathCatalogue = "/nfs/farm/g/glast/u/mtakahas/data/catalogue/gll_psch_v13.fit" #"/disk/gamma/cta/store/takhsm/FermiData/catalogue/gll_psch_v09.fit"
            NHPSIDE_OFF = 32 #16
            list_galoff_healpxs = find_galoff_healpxs(NHPSIDE_OFF, 0, pathCatalogue)
            print 'Galactic OFF HEALPix pixels:\n{}'.format(list_galoff_healpxs)

            # Input
            filetree = TFile(path_filein, "READ")
            print 'File {} is opened.'.format(filetree.GetName())
            filetree.cd()
            treeevt = filetree.Get("MeritTuple")
            print 'Tree {} is found.'.format(treeevt.GetName())

           # Output
            name_fileout = os.path.basename(meritpath).replace("LPA_", "Hist3D_")
            path_fileout = '/'.join([outdir, name_fileout])
            fileout = TFile(path_fileout, "RECREATE")
            print "Output file {} is recreated".format(path_fileout)
            fileout.cd()

          # Event filling
            hist = TH3F("Events",  "Events;log_{10}E;BDT;#cos(#theta)", 28, 4.35, 5.75, 2000, 0, 5, 50 if nobk==True else 100, 0 if nobk==True else -1., 1.)
            fill_flighteventBDT(hist3D=hist, tree=treeevt, list_healpxs=list_galoff_healpxs, nhpside=NHPSIDE_OFF, zmax=zmax, 
                                dict_var_bdt={bdtcth: dict_var_bdt_default[bdtcth]} if bdtcth!=None else dict_var_bdt_default, 
                                var_updown=None if nobk==True else var_updown_default)
            hist.Write()


if __name__ == '__main__':
    main()

