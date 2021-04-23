#!/usr/bin/env python2

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
        self.NAME = str(name)
        self.path_file = path_file
        self.array = array('f',[0.0])


    def get_value(self):
        return self.array[0]


dict_var_bdt_default = {'ZDirDown020100': Variable('G01_S050B050GO08to09_catFourRWZDirRWCALE_16_BDTG5e2NEvMi250'), #'S20_S050B050_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIRcatFourRWZDirRWCALE_16_UpDownZDir_NTree1000_BDTG1000D06NEvtMin250'),
                        'ZDirDown060100': Variable('S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06'),
                        'ZDirDown020060': Variable('S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06'),
                        'ZDirUp020050': Variable('T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR020050_BDTG500D06'),
                        'ZDirUp050100': Variable('T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR050100_BDTG500D06')}

dict_var_bdt_evod = {'ZDirDown020100ev': Variable('G01_S050B050GO08to20AcdSigma0ex3FHL_sepSixCross714595EvOdRWZDirRWCALE_16_RealAda_PruCC60_od_RABD10nC70PruCC60'), #Trained by FT1EventId%2==1 -> Evaluate FT1EventId%2=0
                     'ZDirDown020100od': Variable('G01_S050B050GO08to20AcdSigma0ex3FHL_sepSixCross714595EvOdRWZDirRWCALE_16_RealAda_PruCC60_ev_RABD10nC70PruCC60')} #Trained by FT1EventId%2==0 -> Evaluate FT1EventId%2=1

var_updown_default = Variable('U01_V200909_62IRFBK_020RAWE20woRW_09ZDir_MinNode2_D06_MinNode00005')
        

def convert_bdt_to_fill(bdt, conversion='log'):
    if conversion=='antilog':
        antilogbdt = 1.-(1.+bdt.get_value())/2.
        if antilogbdt<=0:
            print "BDT = {0}".format(antilogbdt)
            if antilogbdt==0:
                logbdt = 10
            else:
                logbdt = -log10(antilogbdt)
        return logbdt
    elif conversion=='x5':
        return bdt.get_value()*5.
    else:
        return bdt.get_values()


def get_healpix_pixel(ra_rad, dec_rad, sign_arrival, nhpside):
    phi_unsign = ra_rad
    theta_unsign = pi/2.-dec_rad
    phi_sign = phi_unsign if sign_arrival>0 else phi_unsign - pi
    theta_sign = theta_unsign if sign_arrival>0 else pi - theta_unsign
    return hppf.ang2pix(nside=nhpside, theta=theta_sign, phi=phi_sign)


def fill_flighteventBDT(hist3D, tree, list_healpxs, nhpside, zmax, dict_var_bdt=dict_var_bdt_default, var_updown=var_updown_default, met_start=247017601., met_stop=sys.maxint, evod=None, bdt_conversion='log'):
    """hist3D: X) log10(Energy) Y) BDT cut Z) Inclination
"""
    print 'Time range: MET {0}-{1}'.format(met_start, met_stop)

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

        kpix = get_healpix_pixel(radians(ev.FT1CalRa), radians(ev.FT1CalDec), sign_arrival=sign_arrival, nhpside=nhpside)
        # ra_rad = radians(ev.FT1CalRa)
        # dec_rad = radians(ev.FT1CalDec)
        # phi_unsign = ra_rad
        # theta_unsign = pi/2.-dec_rad
        # phi_sign = phi_unsign if sign_arrival>0 else phi_unsign - pi
        # theta_sign = theta_unsign if sign_arrival>0 else pi - theta_unsign
        # kpix = hppf.ang2pix(nside=nhpside, theta=theta_sign, phi=phi_sign)

        if sign_arrival>0 and ev.Cal1MomZCrossSide840>=0 and ev.FT1CalZenithTheta<=zmax and (-1.56245e+04+3.80384e+03*log10(ev.WP8CalOnlyEnergy)-1.28404e+03*ev.Cal1MomZDir>=100):
            if evod!=None:
                if ev.FT1EventId%2==0:
                    bdt = dict_var_bdt['ZDirDown020100ev']
                elif ev.FT1EventId%2==1:
                    bdt = dict_var_bdt['ZDirDown020100od']
                else:
                    print "{0} is NOT even or odd!!".format(ev.FT1EventId)
                    sys.exit(1)
            else:
                bdt = dict_var_bdt['ZDirDown020100']
                # if ev.Cal1MomZDir>=0.6:
                #     bdt = dict_var_bdt_default['ZDirDown060100']
                # elif ev.Cal1MomZDir>=0.2:
                #     bdt = dict_var_bdt_default['ZDirDown020060']
        elif sign_arrival<0 and (180.-ev.FT1CalZenithTheta)<=zmax:
            if ev.Cal1MomZDir>=0.5:
                bdt = dict_var_bdt['ZDirUp050100']
            elif ev.Cal1MomZDir>=0.2:
                bdt = dict_var_bdt['ZDirUp020050']

        if kpix in list_healpxs:
            if bdt is not None:
                if not (ev.WP8CalOnlyEnergy>=10**4.35 and ev.WP8CalOnlyEnergy<=10**5.75):
                    continue
                # antilogbdt = 1.-(1.+bdt.get_value())/2.
                # bdt_x5 = 5.*bdt
                # if antilogbdt<=0:
                #     print "BDT = {0}".format(antilogbdt)
                #     if antilogbdt==0:
                #         logbdt = 100
                #     continue
                # else:
                #     logbdt = -log10(antilogbdt)
                fillbdt = convert_bdt_to_fill(bdt, bdt_conversion)
                hist3D.Fill(log10(ev.WP8CalOnlyEnergy), fillbdt, ev.Cal1MomZDir*sign_arrival)
                if iev%10000==0:
                    print "(l, b) = {s}({l}, {b}); Z = ({s}){z}; log10(BDT) = {t}; {p}% done.".format(s=sign_arrival, l=ev.FT1CalL, b=ev.FT1CalB, z=ev.FT1CalZenithTheta, t=fillbdt, p=int(iev*100/tree.GetEntries()))


@click.command()
@click.argument('meritpath', type=str)
@click.option('--suffix', '-s', default='')
@click.option('--outdir', '-o', default='.')
@click.option('--metmin', default=247017601, type=float)
@click.option('--metmax', default=613872005, type=float)
@click.option('--zmax', '-z', default=90, type=float)
@click.option('--nhpside', default=32, type=int)
@click.option('--bdt', '-b', type=str, default=None)
@click.option('--friends', '-f', multiple=True, type=str, default=None, help='Path of your friend MeritTuple file')
@click.option('--bdtcth', '-c', default=None, type=click.Choice([None, 'ZDirDown020100', 'ZDirDown060100', 'ZDirDown020060', 'ZDirUp020050', 'ZDirUp050100']))
@click.option('--nobk', is_flag=True, default=False, help='Use when ignoring input events come from the backside.')
@click.option('--evod', type=(str, str), default=None, help='Use when different variables for even/odd events.')
#@click.option('--evod', is_flag=True, default=False, help='Use when different variables for even/odd events.')
@click.option('--conversion', type=str, default='log')
@click.option('--bsub', '-b', is_flag=True, default=False)
def main(meritpath, suffix, outdir, metmin, metmax, zmax, nhpside, bdt, friends, bdtcth, nobk, evod, conversion, bsub):
    print "Please check the variables in this script!!"
    list_path_filein = ls_list(meritpath)
    if len(suffix)>0:
        if suffix[0]!="_":
            suffix = '_'+suffix

    for path_filein in list_path_filein:
        name_log = os.path.basename(path_filein).replace("LPA_", "Hist3D_").replace(".root","{0}.log".format(suffix))
        path_dir_log = '{d}/logs'.format(d=outdir)
        if not os.path.isdir(path_dir_log):
            os.makedirs(path_dir_log)
        path_log = '{d}/{n}'.format(d=path_dir_log, n=name_log)
        if bsub is True:
            cmd = ['bsub', '-o'+path_log, '-JFill'+name_log[-7:-4], '-W600', 'python2', '/nfs/farm/g/glast/u/mtakahas/PythonModuleMine/Fermi/FillGalOffEventBDT.py', path_filein, '--suffix', suffix, '--outdir', outdir, '--metmin', str(metmin), '--metmax', str(metmax), '--zmax', str(zmax), '--nhpside', str(nhpside), '--conversion', conversion]
            if bdt!=None:
                cmd += ['--bdt', bdt]
            if bdtcth!=None:
                cmd += ['--bdtcth', bdtcth]
            if friends!=None:
                for fr in friends:
                    cmd += ['--friends', fr]
            if nobk==True:
                cmd += ['--nobk']
            if evod!=None:
                cmd += ['--evod', evod[0], evod[1]]
            print cmd
            subprocess.call(cmd)
        else:
            print 'Input file: {}'.format(path_filein)
            # HEALPix
            pathCatalogue = "/nfs/farm/g/glast/u/mtakahas/data/catalogue/gll_psch_v13.fit" #"/disk/gamma/cta/store/takhsm/FermiData/catalogue/gll_psch_v09.fit"
            NHPSIDE_OFF = nhpside #32 #16
            list_galoff_healpxs = find_galoff_healpxs(NHPSIDE_OFF, pathCatalogue, dict_flux_cutangle={2.5E-10:4., 5E-9:10.})
            #list_galoff_healpxs = find_galoff_healpxs(NHPSIDE_OFF, 0, pathCatalogue)
            print 'Galactic OFF HEALPix pixels:\n{}'.format(list_galoff_healpxs)

            # Input
            filetree = TFile(path_filein, "READ")
            print 'File {} is opened.'.format(filetree.GetName())
            filetree.cd()
            treeevt = filetree.Get("MeritTuple")
            print 'Tree {} is found.'.format(treeevt.GetName())

            # Friend tree
            if friends!=None:
                print friends
                for friend in friends:
                    treeevt.AddFriend('MeritTuple', friend)

           # Output
            name_fileout = os.path.basename(meritpath).replace("LPA_", "Hist3D_").replace('.root', '{0}.root'.format(suffix))
            path_fileout = '/'.join([outdir, name_fileout])
            fileout = TFile(path_fileout, "RECREATE")
            print "Output file {} is recreated".format(path_fileout)
            fileout.cd()

          # Event filling
            hist = TH3F("Events",  "Events;log_{10}E;BDT;#cos(#theta)", 28, 4.35, 5.75, 2000, 0, 5, 50 if nobk==True else 100, 0 if nobk==True else -1., 1.)
            if bdtcth!=None:
                dict_var = {bdtcth: dict_var_bdt_default[bdtcth]} 
            else: 
                if evod==None:
                    dict_var = dict_var_bdt_default
                else:
                    dict_var = {'ZDirDown020100ev': Variable(evod[0]), #Trained by FT1EventId%2==1 -> Evaluate FT1EventId%2=0
                     'ZDirDown020100od': Variable(evod[1])} #Trained by FT1EventId%2==0 -> Evaluate FT1EventId%2=1 #dict_var_bdt_evod
            
            fill_flighteventBDT(hist3D=hist, tree=treeevt, list_healpxs=list_galoff_healpxs, nhpside=NHPSIDE_OFF, met_start=metmin, met_stop=metmax, zmax=zmax, 
                                dict_var_bdt=dict_var, 
                                var_updown=None if nobk==True else var_updown_default, evod=evod, bdt_conversion=conversion)
            hist.Write()


if __name__ == '__main__':
    main()

