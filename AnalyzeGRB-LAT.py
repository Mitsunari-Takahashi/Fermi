#!/usr/bin/env python

import os
import sys
import ROOT
from ROOT import TTree, TChain, TH1, TH2, TH3
import numpy as np
import yaml
import xml.etree.ElementTree as ET
#import pandas
import datetime
from ctypes import *
import click
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, log10
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

#from pCutBDT import cutBDT
from pAnalysisConfig import *
from pFindGrbInfo import *
from pLsList import ls_list


def shifttimebox(lst_box, newtime, maxnum=1):
    if len(lst_box)==maxnum:
        for iwait in range(maxnum-1):
            lst_box[iwait] = lst_box[iwait+1]
        lst_box[maxnum-1] = newtime
    else:
        print 'Wrong array length!!!'


@click.command()
@click.option('--suffix', '-s',default="", help="Suffix for name of output product. GRB name is added automatically..")
@click.argument('grbid', type=str)
@click.argument('evtfiles', type=str)
#@click.argument('datafiles', type=str)
@click.argument('start', type=float)
@click.argument('stop', type=float)
@click.option('--fixpsfenergy', '-e', type=float, default=0.0, help="Set energy in log scale if you will fix the PSF cut on energy.")
@click.option('--fixpsfinclin', '-i', type=float, default=0.0, help="Set cos(inclination angle) if you will fix the PSF cut on inclination.")
@click.option('--exclude', type=(float, float), default=(float(sys.maxsize),0.0), help="This time domain is excluded from the data which is analyzed. Assign the start and stop in MET.")
@click.option('--thresholdenergywaitingtime', type=float, default=0.0, help="Set energy threshold for waiting time plot in log scale.")
@click.option('--waitingnumber', '-w', type=int, default=1, help="Set number of consective events you wait.")
def main(grbid, evtfiles, start, stop, suffix, fixpsfenergy, fixpsfinclin, exclude, thresholdenergywaitingtime, waitingnumber):
    # ----- Event class setup -----
    cfg = ClassConfig('Both', [10, 3, 1, 0.3], 1)
    aCutEGB = cfg.aCutEGB
    aaStrSelect = cfg.aaStrSelect
    aStrSelect = cfg.aStrSelect

    #IRF
    listPathFilePerf = [['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                        ['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_R100_perf.root', 
                         '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_R030_perf.root', 
                         '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_R010_perf.root',
                         '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_R003_perf.root']]
#                        ['/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root', '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root']]
    htgPerf = CutPerformanceHtg(listPathFilePerf)

    trCatalogue = ROOT.TTree("trGRB", "GBM burst catalogue")
    # Target
    listTgtGRB = [grbid] #[par[1]]

    # Catalogue Data
    pathList = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/PublicTableGRBs.xml" #"/nfs/farm/g/glast/u/mtakahas/data/catalogue/PublicTableGRBs.xml"
    fileList = ET.parse(pathList)
    rtXml = fileList.getroot()

    print "===================="
    #nFile = (len(par)-5)/2
    listFileIn = ls_list(evtfiles) #par[6:6+nFile]
    #listFileDat = ls_list(datafiles) #par[6+nFile:6+2*nFile]
    #print listFileIn

    aliasSelections = yaml.load(open("{0}/config/pass8_event_selections.yaml".format(os.environ.get("EVENTSELECT")),'r'))

    chInput = ROOT.TChain('EVENTS')
    for fileIn in listFileIn:
        chInput.Add(fileIn)
    if suffix!="":
        suffix = "_" + suffix
    aCutPsf = [95, 68]

    for nameGrb in listTgtGRB:
        fileOut = ROOT.TFile('Plot_GRB{0}{1}.root'.format(nameGrb, suffix), 'UPDATE')
        dict_grb = find_grb_info(nameGrb, rtXml)
        if float(stop)>0 and float(start)<0:
            metStart = dict_grb["TRIGGER_MET"]+float(start)
            metStop = dict_grb["TRIGGER_MET"]+float(stop)
        elif float(stop)>0 and float(start)==0:
            metStart = dict_grb["TRIGGER_MET"]
            metStop = dict_grb["TRIGGER_MET"]+float(stop)
        elif float(start)>0 and float(stop)>float(start):
            metStart = float(start)
            metStop = float(stop)
        else:
            metStart = chIn.GetMinimum("t")-1
            metStop = chIn.GetMaximum("t")+1
        chIn = chInput.CopyTree('TIME>={0} && TIME<{1} && (TIME<{2} || TIME>={3})'.format(metStart, metStop, exclude[0], exclude[1]))
        chIn.Write()

        print 'Analysis time domain: MET', metStart, '-', metStop
        if not (exclude[0]==sys.maxsize and exclude[1]==0):
            print 'Excluded time domain: MET', exclude[0], '-', exclude[1]
        nEventChain = chIn.GetEntries()
        #nEventTime = chIn.GetEntries('TIME>={0} && TIME<{1} && (TIME<{2} || TIME>={3})'.format(metStart, metStop, exclude[0], exclude[1]))
        print "Total number of events in the time domain:", nEventChain #nEventTime
        vecTgt = np.array([cos(radians(dict_grb["DEC"]))*cos(radians(dict_grb["RA"])), cos(radians(dict_grb["DEC"]))*sin(radians(dict_grb["RA"])), sin(radians(dict_grb["DEC"]))])

    # TTree
        fileOut.cd()
        trGRB = ROOT.TTree("EVENTS_GRB{0}".format(nameGrb), "Friend TTree for GRB{0}".format(nameGrb))
        cdTimeGRB = c_double()
        trGRB.Branch('TIME_GRB', cdTimeGRB, 'TIME_GRB/D')
        cdAngSep = c_double()
        trGRB.Branch('ANG_SEP', cdAngSep, 'ANG_SEP/D')
        cbFlagPSF68 = c_bool()
        trGRB.Branch('FLAG_PSF68', cbFlagPSF68, 'FLAG_PSF68/O')
        cbFlagPSF95 = c_bool()
        trGRB.Branch('FLAG_PSF95', cbFlagPSF95, 'FLAG_PSF95/O')

        mgr = ROOT.TMultiGraph("mgr", "Gamma-like events around GRB{0}".format(nameGrb))
        greOn = []
        mgrZenith = ROOT.TMultiGraph("mgrZenith", "Zenith angle of events around GRB{0}".format(nameGrb))
        greZenith = []
        htgRADEC = ROOT.TH2D("htgRADEC", "DEC vs. RA of flagged events", 360, 0, 360, 180, -90, 90)

        aHtgEvt = []
        NBIN_CTH = 50
        EDGE_CTH_LOW = 0.0
        EDGE_CTH_UP = 1.0
        NBIN_ZEN = 180
        EDGE_ZEN_LOW = 0
        EDGE_ZEN_UP = 180
        NBIN_ENE = 28
        EDGE_ENE_LOW =  4.35
        EDGE_ENE_UP =  5.75
        if fixpsfenergy!=0 and fixpsfinclin!=0:
            str_fix_psf = "(fixed with log_{10}E={0}, cos#theta={1})".format(fixpsfenergy, fixpsfinclin)
        elif fixpsfenergy!=0:
            str_fix_psf = "(fixed with log_{10}E={0})".format(fixpsfenergy)
        elif fixpsfinclin!=0:
            str_fix_psf = "(fixed with cos#theta={0})".format(fixpsfinclin)
        else:
            str_fix_psf = ""

        aHtgInterval = []
        EDGE_ITV_LOW = 0
        EDGE_ITV_UP = 10000000 #8
        NBIN_ITV = 10000 #100*(EDGE_ITV_UP-EDGE_ITV_LOW)
        EDGE_MET_LOW = metStart
        EDGE_MET_UP = metStop
        NBIN_MET = int((EDGE_MET_UP-EDGE_MET_LOW)/(365.25/12.*86400))
        metPrevious = []

        for cutPsf in aCutPsf:
            print 'PSF cut:', cutPsf, '%'
            greOn.append([])
            greZenith.append([])
            aHtgEvt.append([])
            aHtgInterval.append([])
            metPrevious.append([])
            for pC in range(len(aaStrSelect)):
                greOn[-1].append([])
                greZenith[-1].append([])
                aHtgEvt[-1].append([])
                aHtgInterval[-1].append([])
                metPrevious[-1].append([])
                for qC in range(len(aaStrSelect[pC])):
                    greOn[-1][-1].append(ROOT.TGraphErrors())
                    greOn[-1][-1][-1].SetName("greOn_{0}_{1}".format(pC, qC))
                    greOn[-1][-1][-1].SetTitle("{0}, PSF{1}%".format(aaStrSelect[pC][qC], cutPsf))
                    greOn[-1][-1][-1].SetMarkerStyle(20+int(cutPsf==68))
                    greZenith[-1][-1].append(ROOT.TGraphErrors())
                    greZenith[-1][-1][-1].SetName("greZenith_{0}_{1}".format(pC, qC))
                    greZenith[-1][-1][-1].SetTitle("{0}, PSF{1}%".format(aaStrSelect[pC][qC], cutPsf))
                    greZenith[-1][-1][-1].SetMarkerStyle(20+int(cutPsf==68))
                    aHtgEvt[-1][-1].append(ROOT.TH3D("htgEvt_{0}_PSF{1}".format(aaStrSelect[pC][qC], cutPsf), "{0} events within PSF{1}{2} from GRB{3} (MET {4} - {5});Cos(Inclination angle);Zenith angle [deg];log_{{10}}Energy [MeV]".format(aaStrSelect[pC][qC], cutPsf, str_fix_psf, nameGrb, metStart, metStop), NBIN_CTH, EDGE_CTH_LOW, EDGE_CTH_UP, NBIN_ZEN, EDGE_ZEN_LOW, EDGE_ZEN_UP, NBIN_ENE, EDGE_ENE_LOW, EDGE_ENE_UP))
                    aHtgInterval[-1][-1].append(ROOT.TH2D("htgInterval_{0}_PSF{1}".format(aaStrSelect[pC][qC], cutPsf), "Intervals (waiting time) between consective {6} {0} events above {7:.1f} GeV within PSF{1}{2} from GRB{3} (MET {4} - {5});MET [s];Waiting time [s]".format(aaStrSelect[pC][qC], cutPsf, str_fix_psf, nameGrb, metStart, metStop, waitingnumber+1, 10**(thresholdenergywaitingtime-3)), NBIN_MET, EDGE_MET_LOW, EDGE_MET_UP, NBIN_ITV, EDGE_ITV_LOW, EDGE_ITV_UP))
                    metPrevious[-1][-1].append([0.]*waitingnumber)
                    if pC==0:
                        greOn[-1][-1][-1].SetMarkerColor(13-12*qC)
                        greZenith[-1][-1][-1].SetMarkerColor(13-12*qC)
                    elif pC==1:
                        greOn[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                        greZenith[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                    mgr.Add(greOn[-1][-1][-1])
                    mgrZenith.Add(greZenith[-1][-1][-1])  

        timeStart = datetime.datetime.now()
        for iEvent in range(nEventChain):
            cbFlagPSF68.value = 0
            cbFlagPSF95.value = 0
            chIn.GetEntry(iEvent)
            ngrb = float(nameGrb)
            grbt = chIn.t - dict_grb["TRIGGER_MET"]
            cdTimeGRB.value = grbt
            vecEvt = np.array([cos(radians(chIn.dec))*cos(radians(chIn.ra)), cos(radians(chIn.dec))*sin(radians(chIn.ra)), sin(radians(chIn.dec))])
            if fixpsfenergy==0:
                epsf = chIn.e
            else:
                epsf = fixpsfenergy
            if fixpsfinclin==0:
                cthpsf = chIn.cth
            else:
                cthpsf = fixpsfinclin
            
            dictDistCut = { 'PSF95': (htgPerf.getPSF95_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384)+3*(chIn.s==32768), epsf, cthpsf) + dict_grb["ERROR_RADIUS"]), 'PSF68': (htgPerf.getPSF68_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384)+3*(chIn.s==32768), epsf, cthpsf) + dict_grb["ERROR_RADIUS"]) }
            radTheta = acos(np.dot(vecTgt, vecEvt))
            degDist = degrees(radTheta)
            cdAngSep.value = degDist
            #if chIn.evid == 6500524:
             #   print "Distance:", degDist, "PSF95:", dictDistCut['PSF95']
            if degDist<dictDistCut['PSF95'] and chIn.t>=metStart and chIn.t<metStop and (chIn.t<exclude[0] or chIn.t>=exclude[1]):
                htgRADEC.Fill(chIn.ra, chIn.dec)
                print ""
                print "== ON photon candidate!!! =="
                cbFlagPSF95.value = 1
                bPSF68=0
                if degDist<dictDistCut['PSF68']:
                    bPSF68=1
                    cbFlagPSF68.value = 1
                if chIn.c == 1:
                    if chIn.s == 4:
                        if chIn.z<90:
                            greOn[bPSF68][0][0].SetPoint(greOn[bPSF68][0][0].GetN(), chIn.t-dict_grb["TRIGGER_MET"], pow(10, chIn.e-3))
                        greZenith[bPSF68][0][0].SetPoint(greZenith[bPSF68][0][0].GetN(), chIn.t-dict_grb["TRIGGER_MET"], chIn.z)
                    elif chIn.s == 128:
                        if chIn.z<90:
                            greOn[bPSF68][0][1].SetPoint(greOn[bPSF68][0][1].GetN(), chIn.t-dict_grb["TRIGGER_MET"], pow(10, chIn.e-3))
                        greZenith[bPSF68][0][1].SetPoint(greZenith[bPSF68][0][1].GetN(), chIn.t-dict_grb["TRIGGER_MET"], chIn.z)
                    if chIn.s >= 4:
                        aHtgEvt[0][0][0].Fill(chIn.cth, chIn.z, chIn.e)
                        if chIn.e>=thresholdenergywaitingtime:
                            #aHtgInterval[0][0][0].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[0][0][0][0])))
                            if metPrevious[0][0][0][0]>0:
                                aHtgInterval[0][0][0].Fill(chIn.t, chIn.t-metPrevious[0][0][0][0])
                            shifttimebox(metPrevious[0][0][0], chIn.t, waitingnumber)
                        #for iwait in range(waitingnumber-1):
                        #    metPrevious[0][0][0][iwait] = metPrevious[0][0][0][iwait+1]
                        #metPrevious[0][0][waitingnumber-1] = chIn.t
                        if bPSF68==1:
                            aHtgEvt[bPSF68][0][0].Fill(chIn.cth, chIn.z, chIn.e)
                            if chIn.e>=thresholdenergywaitingtime:
                                #aHtgInterval[bPSF68][0][0].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[bPSF68][0][0][0])))
                                if metPrevious[bPSF68][0][0][0]>0:
                                    aHtgInterval[bPSF68][0][0].Fill(chIn.t, chIn.t-metPrevious[bPSF68][0][0][0])
                                shifttimebox(metPrevious[bPSF68][0][0], chIn.t, waitingnumber)
                            #metPrevious[bPSF68][0][0] = chIn.t
                    if chIn.s >= 128:
                        aHtgEvt[0][0][1].Fill(chIn.cth, chIn.z, chIn.e)
                        if chIn.e>=thresholdenergywaitingtime:
                            #aHtgInterval[0][0][1].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[0][0][1][0])))
                            if metPrevious[0][0][1][0]>0:
                                aHtgInterval[0][0][1].Fill(chIn.t, chIn.t-metPrevious[0][0][1][0])
                            shifttimebox(metPrevious[0][0][1], chIn.t, waitingnumber)
                        #metPrevious[0][0][1] = chIn.t
                        if bPSF68==1:
                            aHtgEvt[bPSF68][0][1].Fill(chIn.cth, chIn.z, chIn.e)
                            if chIn.e>=thresholdenergywaitingtime:
                                #aHtgInterval[bPSF68][0][1].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[bPSF68][0][1][0])))
                                if metPrevious[bPSF68][0][1][0]>0:
                                    aHtgInterval[bPSF68][0][1].Fill(chIn.t, chIn.t-metPrevious[bPSF68][0][1][0])
                                shifttimebox(metPrevious[bPSF68][0][1], chIn.t, waitingnumber)
                            #metPrevious[bPSF68][0][1] = chIn.t

                elif chIn.c == 2:
                    if chIn.s == 4096:
                        if chIn.z<90:
                            greOn[bPSF68][1][0].SetPoint(greOn[bPSF68][1][0].GetN(), chIn.t-dict_grb["TRIGGER_MET"], pow(10, chIn.e-3))
                        greZenith[bPSF68][1][0].SetPoint(greZenith[bPSF68][1][0].GetN(), chIn.t-dict_grb["TRIGGER_MET"], chIn.z)
                    elif chIn.s == 8192:
                        if chIn.z<90:
                            greOn[bPSF68][1][1].SetPoint(greOn[bPSF68][1][1].GetN(), chIn.t-dict_grb["TRIGGER_MET"], pow(10, chIn.e-3))
                        greZenith[bPSF68][1][1].SetPoint(greZenith[bPSF68][1][1].GetN(), chIn.t-dict_grb["TRIGGER_MET"], chIn.z)
                    elif chIn.s == 16384:
                        if chIn.z<90:
                            greOn[bPSF68][1][2].SetPoint(greOn[bPSF68][1][2].GetN(), chIn.t-dict_grb["TRIGGER_MET"], pow(10, chIn.e-3))
                        greZenith[bPSF68][1][2].SetPoint(greZenith[bPSF68][1][2].GetN(), chIn.t-dict_grb["TRIGGER_MET"], chIn.z)
                    elif chIn.s == 32768:
                        if chIn.z<90:
                            greOn[bPSF68][1][3].SetPoint(greOn[bPSF68][1][3].GetN(), chIn.t-dict_grb["TRIGGER_MET"], pow(10, chIn.e-3))
                        greZenith[bPSF68][1][3].SetPoint(greZenith[bPSF68][1][3].GetN(), chIn.t-dict_grb["TRIGGER_MET"], chIn.z)
                    if chIn.s >= 4096:
                        aHtgEvt[0][1][0].Fill(chIn.cth, chIn.z, chIn.e)
                        if chIn.e>=thresholdenergywaitingtime:
                            #aHtgInterval[0][1][0].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[0][1][0][0])))
                            print chIn.t, chIn.t-metPrevious[0][1][0][0]
                            if metPrevious[0][1][0][0]>0:
                                aHtgInterval[0][1][0].Fill(chIn.t, chIn.t-metPrevious[0][1][0][0])
                            shifttimebox(metPrevious[0][1][0], chIn.t, waitingnumber)
                        #metPrevious[0][1][0] = chIn.t
                        if bPSF68==1:
                            aHtgEvt[bPSF68][1][0].Fill(chIn.cth, chIn.z, chIn.e)
                            if chIn.e>=thresholdenergywaitingtime:
                                #aHtgInterval[bPSF68][1][0].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[bPSF68][1][0][0])))
                                if metPrevious[bPSF68][1][0][0]>0:
                                    aHtgInterval[bPSF68][1][0].Fill(chIn.t, chIn.t-metPrevious[bPSF68][1][0][0])
                                shifttimebox(metPrevious[bPSF68][1][0], chIn.t, waitingnumber)
                            #metPrevious[bPSF68][1][0] = chIn.t
                    if chIn.s >= 8192:
                        aHtgEvt[0][1][1].Fill(chIn.cth, chIn.z, chIn.e)
                        if chIn.e>=thresholdenergywaitingtime:
                            #aHtgInterval[0][1][1].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[0][1][1][0])))
                            if metPrevious[0][1][1][0]>0:
                                aHtgInterval[0][1][1].Fill(chIn.t, chIn.t-metPrevious[0][1][1][0])
                            shifttimebox(metPrevious[0][1][1], chIn.t, waitingnumber)
                        #metPrevious[0][1][1] = chIn.t
                        if bPSF68==1:
                            aHtgEvt[bPSF68][1][1].Fill(chIn.cth, chIn.z, chIn.e)
                            if chIn.e>=thresholdenergywaitingtime:
                                #aHtgInterval[bPSF68][1][1].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[bPSF68][1][1][0])))
                                if metPrevious[bPSF68][1][1][0]>0:
                                    aHtgInterval[bPSF68][1][1].Fill(chIn.t, chIn.t-metPrevious[bPSF68][1][1][0])
                                shifttimebox(metPrevious[bPSF68][1][1], chIn.t, waitingnumber)
                            #metPrevious[bPSF68][1][1] = chIn.t
                    if chIn.s >= 16384:
                        aHtgEvt[0][1][2].Fill(chIn.cth, chIn.z, chIn.e)
                        if chIn.e>=thresholdenergywaitingtime:
                            #aHtgInterval[0][1][2].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[0][1][2][0])))
                            if metPrevious[0][1][2][0]>0:
                                aHtgInterval[0][1][2].Fill(chIn.t, chIn.t-metPrevious[0][1][2][0])
                            shifttimebox(metPrevious[0][1][2], chIn.t, waitingnumber)
                        #metPrevious[0][1][2] = chIn.t
                        if bPSF68==1:
                            aHtgEvt[bPSF68][1][2].Fill(chIn.cth, chIn.z, chIn.e)
                            if chIn.e>=thresholdenergywaitingtime:
                                #aHtgInterval[bPSF68][1][2].Fill(chIn.t, max(0., log10(chIn.t-metPrevious[bPSF68][1][2][0])))
                                if metPrevious[bPSF68][1][2][0]>0:
                                    aHtgInterval[bPSF68][1][2].Fill(chIn.t, chIn.t-metPrevious[bPSF68][1][2][0])
                                shifttimebox(metPrevious[bPSF68][1][2], chIn.t, waitingnumber)
                            #metPrevious[bPSF68][1][2] = chIn.t
                    if chIn.s >= 32768:
                        aHtgEvt[0][1][3].Fill(chIn.cth, chIn.z, chIn.e)
                        if chIn.e>=thresholdenergywaitingtime:
                            if metPrevious[0][1][3][0]>0:
                                aHtgInterval[0][1][3].Fill(chIn.t, chIn.t-metPrevious[0][1][3][0])
                            shifttimebox(metPrevious[0][1][3], chIn.t, waitingnumber)
                        if bPSF68==1:
                            aHtgEvt[bPSF68][1][3].Fill(chIn.cth, chIn.z, chIn.e)
                            if chIn.e>=thresholdenergywaitingtime:
                                if metPrevious[bPSF68][1][3][0]>0:
                                    aHtgInterval[bPSF68][1][3].Fill(chIn.t, chIn.t-metPrevious[bPSF68][1][3][0])
                                shifttimebox(metPrevious[bPSF68][1][3], chIn.t, waitingnumber)

                print "Run ID:", chIn.run
                print "Event ID:", chIn.evid
                print "Event category:", cfg.aStrSelect[chIn.c-1]
                print "Event class:", chIn.s
                print "Time from the trigger:", chIn.t-dict_grb["TRIGGER_MET"], "s"
                print "Anglular distance:", degDist, "deg"
                print "PSF68:", htgPerf.getPSF68_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384)+3*(chIn.s==32768), chIn.e, chIn.cth), "deg"
                print "PSF95:", htgPerf.getPSF95_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384)+3*(chIn.s==32768), chIn.e, chIn.cth), "deg"
                print "Energy:", pow(10,chIn.e-3), "GeV"
                print "Edisp68:", 100*htgPerf.getEdisp68_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384)+3*(chIn.s==32768), chIn.e, chIn.cth), "%"
                print "Cos( inclination angle ):", chIn.cth
                print "Zenith angle:", chIn.z, "deg"
            trGRB.Fill()
            if iEvent%(nEventChain/20)==0:
                rate = int((iEvent*100.)/nEventChain+0.5)
                if rate>0:
                    nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                    meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
                else:
                    meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                sys.stdout.write(meter)
                sys.stdout.flush()

        #trGRB.Merge(chIn.GetName())
        trGRB.AddFriend(chIn)
        cEvent = ROOT.TCanvas("cEvent", "GRB {0} gamma-like events".format(nameGrb))
        cEvent.cd()
        mgr.Draw("AP")
        mgr.GetXaxis().SetTitle("Time [s]")
        mgr.GetYaxis().SetTitle("Energy [GeV]")
        leg = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
        for pD in range(len(aaStrSelect)):
            for qD in range(len(aaStrSelect[pD])):
                leg.AddEntry(greOn[0][pD][qD], greOn[0][pD][qD].GetTitle(), "p")
                leg.AddEntry(greOn[1][pD][qD], greOn[1][pD][qD].GetTitle(), "p")
        leg.Draw("same")
        cZenith = ROOT.TCanvas("cZenith", "Zenith angle of ON/OFF events")
        cZenith.cd()
        mgrZenith.Draw("AP")
        mgrZenith.GetXaxis().SetTitle("Time [s]")
        mgrZenith.GetYaxis().SetTitle("Zenith angle [deg]")
        legZenith = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
        for pD in range(len(aaStrSelect)):
            for qD in range(len(aaStrSelect[pD])):
                legZenith.AddEntry(greZenith[0][pD][qD], greZenith[0][pD][qD].GetTitle(), "p")
                legZenith.AddEntry(greZenith[1][pD][qD], greZenith[1][pD][qD].GetTitle(), "p")

        legZenith.Draw("same")
        print ""
        fileOut.cd()
        htgRADEC.Write()
        cEvent.Write()
        cZenith.Write()
        trGRB.Write()
        for rD in range(2):
            for pD in range(len(aaStrSelect)):
                for qD in range(len(aaStrSelect[pD])):
                    aHtgEvt[rD][pD][qD].Write()
                    aHtgInterval[rD][pD][qD].Write()
        print "Finished!"


if __name__ == '__main__':
    main()
