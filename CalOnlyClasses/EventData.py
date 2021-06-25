#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi, sqrt
import matplotlib as mpl
import matplotlib.pyplot as plt
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TFile, TDirectory, TTree, TChain, TH1F, TH2F, TH3F
ROOT.gROOT.SetBatch()
from pColor import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
from ctypes import *


##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


REGIONTYPES = ('healpixels', 'annuli')


class EventData:
    def __init__(self, name, data_path=None, regiontype='annuli', offset_time=0, intervals={'EVENT_CLASS':[{'>=':128}],'TIME':[],'ENERGY':[{'>=':10**4.35, '<':10**5.75}],'ZENITH_ANGLE':[{'<':90.}]}, event_tree_name='EVENTS'):
        self.name = name
        self.offset_time = offset_time
        if data_path!=None:
            self.data_path = data_path
        else:
            self.data_path = './{0}.root'.format(name)
        self.data_file = TFile(self.data_path, "UPDATE")
        if regiontype in REGIONTYPES:
            self.regiontype = regiontype
        else:
            logger.error('regiontype must be one of {0}'.format(REGIONTYPES))
            return 1
        self.eventchain = TChain(event_tree_name)

        self.intervals = intervals
        self.spatical_cuts = []
        

    def add_events(self, event_path):
        self.eventchain.Add(event_path)
        logger.info('Total entries: {0} events'.format(self.eventchain.GetEntries()))

        
    def add_angular_distance_friend(self, origin_name, origin_coordinates):
        ang_dist = c_float()
        self.data_file.cd()
        ang_dist_branch = self.friends[origin_name].Branch('ANG_DIST', ang_dist, 'ANG_DIST/F')

        origin_vec = np.array([cos(radians(origin_coordinates['dec']))*cos(radians(origin_coordinates['ra'])),
                               cos(radians(origin_coordinates['dec']))*sin(radians(origin_coordinates['ra'])),
                               sin(radians(origin_coordinates['ra']))])
        for evt in self.eventchain:
            evt_vec = np.arrray([cos(radians(evt.DEC))*cos(radians(evt.RA)),
                                 cos(radians(evt.DEC))*sin(radians(evt.RA)),
                                 sin(radians(evt.RA))])
            ang_dist.value = degrees(acos(np.dot(origin_vec, evt_vec)))
            ang_dist_branch.Fill()
        self.friends[origin_name].Write("", ROOT.kOverwrite)

            
    def get_cut_string(self, branches=None):
        if branches==None:
            branches = set(self.intervals.kyes())
        branch_cut_list = []
        for branch in branches:
            if branch in self.intervals:
                if len(self.intervals[branch])>0:
                    condition_cut_list = []
                    for condition in self.intervals[branch]:
                        condition_cut_list.append(' && '.join(['{branch}{toffset}{ineq}{vlim}'.format(branch=branch, toffset='-{}'.format(self.offset_time) if branch=='TIME' else '', ineq=ineq, vlim=vlim)
                                                               for ineq, vlim in condition.items()]))
                    branch_cut_list.append('({})'.format(' || '.join(condition_cut_list)))
        return ' && '.join(branch_cut_list)


    def cut_events(self, brnaches=None):
        self.data_file.cd()
        cut_string = self.get_cut_string(branches=branches)
        cuttree = self.eventchain.CopyTree(cut_string)
        cuttree.SetTitle('EVENTS ({})'.format(cut_string))
        return cuttree

        
class EventDataAnnuli(EventData):
    def __init__(self, center_coordinates, rinterval=[{'<':'PSF95'}], data_path=None, offset_time=0, intervals={'EVENT_CLASS':[{'>=':128}],'TIME':[],'ENERGY':[{'>=':10**4.35, '<':10**5.75}],'ZENITH_ANGLE':[{'<':90.}]}, event_tree_name='EVENTS'):
        EventData.__init__(self, name=name, data_path=data_path, regiontype='annuli', offset_time=offset_time, intervals=intervals, event_tree_name=event_tree_name)
        self.center_coordinates = center_coordinates
        
        # if 'ra' in center_coodinates.keys() and 'dec' in center_coodinates.keys():
        #     self.angular_distance_string = 'TMath::Cos(DEC*TMath::DegToRad())*TMath::Cos(RA*TMath::DegToRad()) * TMath::Cos({dec}*TMath::DegToRad())*TMath::Cos({ra}*TMath::DegToRad()) + TMath::Cos(DEC*TMath::DegToRad())*TMath::Sin(RA*TMath::DegToRad()) * TMath::Cos({dec}*TMath::DegToRad())*TMath::Sin({ra}*TMath::DegToRad()) + TMath::Sin(RA*TMath::DegToRad())*TMath::Sin({ra}*TMath::DegToRad())'.format(ra=center_coodinates['ra'], dec=center_coodinates['dec'])
        # elif 'l' in center_coodinates.keys() and 'b' in center_coodinates.keys():
        #     self.angular_distance_string = 'TMath::Cos(B*TMath::DegToRad())*TMath::Cos(L*TMath::DegToRad()) * TMath::Cos({b}*TMath::DegToRad())*TMath::Cos({l}*TMath::DegToRad()) + TMath::Cos(B*TMath::DegToRad())*TMath::Sin(L*TMath::DegToRad()) * TMath::Cos({b}*TMath::DegToRad())*TMath::Sin({l}*TMath::DegToRad()) + TMath::Sin(L*TMath::DegToRad())*TMath::Sin({l}*TMath::DegToRad())'.format(l=center_coodinates['l'], b=center_coodinates['b'])
        # else:
        #     logger.error('The central position should be given by (ra, dec) or (l, b)!!')
        self.friends['CENTER'] = TTree('CENTER')
        self.eventchain.AddFriend(self.friends['CENTER'])
        self.add_angular_distance_friend(self, origin_name='CENTER', origin_coordinates=center_coordinates)
        self.eventchain.Write("", ROOT.kOverwrite)
        self.intervals['CENTER.ANG_DIST'] = rinterval
        
