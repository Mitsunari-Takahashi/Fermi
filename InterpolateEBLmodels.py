#!/usr/bin/env python

import sys
#import os
#import numpy
#import yaml
#import datetime
#from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
ROOT.gROOT.SetBatch()
from pColor import *
#from pMETandMJD import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL

##### Logger #####
logger = getLogger(name)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


def make_graph2D():
    
    

@click.command()
@click.argument('name', type=str)
@click.argument('dst', nargs=-1)
@click.option('--suffix', type=str, default='')
@click.option('--values', type=(str, int))
@click.option('--values', multiple=True)
@click.option('--language', type=click.Choice(['Japanese', 'English']))
@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(name, sed, ra, dec, king, acceptance, livetime, suffix, nside):
        ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    
if __name__ == '__main__':
    main()
