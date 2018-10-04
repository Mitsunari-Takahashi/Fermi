#!/usr/bin/env python

import sys
#import os
import numpy as np
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from astropy.table import Table, Column
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
from pColor import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import pickle_utilities
import InterpolateEBLmodels
from pMatplot import TPL_LINE

##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)

##### d-loglike value for certain signicicance cuts #####
# Significance(sigma): [ Doubled d-loglike for NDF=1,2,3 ]
TABLE_DLOGLIKE_SIGNIFICANCE = Table({'1.0': [1.00, 2.30, 3.53],
                                     '2.0': [4.00, 6.18, 8.03],
                                     '3.0': [9.00, 11.83, 14.16]
                                     })

##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6


def convert_axis(orig_points, orig_range, plot_range): #X:(136, 804), Y:(485, 26)
    log_orig_points = np.log10(orig_points)
    log_orig_range = np.log10(orig_range)
    return (log_orig_points-log_orig_range[0])/(log_orig_range[1]-log_orig_range[0])*(plot_range[1]-plot_range[0])+plot_range[0]


@click.command()
@click.argument('grb', type=str)
@click.argument('latcurve', type=str)
@click.option('--redshift', '-z', type=float, default=1.)
@click.option('--index', '-i', type=float, default=None)
@click.option('--tmin', type=float, default=100.)
@click.option('--ctaplot', '-c', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/CompareCTA/FIG_CTA_DiffSensitivities_for_ObsTimes.png')
@click.option('--outdir', '-d', type=str, default='.')
@click.option('--outname', '-n', type=str, default='FIG_CTA_Sensitivity')
@click.option('--suffix', type=str, default='')
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(grb, latcurve, redshift, index, tmin, ctaplot, outdir, outname, suffix, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    fig = plt.figure(figsize=(8.91, 5.56))
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))

    if index is not None:
        str_index = '_Index{0:0>2.0f}'.format(index*10)
        index = -index
    else:
        str_index = ''
    #EBL model
    grEBL = InterpolateEBLmodels.read_model(InterpolateEBLmodels.DICT_PATH_MODEL['Franceschini08'])[0]
    
    xran = np.array([1E-2, 200])
    xdata = 10**np.linspace(np.log10(xran[0]), np.log10(xran[1]), 101)
    yran = np.array([3E-14, 9E-10])

    latresults = pickle_utilities.load(latcurve)
    eref_lat = latresults['config']['energy']['ref']
    periods = latresults['results']
    line_counts = {'b':0, 'g':0, 'r':0, 'k':0}
    for period in periods:
        t0 = max(1, period['time']['min'])
        t1 = period['time']['max']
        tref = 10**((np.log10(t0)+np.log10(t1))/2.0) 
        tlabel = '{tmin:.1f} - {tmax:.1f} s'.format(tmin=t0, tmax=t1)
        logger.info('----- {0} -----'.format(tlabel))
        if period['dloglike']['dloglike'][0][0]*2>TABLE_DLOGLIKE_SIGNIFICANCE['2.0'][1]:
            if t0>tmin:
                if t1<1800: #t0<100 and t1>100:
                    color = 'b'
                    lstyle = TPL_LINE[line_counts['b']]
                    line_counts['b'] +=1
                elif t1<18000: #t0<1800 and t1>1800:
                    color = 'g'
                    lstyle = TPL_LINE[line_counts['g']]
                    line_counts['g'] +=1
                elif t1<180000: #t0<18000 and t1>18000:
                    color = 'r'
                    lstyle = TPL_LINE[line_counts['r']]
                    line_counts['r'] +=1
                else:
                    color = 'k'
                    lstyle = TPL_LINE[line_counts['k']]
                    line_counts['k'] +=1
            else:
                color = None
        else:
            color = None
        if color is not None:                
            e2dnde = period['dloglike']['best']['e2dnde'] * MEVtoERG
            logger.debug('e2dnde: {0} erg/cm^2/s'.format(e2dnde))
            if index is None:
                index = period['dloglike']['best']['Index']
            logger.debug('Index: {0}'.format(index))
            ydata = np.array([e2dnde * pow(ene*1e6/eref_lat, index+2) for ene in xdata ])
            tau = np.array([grEBL.Interpolate(redshift, ebin) for ebin in xdata])
            logger.debug("""Tau:
{0}""".format(tau))
            absorp = np.exp(-tau)
            ydata_absorbed = ydata * absorp
            logger.debug("""X values:
{0}""".format(xdata))
            logger.debug("""Y values:
{0}""".format(ydata_absorbed))
            xdata_masked = np.ma.masked_where((xdata < xran[0]) + (xdata >= 3), xdata)
            ydata_masked = np.ma.masked_where((ydata_absorbed < yran[0]) + (ydata_absorbed > yran[1]), ydata_absorbed)
            xplot = convert_axis(xdata_masked, xran, np.array((135, 804)))
            yplot = convert_axis(ydata_masked, yran, np.array((484, 26)))
            ax.plot(xplot, yplot, c=color, label=tlabel, ls=lstyle)  
        else:
            logger.info('Skipping...')

    im_cta = Image.open(ctaplot)
    array_cta = np.asarray(im_cta)
    #logger.debug("""Image array:
#{0}""".format(array_cta.shape))
    plt.imshow(array_cta)
    ax.legend(loc='upper right', fontsize='large', framealpha=1)

    for ff in ('png',):
        fig.savefig('{0}/{1}_GRB{2}{3}{4}.{5}'.format(outdir, outname, grb, str_index, suffix if suffix=='' else '_'+suffix, ff))


if __name__ == '__main__':
    main()
