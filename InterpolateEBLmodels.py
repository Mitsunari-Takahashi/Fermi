#!/usr/bin/env python

import sys
#import os
import numpy as np
import csv
from scipy import interpolate, ndimage
#import pandas as pd
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import ROOT


##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


DICT_PATH_MODEL = {'Franceschini2008':'/Users/mitsunari/FermiAnalysis/EBL/Model_Franceschini2008.dat',
                   'Franceschini2016':'/Users/mitsunari/FermiAnalysis/EBL/Model_Franceschini2016.dat',
                   #'Gilmore2012': '/Users/mitsunari/FermiAnalysis/EBL/Gilmore2012_opdep_fiducial.dat'}
                   'Dominguez2011': '/Users/mitsunari/FermiAnalysis/EBL/tau_dominguez11.dat'}


def read_model(pathin):
    redshifts = []
    taus = []
    with open(pathin, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            redshifts.append(float(row[0]))
            taus.append([float(tau) for tau in row[1:]])
                            
    #df = pd.read_csv(pathin, sep=' ', comment='#')
    energies = np.array([float(estr) for estr in header[1:]]) #df.iloc[:,0])
    redshifts = np.array(redshifts) #[float(zstr) for zstr in df.columns[1:]])
    tau = np.array(taus)
    # zmesh, emesh = np.meshgrid(redshifts, energies)
    # tau = df.iloc[:, 1:]
    # tau_interp = interpolate.interp2d(zmesh, emesh, tau, kind='cubic', bounds_error=False, fill_value=np.nan)

    gr2d = ROOT.TGraph2D()
    for iz,z in enumerate(redshifts):
        for ie,e in enumerate(redshifts):
            gr2d.SetPoint(z,e,tau[iz, ie])
    return gr2d #tau_interp


def plot(ax, x, y, z, title):
    ax.set_title(title)
    ax.set_xlabel('Redshift')
    ax.set_ylabel('Energy [GeV]')
    ax.set_yscale('log')
    ax.grid()
    ax.contour(x, y, z, levels=[1,2,3])

    
@click.command()
#@click.argument('name', type=str)
#@click.argument('dst', nargs=-1)
@click.option('--energy', '-e', type=float, default=1.)
@click.option('--redshift', '-z', type=float, default=1.)
#@click.option('--values', type=(str, int))
#@click.option('--values', multiple=True)
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
#@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main():
    energy=1.
    redshift=1.
    loglevel='INFO'
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    fig, ax = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(10, 10))
    logeplot = np.linspace(1.5, 3., num=151) #GeV
    eplot = 10**logeplot
    zplot = np.linspace(0.01, 2., num=200)
    zplot_mesh, eplot_mesh = np.meshgrid(zplot, eplot)
    
    for iax, (model, pathfile) in enumerate(DICT_PATH_MODEL.items()):
        logger.info("===== {0} =====".format(model))
        tau_interp = read_model(pathfile)
        tau_value = tau_interp.Interpolate(redshift, energy*1E6 if model=='Gilmore2012' else energy) #(redshift, energy*1E6 if model=='Gilmore2012' else energy)
        logger.info("""tau = {0}
""".format(tau_value))
        nxax = iax%2
        nyax = iax/2
        tau_plot = np.zeros_like(zplot_mesh) #tau_interp(zplot, eplot*1E3 if model=='Gilmore2012' else eplot*1E-3)
        for iz,z in enumerate(zplot):
            for ie,e in enumerate(eplot*1E3 if model=='Gilmore2012' else eplot*1E-3):
                tau_plot[iz][ie] = gr2d.Interpolate(z,e)
        print tau_plot
        plot(ax[nxax][nyax], zplot_mesh, eplot_mesh, tau_plot, model)
        
    for ff in ['png', 'pdf']:
        fig.savefig('./tau_EBLmodels.{form}'.format(form=ff))

        
if __name__ == '__main__':
    main()
