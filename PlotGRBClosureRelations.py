#!/usr/bin/env python
"""Module for showing closure relations of GRBs.
The authour: Mitsunari Takahashi
"""
import sys
import os
import os.path
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import logging
#import pickle
#import datetime
import numpy as np
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
#from astropy.io import fits
import click
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from matplotlib.ticker import FormatStrFormatter
#import pickle_utilities
import pMatplot
#import pMETandMJD

# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams["font.size"] = 12
# pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
# mpl.rcParams.update(pgf_with_rc_fonts)
NMARKER_STYLE = 10

##### VERSION OF THIS MACRO #####
VERSION = 0.1


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6

##### Refered information #####
GRB_CATALOGUE = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'


DICT_LABEL = {'Fast':{'1st HE': '1',
                  '2nd HE': '2',
                  '1st HE (IC-dom)': '1*'},
              'Slow':{'1st HE': '1',
                      '2nd HE': '2',
                      '1st HE (IC-dom)': '1*'},
              'Radiative':{'1st HE': '1R',
                           '2nd HE': '2R'}
              }
DICT_ALPHA = {'Synchrotron':{'ISM':{'Fast':{'1st HE': lambda p:(3.*p-2.)/4., # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p:1./4.}, # Panaitescu & Kumar, 2000
                                    'Slow':{'1st HE': lambda p:(3.*p-2.)/4., # Panaitescu & Kumar, 2000
                                            '1st HE (IC-dom)': lambda p:3.*p/4.-1./(4.-p), # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p:3.*(p-1.)/4.}, # Panaitescu & Kumar, 2000
                                    'Radiative':{'1st HE': lambda p:2.*(3.*p-1.)/7., #Panaitescu et al., 2006
                                                 '2nd HE': lambda p:4./7.} #Panaitescu et al., 2006
                                    },
                             'Wind':{'Fast':{'1st HE': lambda p:(3.*p-2.)/4., # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p:1./4.}, # Panaitescu & Kumar, 2000
                                    'Slow':{'1st HE': lambda p:(3.*p-2.)/4.,  # Panaitescu & Kumar, 2000
                                            '1st HE (IC-dom)': lambda p:3.*p/4.-p/2./(4.-p), # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p: (3.*p-1.)/4.}, # Panaitescu & Kumar, 2000
                                    'Radiative':{'1st HE': lambda p:(5.*p-2.)/6., #Panaitescu et al., 2006
                                                 '2nd HE': lambda p:1./2.} #Panaitescu et al., 2006
                                 }
                         },
              'SSC':{'ISM':{'Fast':{'1st HE': lambda p:(9.*p-10.)/8., # Panaitescu & Kumar, 2000
                                    '2nd HE': lambda p:-1./8.}, # Panaitescu & Kumar, 2000
                            'Slow':{'1st HE': lambda p:(9.*p-10.)/8., # Panaitescu & Kumar, 2000
                                    '1st HE (IC-dom)': lambda p:(9.*p-10.)/8.-(p-2.)/(4.-p), #Sari & Esin (2001)
                                    '2nd HE': lambda p:(9.*p-11.)/8.} # Panaitescu & Kumar, 2000. # IC-dom is same as this according to Sari & Esin (2001)
                        },
                     'Wind':{'Fast':{'1st HE': lambda p:(p-1.), # Panaitescu & Kumar, 2000
                                     '2nd HE': lambda p:0}, # Panaitescu & Kumar, 2000
                             'Slow':{'1st HE': lambda p:(p-1.), # Panaitescu & Kumar, 2000
                                     '1st HE (IC-dom)': lambda p:p*(3.-p)/(4.-p),
                                     '2nd HE': lambda p: p} # Panaitescu & Kumar, 2000
                         } # Panaitescu & Kumar, 2000
                 }
              }
                            
    
DICT_BETA = {'Synchrotron':{'ISM':{'Fast':{'1st HE': lambda p:p/2.,  # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p:1./2.}, # Panaitescu & Kumar, 2000
                                    'Slow':{'1st HE': lambda p:p/2., # Panaitescu & Kumar, 2000
                                            '1st HE (IC-dom)': lambda p:p/2., # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p:(p-1.)/2.}, # Panaitescu & Kumar, 2000
                                    'Radiative':{'1st HE': lambda p:p/2.,
                                                 '2nd HE': lambda p:1./2.}
                                   },
                            'Wind':{'Fast':{'1st HE': lambda p:p/2., # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p:1./2.}, # Panaitescu & Kumar, 2000
                                    'Slow':{'1st HE': lambda p:p/2.,  # Panaitescu & Kumar, 2000
                                            '1st HE (IC-dom)': lambda p:p/2., # Panaitescu & Kumar, 2000
                                            '2nd HE': lambda p: (p-1.)/2.}, # Panaitescu & Kumar, 2000
                                    'Radiative':{'1st HE': lambda p:p/2.,
                                                 '2nd HE': lambda p:1./2.}
                                }
                        },
             'SSC':{'ISM':{'Fast':{'1st HE': lambda p:p/2., #Sari & Esin (2001)
                                   '2nd HE': lambda p:1./2.}, #Sari & Esin (2001)
                           'Slow':{'1st HE': lambda p:p/2., #Sari & Esin (2001)
                                   '1st HE (IC-dom)': lambda p:p/2., #Sari & Esin (2001) # Same as Synchrotron-dominant
                                   '2nd HE': lambda p:(p-1.)/2.} #Sari & Esin (2001)
                       },
                    'Wind':{'Fast':{'1st HE': lambda p:p/2., #Sari & Esin (2001)
                                    '2nd HE': lambda p:1./2.}, #Sari & Esin (2001)
                            'Slow':{'1st HE': lambda p:p/2., #Sari & Esin (2001)
                                    '1st HE (IC-dom)': lambda p:p/2.,
                                    '2nd HE': lambda p: (p-1.)/2.} #Sari & Esin (2001)
                        }
                }
         }


class ObservedIndices:
    def __init__(self, alpha=None, beta=None, alpha_err=0., beta_err=0., name=None):
        self.name = name
        self.alpha = alpha
        self.alpha_err = alpha_err
        self.beta = beta
        self.beta_err = beta_err


    def read(self, indata, instruments=None): 
        ndata = sum(1 for line in open(indata))-1
        with open(indata, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            dict_ncol = {}
            for irow, row in enumerate(header):
                dict_ncol[row] = irow
            alpha = []
            alpha_err_hi = []
            alpha_err_lo = []
            beta = []
            beta_err_hi = []
            beta_err_lo = []
            name_data = None
            for irow, row in enumerate(reader):
                if instruments is None or row[dict_ncol['Instrument']] in instruments:
                    alpha.append(float(row[dict_ncol['alpha']]))
                    alpha_err_hi.append(float(row[dict_ncol['alpha_err_hi']]))
                    alpha_err_lo.append(float(row[dict_ncol['alpha_err_lo']]))
                    print '{v} + {eh} -{el}'.format(v=alpha[irow], eh=alpha_err_hi[irow], el=alpha_err_lo[irow])
                    beta.append(float(row[dict_ncol['beta']]))
                    beta_err_hi.append(float(row[dict_ncol['beta_err_hi']]))
                    beta_err_lo.append(float(row[dict_ncol['beta_err_lo']]))
                    print '{v} + {eh} -{el}'.format(v=beta[irow], eh=beta_err_hi[irow], el=beta_err_lo[irow])
                if irow==0:
                    name_data = row[dict_ncol['GRB']]

            alpha = np.array(alpha)
            alpha_err_hi = np.array(alpha_err_hi)
            alpha_err_lo = np.array(alpha_err_lo)
            beta = np.array(beta)
            beta_err_hi = np.array(beta_err_hi)
            beta_err_lo = np.array(beta_err_lo)

        self.name = name_data
        self.alpha = alpha
        self.alpha_err = {'err_hi':alpha_err_hi,'err_lo':alpha_err_lo}
        self.beta = beta
        self.beta_err = {'err_hi':beta_err_hi,'err_lo':beta_err_lo}


    def draw(self, ax):
        print 'alpha:', self.alpha
        print 'alpha_err:', self.alpha_err
        print 'beta:', self.beta
        print 'beta_err:', self.beta_err
        ax.errorbar(x=self.alpha, y=self.beta, xerr=[self.alpha_err['err_lo'],self.alpha_err['err_hi']], yerr=[self.beta_err['err_lo'],self.beta_err['err_hi']], label=self.name, c='k', lw=1, fmt='.')

        

class ClosureRelation:
    def __init__(self, alpha, beta, emission, cbmprof, cooling, esegment, name='', expression=''):
        """Input alpha and beta as lambda functions and name and expression as string
"""
        self.name = name #str
        self.alpha = alpha #lambda
        self.beta = beta #lambda
        self.expression = expression #str
        self.emission = emission
        self.cbmprof = cbmprof
        self.cooling = cooling
        self.esegment = esegment


    def eval_alpha(self, p):
        return self.alpha(p)


    def eval_beta(self, p):
        return self.beta(p)


    def draw(self, ax, p_range=(2, 3.0), npoint=200, **kwargs):
        p_indices = np.linspace(p_range[0], p_range[1], int((p_range[1]-p_range[0])*npoint)+1)
        alpha_indices = np.zeros_like(p_indices)
        beta_indices = np.zeros_like(p_indices)
        for ip,p in enumerate(p_indices):
            alpha_indices[ip] = self.eval_alpha(p)
            beta_indices[ip] = self.eval_beta(p)
        pointlike = True if all([a==alpha_indices[0] for a in alpha_indices]) and all([b==beta_indices[0] for b in beta_indices]) else False
        self.im = ax.scatter(alpha_indices, beta_indices, c=p_indices if not pointlike else 'grey', cmap=cm.rainbow, marker='o', vmin=min(p_indices), vmax=max(p_indices), edgecolors=None, linewidths=0, s=30 if pointlike else 3, **kwargs)
        xtext = alpha_indices[-1] -0.03 
        if xtext < np.mean(alpha_indices):
            xtext-=0.05
        ytext = beta_indices[-1] + 0.03
        if self.name[:3]=='Syn':
            stext = 'S' +DICT_LABEL[self.cooling][self.esegment]
            #stext = 'S' +DICT_LABEL[self.name[4:]]
        elif self.name[:3]=='SSC':
            stext = 'C'+DICT_LABEL[self.cooling][self.esegment]
            #stext = 'C'+DICT_LABEL[self.name[4:]]
        ax.text(x=xtext, y=ytext, s=stext, fontsize=8)


@click.command()
@click.argument('name', type=str)
@click.argument('indata', type=str)
@click.option('--emission', type=click.Choice(['Synchrotron', 'SSC', 'both']))
@click.option('--cbm', type=click.Choice(['ISM', 'Wind', 'both']))
@click.option('--cooling', type=click.Choice(['Fast', 'Slow', 'both']))
@click.option('--suffix', '-s', type=str, default='')
@click.option('--pathout', type=str, default='./ClosureRelations')
@click.option('--figform', type=str, default=('png',), multiple=True)
def main(name, indata, emission, cbm, cooling, suffix, pathout, figform):
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_axes((0.15, 0.15, 0.75, 0.75))
    #im = ax.imshow(data, cmap='rainbow')
    
    suffix = suffix if suffix=='' else '_'+suffix
    flag_cbar = False
    for em in (emission,) if emission!='both' else ('Synchrotron', 'SSC'):
        for cb in (cbm,) if cbm!='both' else ('ISM', 'Wind'):
            for coo in (cooling,) if cooling!='both' else ('Fast', 'Slow'):
                for eseg, formula in DICT_ALPHA[em][cb][coo].items():
                    str_name = """{em} {cb} 
{coo} {eseg}""".format(em=em if emission=='both' else '', cb=cb if cbm=='both' else '', coo=coo if cooling=='both' else '', eseg=eseg)
                    clrel = ClosureRelation(alpha=formula, beta=DICT_BETA[em][cb][coo][eseg], name=str_name)
                    clrel.draw(ax)
                    if flag_cbar==False:
                        cbar = fig.colorbar(clrel.im)
                        cbar.set_label('p of the electrons')
                        flag_cbar = True
    # ndata = sum(1 for line in open(indata))-1
    # with open(indata, 'r') as f:
    #     reader = csv.reader(f)
    #     header = next(reader)
    #     dict_ncol = {}
    #     for irow, row in enumerate(header):
    #         dict_ncol[row] = irow
    #     alpha = np.zeros(ndata)
    #     alpha_err_hi = np.zeros(ndata)
    #     alpha_err_lo = np.zeros(ndata)
    #     beta = np.zeros(ndata)
    #     beta_err_hi = np.zeros(ndata)
    #     beta_err_lo = np.zeros(ndata)
    #     name_data = None
    #     for irow, row in enumerate(reader):
    #         alpha[irow] = float(row[dict_ncol['alpha']])
    #         alpha_err_hi[irow] = float(row[dict_ncol['alpha_err_hi']])
    #         alpha_err_lo[irow] = float(row[dict_ncol['alpha_err_lo']])
    #         print '{v} + {eh} -{el}'.format(v=alpha[irow], eh=alpha_err_hi[irow], el=alpha_err_lo[irow])
    #         beta[irow] = float(row[dict_ncol['beta']])
    #         beta_err_hi[irow] = float(row[dict_ncol['beta_err_hi']])
    #         beta_err_lo[irow] = float(row[dict_ncol['beta_err_lo']])
    #         print '{v} + {eh} -{el}'.format(v=beta[irow], eh=beta_err_hi[irow], el=beta_err_lo[irow])
    #         if irow==0:
    #             name_data = row[dict_ncol['GRB']] #','.join([row[dict_ncol['GRB']],row[dict_ncol['Band']],
                #row[dict_ncol['Time']]])
    #print name_data
    # obs = ObservedIndices(alpha=alpha,
    #                       beta=beta,
    #                       alpha_err={'err_hi':alpha_err_hi,'err_lo':alpha_err_lo},
    #                       beta_err={'err_hi':beta_err_hi,'err_lo':beta_err_lo},
    #                       name=name_data)
    obs = ObservedIndices()
    obs.read(indata)
    obs.draw(ax)
    ax.legend()
    ax.grid()
    str_title = ''
    str_add = ''
    li_title = []
    str_label = ''
    dict_title_suffix = {emission:'', cbm:'-like', cooling:'-cooling'}
    for cha in (emission, cbm, cooling):
        if cha is not 'both':
            li_title.append(cha + dict_title_suffix[cha])
            str_add = '_' + cha + str_add
    str_title = ', '.join(li_title)
    ax.set_title(str_title)
    ax.set_xlabel('alpha')
    ax.set_ylabel('beta')
    ax.set_xlim((0.5, 2.0))
    ax.set_ylim((0.25, 1.75))
    for ff in figform:
        fig.savefig('{0}{1}{2}.{3}'.format(pathout, str_add, suffix, ff))
    

if __name__ == '__main__':
    main()
