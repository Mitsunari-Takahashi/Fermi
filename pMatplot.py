#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.odr import Model, Data, RealData, ODR
import pickle_utilities
import logging


TPL_MARKER = ('.', 'D', 's', 'o', 'x', '*', 'p', 'h', '8')
TPL_COLOR = ("r", "g", "b", "c", "m")
TPL_LINE = ('-', '--', ':', '-.')


def find_range_shown(x, y, f):
    nxmax_shown = 0 
    nymax_shown = 0
    nxmin_shown = x.shape[1]-1 #len(x)-1
    nymin_shown = y.shape[0]-1 #len(y)-1
    #print 'Shape of X: {0}'.format(x.shape)
    #print 'Shape of Y: {0}'.format(y.shape)
    for ix, iy in itertools.product(range(x.shape[1]), range(y.shape[0])):
        if f(iy, ix)==True:
            if ix>nxmax_shown:
                nxmax_shown = ix
            if ix<nxmin_shown:
                nxmin_shown = ix
            if iy>nymax_shown:
                nymax_shown = iy
            if iy<nymin_shown:
                nymin_shown = iy
    return ((x[0][nxmin_shown], x[0][nxmax_shown]), (y[nymin_shown][0], y[nymax_shown][0]))



class Data_plotted():
    def __init__(self, label, gr_type, xdata, ydata=None, zdata=None, wdata=None, xdata_err=None, ydata_err=None, ul=False, ll=False):
        self.gr_type = gr_type
        self.xdata = xdata
        self.ydata = ydata
        self.zdata = zdata
        self.wdata = wdata
        self.label = label
        self.xdata_err = xdata_err
        self.ydata_err = ydata_err
        self.ul = ul
        self.ll = ll


class Plot_single():
    """lst_data is a list of lists of class "Data_plotted" objects.
Returns a tuple of Figure and Axes.
"""
    def __init__(self, lst_data, title='', xtitle='', ytitle='', xlog=False, ylog=False, xfigsize=16, yfigsize=10, fontsize=12, zlevels=[1.0]):
        self.datasets = lst_data
        print '{0} plots.'.format(len(self.datasets))
        self.fig = plt.figure(12, 5)
        self.ax = self.fig.add_axes((0.125, 0.125, 0.8, 0.8))
        self.title = title
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.xlog = xlog
        self.ylog = ylog
        if self.xlog==True and self.__class__.__name__!='hist_plot':
            self.ax.set_xscale("log", nonposx='clip')
        if self.ylog==True and self.__class__.__name__!='hist_plot':
            self.ax.set_yscale("log", nonposy='clip')
        self.xfigsize = xfigsize
        self.yfigsize = yfigsize
        self.fontsize = fontsize

        if self.title is not None and len(self.title)>0:
            self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xtitle)
        self.ax.set_ylabel(self.ytitle)
        self.zlevels = zlevels


    def plot(self, xlim=None, ylim=None, path_save=None, figforms=['png', 'pdf']):
        for idset, dset in enumerate(self.datasets):
            print 'No.', idset
            print 'X-axis:', len(dset.xdata)
            print 'Y-axis:', len(dset.ydata)
            if len(self.datasets)<len(TPL_MARKER):
                marker = TPL_MARKER[idset+1]
            if dset.ul==True:
                marker = 'v'
            elif dset.ll==True:
                marker = '^'
            logging.info(dset.label)
            if dset.gr_type=='plot':
                self.ax.plot(dset.xdata, dset.ydata, label=dset.label, fmt=marker)
            elif dset.gr_type=='errorbar':
                self.ax.errorbar(dset.xdata, dset.ydata, xerr=dset.xdata_err, yerr=dset.ydata_err, label=dset.label, fmt=marker)
            elif dset.gr_type=='scatter':
                
                array_size = np.ones_like(dset.zdata)+10*dset.zdata if dset.zdata is not None else 10
                im = plt.scatter(dset.xdata, dset.ydata, s=array_size, c=dset.wdata, label=dset.label, linewidths=2, cmap=cm.nipy_spectral, edgecolors='face')#, vmin=-3.2, vmax=-1.5)#array_color)
                #im = plt.scatter(iris.data[:,0], iris.data[:,1], c=iris.target, linewidths=0, alpha=1, cmap=cm.Accent)
                #plt.colorbar(im.cmap)
               # self.fig.colorbar(im) #, cax=im.cmap, ax=self.ax)#im)
            elif dset.gr_type=='stackplot':
                if len(dset.ydata)==2:
                    self.ax.stackplot(dset.xdata, dset.ydata[0], dset.ydata[1], label=dset.label)
                elif len(dset.ydata)==3:
                    self.ax.stackplot(dset.xdata, dset.ydata[0], dset.ydata[1], dset.ydata[2], label=dset.label)
                elif len(dset.ydata)==4:
                    self.ax.stackplot(dset.xdata, dset.ydata[0], dset.ydata[1], dset.ydata[2], dset.ydata[3], label=dset.label)
            elif dset.gr_type=='contour':
                cont = self.ax.contour(dset.xdata, dset.ydata, dset.zdata, levels=self.zlevels, linewidths=0.5, colors=TPL_COLOR[idset%len(TPL_COLOR)])
                #cont.clabel(fmt='%1.0f'.format(), fontsize=14)
            else:
                logging.critical('{0} has NOT been implmented!!'.format(dset.gr_type))
                sys.exit(1)
        self.ax.legend(loc=0, fancybox=True, framealpha=0.5)
        if xlim is not None:
            self.ax.set_xlim(xlim[0], xlim[1])
        if ylim is not None:
            self.ax.set_ylim(ylim[0], ylim[1])
        for iff in figforms:
            self.fig.savefig('.'.join((path_save, iff)))
        return self.fig, self.ax
                    

class Hist_single(Plot_single):
    """lst_data is a list of 1-D data.
Returns a tuple of Figure and Axes.
"""
    def __init__(self, lst_data, title='', xtitle='', ytitle='', xlog=False, ylog=False, xfigsize=16, yfigsize=10, fontsize=12):
        Plot_single.__init__(self, lst_data, title=title, xtitle=xtitle, ytitle=ytitle, xlog=xlog, ylog=ylog, xfigsize=xfigsize, yfigsize=yfigsize, fontsize=fontsize)


    def plot(self, bins=None, range=None, stacked=True, path_save=None, figforms=['png', 'pdf']):
        if self.xlog==True:
            for idset, dset in enumerate(self.datasets):
                self.datasets[idset].xdata = np.log10(dset.xdata)
        data_series = []
        label_series = []
        for idset, dset in enumerate(self.datasets):
            data_series.append(dset.xdata)
            label_series.append(dset.label)
        #for idset, dset in enumerate(self.datasets):                
        self.ax.hist(data_series, label=label_series, bins=bins, range=range, stacked=stacked, log=self.ylog)
        self.ax.legend(loc=0, fancybox=True, framealpha=0.5)
        if path_save is not None:
            for iff in figforms:
                self.fig.savefig('.'.join((path_save, iff)))
        return self.fig, self.ax


class Curve:
    def __init__(self, quantity, xlabel='time [s]', ylabel='', ul=False, xerr_asym=False, yerr_asym=False):
        self.quantity = quantity
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.lst_xdata = []
        self.lst_ydata = []
        self.xerr_asym = xerr_asym
        self.yerr_asym = yerr_asym
        if self.xerr_asym:
            self.lst_xerr = {'hi':[], 'lo':[]}
        else:
            self.lst_xerr = []
        if self.yerr_asym:
            self.lst_yerr = {'hi':[], 'lo':[]}
        else:
            self.lst_yerr = []
        self.ul = ul
        self.fmt = '.' if self.ul==False else 'v'


    def set_point(self, x, y, xerr=0, yerr=0):
        self.lst_xdata.append(x)
        self.lst_ydata.append(y)
        if True: #self.xerr_asym:
            self.lst_xerr['hi'].append(xerr['hi'])
            self.lst_xerr['lo'].append(xerr['lo'])
        else:
            self.lst_xerr.append(xerr)
        if self.yerr_asym:
            self.lst_yerr['hi'].append(yerr['hi'])
            self.lst_yerr['lo'].append(yerr['lo'])
        else:
            self.lst_yerr.append(yerr)


    def get_n(self):
        n = len(self.lst_xdata)
        if n == len(self.lst_ydata):
            return n
        else:
            logging.error('Numbers of data points are mismatched!!!')
            logging.error('X points: {0}'.format(len(self.lst_xdata)))
            logging.error('Y points: {0}'.format(len(self.lst_ydata)))
            sys.exit(1)


    def get_xdata(self):
        return np.array(self.lst_xdata)


    def get_ydata(self):
        return np.array(self.lst_ydata)


    def get_xerr(self):
        if self.xerr_asym:
            return [np.array(self.lst_xerr['lo']), np.array(self.lst_xerr['hi'])]
        else:
            return np.array(self.lst_xerr)


    def get_maximum(self):
        if self.get_n()>0:
            idx = np.argmax(self.get_ydata())
            return (idx, self.get_xdata()[idx], self.get_ydata()[idx])
        else:
            return 0


    def get_yerr(self):
        logging.info('{0} data points.'.format(self.get_n()))
        if self.get_n()<=2:
            logging.warning('Fitting is skipped!')
        if self.yerr_asym:
            return [np.array(self.lst_yerr['lo']), np.array(self.lst_yerr['hi'])]
        else:
            return np.array(self.lst_yerr)

    def fit(self, t_scale, xmin=-sys.maxint, xmax=sys.maxint):
        # initial guess for the parameters
        params_initial = np.array([1e-4, 1.0]) #a, b
        # function to fit
        def powerlaw(x, a, b):
            return a*(x/t_scale)**(-b)

        def powerlaw_err(x, a, b, aerr, berr, cov):
            val = powerlaw(x, a, b)
            err = np.sqrt(pow(aerr/a, 2)+pow(np.log(x/t_scale), 2)*cov)#*val
            return err

        xdata = []
        ydata = []
        yerr = []
        for ipoint, xpoint in enumerate(self.get_xdata()):
            if xpoint>=xmin and xpoint<xmax:
                xdata.append(xpoint)
                ydata.append(self.get_ydata()[ipoint])
                yerr.append(np.sqrt(self.get_yerr()[0][ipoint]*self.get_yerr()[1][ipoint]))
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        yerr = np.array(yerr)
        if len(xdata)<3:
            logging.warning('Only {0} data points. Fitting is skipped.'.format(len(xdata)))
            return 1
        params_optimal, cov = curve_fit(powerlaw, xdata, ydata, sigma=yerr/ydata, absolute_sigma=False, p0=params_initial)
        logging.info("""Opimized parameters: {0}
Covariance: {1}""".format(params_optimal, cov))
        params_err = np.sqrt(np.diag(cov))
        for iparam, param, in enumerate(params_optimal):
            logging.info("""Parameter {0}: {1} +/- {2}""".format(iparam, params_optimal[iparam], params_err[iparam]))

        x_draw = np.linspace(min(xdata), max(xdata), 100)
        y_draw = np.zeros_like(x_draw)
        yerr_draw = np.zeros_like(y_draw)
        for ix, x in enumerate(x_draw):
            y_draw[ix] = powerlaw(x, params_optimal[0], params_optimal[1])
            yerr_draw[ix] = powerlaw_err(x, params_optimal[0], params_optimal[1], params_err[0], params_err[1], cov[1][1])

        return ((params_optimal, params_err), (x_draw, y_draw, yerr_draw))

    def fit_odr(self, t_scale, xmin=-sys.maxint, xmax=sys.maxint):
        # initial guess for the parameters
        params_initial = [1e-4, 1.0] #np.array([1e-4, 1.0]) #a, b
        # function to fit
        def f(B, x):
            return B[0]*pow(x/t_scale, -B[1])
        powerlaw = Model(f)

        def powerlaw_err(x, a, b, aerr, berr, cov):
            val = f([a,b], x)#f(x, a, b)
            err = np.sqrt(pow(aerr/a, 2)+pow(np.log(x/t_scale), 2)*cov)#*val
            return err

        xdata = []
        ydata = []
        yerr = []
        for ipoint, xpoint in enumerate(self.get_xdata()):
            if xpoint>=xmin and xpoint<xmax:
                xdata.append(xpoint)
                ydata.append(self.get_ydata()[ipoint])
                yerr.append(np.sqrt(self.get_yerr()[0][ipoint]*self.get_yerr()[1][ipoint]))
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        yerr = np.array(yerr)
        if len(xdata)<3:
            logging.warning('Only {0} data points. Fitting is skipped.'.format(len(xdata)))
            return 1
        mydata = RealData(x=xdata, y=ydata, sy=yerr)
        myodr = ODR(mydata, powerlaw, beta0=params_initial)
        myoutput = myodr.run()
        myoutput.pprint()
        params_optimal = myoutput.beta
        cov = myoutput.cov_beta
        params_err = myoutput.sd_beta

        #logging.info("""Opimized parameters: {0}
#Error: {1}""".format(params_optimal, params_err))
        #params_err = np.sqrt(np.diag(cov))
        for iparam, param, in enumerate(params_optimal):
            logging.info("""Parameter {0}: {1} +/- {2}""".format(iparam, params_optimal[iparam], params_err[iparam]))

        x_draw = np.linspace(min(xdata), max(xdata), 100)
        y_draw = np.zeros_like(x_draw)
        yerr_draw = np.zeros_like(y_draw)
        for ix, x in enumerate(x_draw):
            y_draw[ix] = f(params_optimal, x) #f(x, params_optimal[0], params_optimal[1])
            yerr_draw[ix] = powerlaw_err(x, params_optimal[0], params_optimal[1], params_err[0], params_err[1], cov[1][1])

        return ((params_optimal, params_err), (x_draw, y_draw, yerr_draw))


    def fit_lin(self, t_scale, xmin=-sys.maxint, xmax=sys.maxint):
        xdata = []
        ydata = []
        yerr = []
        for ipoint, xpoint in enumerate(self.get_xdata()):
            if xpoint>=xmin and xpoint<xmax:
                xdata.append(xpoint)
                ydata.append(self.get_ydata()[ipoint])
                yerr.append(np.sqrt(self.get_yerr()[0][ipoint]*self.get_yerr()[1][ipoint]))
        xdata = np.array(xdata)
        ydata = np.array(ydata)
        yerr = np.array(yerr)
        if len(xdata)<3:
            logging.warning('Only {0} data points. Fitting is skipped.'.format(len(xdata)))
            return 1
        logx = np.log10(xdata)
        logy = np.log10(ydata)
        logyerr = yerr / ydata

        # define our (line) fitting function
        fitfunc = lambda p, x: p[0] + p[1] * (x-np.log10(t_scale))
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

        pinit = [1.0e-4, -1.0]
        out = optimize.leastsq(errfunc, pinit,
                               args=(logx, logy, logyerr), full_output=1)

        pfinal = out[0]
        covar = out[1]
        print pfinal
        print covar

        index = pfinal[1]
        amp = 10.0**pfinal[0]
        params_optimal = [amp, index]

        indexErr = np.sqrt( covar[1][1] )
        ampErr = np.sqrt( covar[0][0] ) * amp
        params_err = [ampErr, indexErr]

        for iparam, param, in enumerate(params_optimal):
            logging.info("""Parameter {0}: {1} +/- {2}""".format(iparam, params_optimal[iparam], params_err[iparam]))

        logx_draw = np.linspace(np.log10(min(xdata)), np.log10(max(xdata)), 100)
        logy_draw = np.zeros_like(logx_draw)
        logyerr_draw = np.zeros_like(logy_draw)
        for ix, logx in enumerate(logx_draw):
            logy_draw[ix] = fitfunc(pfinal, logx) #fitfunc(params_optimal, logx)
            #logyerr_draw[ix] = powerlaw_err(x, params_optimal[0], params_optimal[1], params_err[0], params_err[1], cov[1][1])

        return ((params_optimal, params_err), (10**logx_draw, 10**logy_draw, logyerr_draw))
