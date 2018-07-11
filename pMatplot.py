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
TPL_COLOR = ("k", "r", "g", "b", "c", "hotpink", "brown", "purple", "olive", "orange")
TPL_COLOR_WO_BLACK = ("r", "g", "b", "c", "hotpink", "brown", "purple", "olive", "orange")
TPL_LINE = ('-', '--', '-.', ':')


def find_range_shown(y, x, f, ny_restricted=None, nx_restricted=None):
    yindexrange, xindexrange = find_indexrange(y, x, f, ny_restricted=ny_restricted, nx_restricted=nx_restricted)
    print 'Y index range: {0}'.format(yindexrange)
    print 'X index range: {0}'.format(xindexrange)
    #print 'Y shape: {0}'.format(y.shape)
    #print 'X shape: {0}'.format(x.shape)
    return ((y[yindexrange[0][0], yindexrange[0][1]], y[yindexrange[1][0], yindexrange[1][1]]), 
            (x[xindexrange[0][0], xindexrange[0][1]], x[xindexrange[1][0], xindexrange[1][1]]))


def find_indexrange(y, x, f, ny_restricted=None, nx_restricted=None):
    nxmax_shown = None
    nymax_shown = None
    nxmin_shown = None
    nymin_shown = None
    for iy, ix in itertools.product(range(y.shape[0]), range(x.shape[1])):
        if f(iy, ix)==True and (nx_restricted is None or ix in nx_restricted) and (ny_restricted is None or iy in ny_restricted):
            if nxmax_shown is None or x[iy,ix]>x[nxmax_shown[0],nxmax_shown[1]]:
                nxmax_shown = (iy,ix)
            if nymax_shown is None or y[iy,ix]>y[nymax_shown[0],nymax_shown[1]]:
                nymax_shown = (iy,ix)
            if nxmin_shown is None or x[iy,ix]<x[nxmin_shown[0],nxmin_shown[1]]:
                nxmin_shown = (iy,ix)
            if nymin_shown is None or y[iy,ix]<y[nymin_shown[0],nymin_shown[1]]:
                nymin_shown = (iy,ix)
    return ((nymin_shown, nymax_shown), (nxmin_shown, nxmax_shown))    


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
        self.fig = plt.figure(figsize=(12, 5))
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
    def __init__(self, quantity, xlabel='time [s]', ylabel='', ul=False, ll=False, xerr_asym=False, yerr_asym=False):
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
        self.ll = ll
        if self.ul==True:
            self.fmt =  'v'
        elif self.ul==True:
            self.fmt =  '^'
        else:
            self.fmt = '.' 


    def set_point(self, x, y, xerr=0, yerr=0):
        self.lst_xdata.append(x)
        self.lst_ydata.append(y)
        if self.xerr_asym:
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
            logging.info('Fitting is skipped!')
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
            logging.info('Only {0} data points. Fitting is skipped.'.format(len(xdata)))
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
            logging.info('Only {0} data points. Fitting is skipped.'.format(len(xdata)))
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
            logging.info('Only {0} data points. Fitting is skipped.'.format(len(xdata)))
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

        return ((params_optimal, params_err), (10**logx_draw, 10**logy_draw, logyerr_draw))


class ColorMeshCurve(Curve):
    """"Light curve represented by pcolormesh of matplotlib
https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.pcolormesh.html
"""
    def __init__(self, quantity, xlabel='time [s]', ylabel='', zlabel=''):
        Curve.__init__(self, quantity, xlabel=xlabel, ylabel=ylabel, ul=False, ll=False, xerr_asym=True, yerr_asym=False)
        self.zlabel = zlabel
        self.lst_zdata = []


    def set_point(self, x, y, z, xerr=0, yerr=0):
        Curve.set_point(self, x=x, y=y, xerr=xerr, yerr=yerr)
        self.lst_zdata.append(z)


    def get_zdata(self):
        return np.array(self.lst_zdata)


    def make_meshes(self):
        # X-axis
        xcenters = self.get_xdata()
        # Z-axis
        zvalues = self.get_zdata()
        print 'z-valus: {0}'.format(zvalues)

        if self.xerr_asym:
            xerr_lo, xerr_hi = self.get_xerr()
        else:
            xerr_lo, xerr_hi = self.get_xerr(), self.get_xerr()
        xlowedges, xupedges = xcenters - xerr_lo, xcenters + xerr_hi
        list_xgrid = [xlowedges[0]]
        list_zvalues = []
        list_zmask = []
        for ix, x in enumerate(xcenters):
            #if ix==0:
            #    list_xgrid.append(xlowedges[ix])
            if abs(xlowedges[ix]-list_xgrid[-1])<0.01:
                list_xgrid.append(xupedges[ix])
            elif xlowedges[ix]>list_xgrid[-1]:
                list_zvalues.append(np.full_like(list_zvalues[0], float(sys.maxint)))
                list_zmask.append(np.full_like(zvalues[0], True, dtype=bool))
                list_xgrid.append(xlowedges[ix])
                list_xgrid.append(xupedges[ix])
            else:
                logging.critical('X bin edge (No.{0}) {1} is smaller than the previous one {2}!!!'.format(ix, xlowedges[ix], list_xgrid[-1]))
                logging.critical('Difference: {0}'.format(xlowedges[ix]-list_xgrid[-1]))
                logging.critical('Filled grids: {0}'.format(list_xgrid))
                logging.critical('Lower edges: {0}'.format(xlowedges))
                logging.critical('Upper edges: {0}'.format(xupedges))
                sys.exit(1)
            list_zvalues.append(zvalues[ix])
            list_zmask.append(np.full_like(zvalues[0], False, dtype=bool))
        xedges = np.array(list_xgrid)
        self.z_mesh = np.ma.array(list_zvalues, mask=list_zmask).T
        print 'Z mesh: {0}'.format(self.z_mesh)

        # Y-axis
        for iys, ys in enumerate(self.lst_ydata):
            if any(ys!=self.lst_ydata[0]):
                logging.critical('Y bin edges do NOT match!!!')
                logging.critical('In the first time bin: {0}'.format(self.lst_ydata[0]))
                logging.critical('In the {0}th time bin: {1}'.format(iys, ys))
                sys.exit(1)
        yedges = self.lst_ydata[0]

        # Y-X mesh
        self.x_mesh, self.y_mesh = np.meshgrid(xedges, yedges)

        print 'X mesh: {0}'.format(self.x_mesh.shape)
        print 'Y mesh: {0}'.format(self.y_mesh.shape)
        print 'Z mesh: {0}'.format(self.z_mesh.shape)



            
            
