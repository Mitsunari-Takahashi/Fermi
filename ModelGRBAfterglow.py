#!/usr/bin/env python

import sys
import numpy as np
from math import pi
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import animation
#from PIL import Image, ImageDraw
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
from scipy.constants import c,h,e
from scipy import optimize
from astropy.cosmology import WMAP9 as cosmo
import astropy.constants as constants
from astropy.table import Table, Column, QTable
from astropy.io import ascii
from astropy import units as u
import pickle_utilities
import CalcCharaFreq
import InterpolateEBLmodels
import pMatplot

rc('text', usetex=True)

##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'DEBUG'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)

DICT_COLOR = {'LAT': 'r',
              'XRT': 'b',
              'c':'r',
              'm':'b'}

DICT_LINE = {'Sync.': '--',
            'SSC': ':'}


class SRC_PARAMETERS:
    def __init__(self, p,z,e_erg,epsilon_e,epsilon_B,n_0,a_s,t_sec):
        self.p = p
        self.z = z
        self.e_erg = e_erg
        self.e_52 = self.e_erg * 1E-52
        self.epsilon_e = epsilon_e
        self.epsilon_B = epsilon_B
        self.n_0 = n_0
        self.a_s = a_s
        self.t_sec = t_sec
        self.t_min = self.t_sec / 60.
        self.t_hr = self.t_sec / 3600.
        self.t_day = self.t_sec / 86400.
        self.epsilon_e_bar = self.epsilon_e*(self.p-2.)/(self.p-1.)
        self.d_L = cosmo.luminosity_distance(self.z) #Weinberg, 1972, pp 420-424; Weedman, 1986, pp 60-62.
        self.d_L_Mpc = self.d_L.value
        self.d_L_cm = self.d_L.to(u.cm).value
        self.d_L28 = self.d_L_cm * 1E-28


class ObservedLightCurve:
    def __init__(self, name, title):
        self.name = name
        self.title = title
        self.xdata = None
        self.ydata = None
        self.xerr = None
        self.yerr = None

    def load(self, qtable, colx, coly, colxerr=None, colyerr=None, selector=None):
        if selector is None:
            qtb_filtered = qtable
        else:
            qtb_filtered = qtable[np.where(qtable[selector]==self.name)]

        self.xdata = qtb_filtered[colx].to(u.s).value
        self.ydata = qtb_filtered[coly].to(u.erg/u.cm/u.cm/u.s).value

        if colxerr is None:
            self.xerr = np.zeros_like(xdata)
        elif isinstance(colxerr, str):
            self.xerr = qtb_filtered[colxerr].to(u.s).value
        elif isinstance(colxerr, list) or isinstance(colxerr, tuple):
            self.xerr = np.array([np.abs(qtb_filtered[colxerr[0]].to(u.s).value), np.abs(qtb_filtered[colxerr[1]].to(u.s).value)])
        else:
            logger.critical('colxerr: {0}'.format(colxerr))
            logger.critical('Input None or str or (str, str) as colxerr')
            sys.exit(1)

        if colyerr is None:
            self.yerr = np.zeros_like(ydata)
        elif isinstance(colyerr, str):
            self.yerr = qtb_filtered[colyerr].to(u.erg/u.cm/u.cm/u.s).value
        elif isinstance(colyerr, list) or isinstance(colyerr, tuple):
            self.yerr = np.array([np.abs(qtb_filtered[colyerr[0]].to(u.erg/u.cm/u.cm/u.s).value), np.abs(qtb_filtered[colyerr[1]].to(u.erg/u.cm/u.cm/u.s).value)])
        else:
            logger.critical('colyerr: {0}'.format(colyerr))
            logger.critical('Input None or str or (str, str) as colyerr')
            sys.exit(1)

    def plot(self, ax):
        ax.errorbar(x=self.xdata, y=self.ydata, xerr=self.xerr, yerr=self.yerr, label=self.title, fmt='.', lw=1, c=DICT_COLOR[self.name])


def read_observation_table(pathin):
    if pathin[-4:]=='.csv':
        table = ascii.read(pathin)
        qtable = QTable(table)
        for key in qtable.keys():
            if key[-3:]=='[s]' or key[-5:]=='[sec]':
                qtable[key] = qtable[key] * u.s
            if key[-5:]=='[min]':
                qtable[key] = qtable[key] * u.min
            if key[-3:]=='[h]' or key[-4:]=='[hr]' or key[-6:]=='[hour]':
                qtable[key] = qtable[key] * u.s
            if key[-3:]=='[d]' or key[-5:]=='[day]':
                qtable[key] = qtable[key] * u.d
            if key[-5:]=='[mag]':
                qtable[key] = qtable[key] * u.Magnitude
            if key[-4:]=='[Jy]':
                qtable[key] = qtable[key] * u.Jy
            if key[-5:]=='[mJy]':
                qtable[key] = qtable[key] * u.mJy
            if key[-12:]=='[erg/cm^2/s]':
                qtable[key] = qtable[key] * u.erg / u.cm / u.cm / u.s
    elif pathin[-7:]=='.pickle':
        pdata = pickle_utilities.load(pathin)['results']
        time = []
        time_err_lo = []
        time_err_hi = []
        eflux = []
        eflux_err_lo = []
        eflux_err_hi = []
        for period in pdata:
            time.append(np.sqrt(period['time']['min']*period['time']['max']))
            time_err_lo.append(time[-1] - period['time']['min'])
            time_err_hi.append(period['time']['max'] - time[-1])
            
            dll_map = period['dloglike']['dloglike']
            eflux_map = period['dloglike']['eflux']

            args_bestlike = zip(np.where(dll_map==dll_map.min())[0], np.where(dll_map==dll_map.min())[1])[0]

            eflux_shown = eflux_map[2.*dll_map<=2.30]
            eflux_min = min(eflux_shown.flatten())
            eflux_max = max(eflux_shown.flatten())
            eflux.append(eflux_map[args_bestlike[0]][args_bestlike[1]])
            eflux_err_hi.append(eflux_max-eflux[-1]) 
            eflux_err_lo.append(eflux[-1]-eflux_min) 
        qtable = QTable()
        qtable['time'] = Column(time, unit=u.s)
        qtable['time_err_lo'] = Column(time_err_lo, unit=u.s)
        qtable['time_err_hi'] = Column(time_err_hi, unit=u.s)
        qtable['eflux'] = Column(eflux, unit=u.MeV/u.cm/u.cm/u.s)
        qtable['eflux_err_lo'] = Column(eflux_err_lo, unit=u.MeV/u.cm/u.cm/u.s)
        qtable['eflux_err_hi'] = Column(eflux_err_hi, unit=u.MeV/u.cm/u.cm/u.s)

    return qtable


DICT_FREQ_BREAK = {'Slow':{'m':{'ISM': lambda sp: 3.73*(sp.p-0.67)*1E15*(1.+sp.z)**0.5*sp.e_52**0.5*sp.epsilon_e_bar**2*sp.epsilon_B**0.5*sp.t_day**-1.5, #Granot and Sari (2002) 2
                                'Wind': lambda sp: 4.02*(sp.p-0.69)*1E15*(1.+sp.z)**0.5*sp.e_52**0.5*sp.epsilon_e_bar**2*sp.epsilon_B**0.5*sp.t_day**-1.5}, #2
                           'c':{'ISM': lambda sp: 6.37*(sp.p-0.46)*1E13*np.exp(-1.16*sp.p)*(1.+sp.z)**-0.5*sp.epsilon_B**-1.5*sp.n_0**-1*sp.e_52**-0.5*sp.t_day**-0.5, #3
                                'Wind': lambda sp: 4.40*(3.45-sp.p)*1E10*np.exp(0.45*sp.p)*(1.+sp.z)**-1.5*sp.epsilon_B**-1.5*sp.a_s**-2*sp.e_52**0.5*sp.t_day**0.5}}, #3
                   'Fast':{'m':{'ISM': lambda sp: 3.94*(sp.p-0.74)*1E15*(1.+sp.z)**0.5*sp.epsilon_e_bar**2*sp.epsilon_B**0.5*sp.e_52**0.5*sp.t_day**-1.5,
                                'Wind': lambda sp: 3.52*(sp.p-0.31)*1E15*(1.+sp.z)**0.5*sp.epsilon_e_bar**2*sp.epsilon_B**0.5*sp.e_52**0.5*sp.t_day**-1.5},
                           'c':{'ISM': lambda sp: 5.86*1E12*(1.+sp.z)**-0.5*sp.epsilon_B**-1.5*sp.n_0**-1*sp.e_52**-0.5*sp.t_day**-0.5,
                                'Wind': lambda sp: np.nan}}}


DICT_FLUX = {'D': {'ISM': lambda nu_14,sp: 27.9*(sp.p-1.)/(3*sp.p-1.)*pow(1.+sp.z, 5./6.)*pow(sp.epsilon_e_bar,-2./3.)*pow(sp.epsilon_B,1./3.)*pow(sp.n_0, 1./2.)*pow(sp.e_52,5./6.)*pow(sp.t_day,1./2.)*pow(sp.d_L28, -2)*pow(nu_14,1./3.),
                   'Wind': lambda nu_14,sp: 211.*(sp.p-1.)/(3*sp.p-1.)*pow(1.+sp.z, 4./3.)*pow(sp.epsilon_e_bar,-2./3.)*pow(sp.epsilon_B,1./3.)*sp.a_s*pow(sp.e_52,1./3.)*pow(sp.d_L28, -2)*pow(nu_14,1./3.)},
             'E': {'ISM': lambda nu_14,sp: 73.*pow(1.+sp.z, 7./6.)*sp.epsilon_B*pow(sp.n_0, 5./6.)*pow(sp.e_52, 7./6.)*pow(sp.t_day,1./6.)*pow(sp.d_L28,-2)*pow(nu_14, 1./3.),
                   'Wind': np.nan},
             'F': {'ISM': lambda nu_14,sp: 6.87*pow(1.+sp.z, 3./4.)*pow(sp.epsilon_B,-1./4.)*pow(sp.e_52, 3./4.)*pow(sp.t_day,-1./4.)*pow(sp.d_L28,-2)*pow(nu_14, -1./2.),
                   'Wind': lambda nu_14,sp: 6.68*pow(1.+sp.z, 3./4.)*pow(sp.epsilon_B,-1./4.)*pow(sp.e_52, 3./4.)*pow(sp.t_day,-1./4.)*pow(sp.d_L28,-2)*pow(nu_14, -1./2.)},
             'G': {'ISM': lambda nu_14,sp: 0.461*(sp.p-0.04)*np.exp(2.53*sp.p)*pow(1.+sp.z, (3.+sp.p)/4.)*pow(sp.epsilon_e_bar,sp.p-1.)*pow(sp.epsilon_B,(1.+sp.p)/4.)*pow(sp.n_0,1./2.)*pow(sp.e_52,(3.+sp.p)/4.)*pow(sp.t_day,3.*(1.-sp.p)/4.)*pow(sp.d_L28,-2)*pow(nu_14,(1.-sp.p)/2.),
                   'Wind': lambda nu_14,sp: 3.82*(sp.p-0.18)*np.exp(2.54*sp.p)*pow(1.+sp.z, (5.+sp.p)/4.)*pow(sp.epsilon_e_bar,sp.p-1.)*pow(sp.epsilon_B,(1.+sp.p)/4.)*sp.a_s*pow(sp.e_52,(1.+sp.p)/4.)*pow(sp.t_day,(1.-3.*sp.p)/4.)*pow(sp.d_L28,-2)*pow(nu_14,(1.-sp.p)/2.)},
             'H': {'ISM': lambda nu_14,sp: 0.855*(sp.p-0.98)*np.exp(1.95*sp.p)*pow(1.+sp.z, (2.+sp.p)/4.)*pow(sp.epsilon_e_bar,sp.p-1.)*pow(sp.epsilon_B,(sp.p-2.)/4.)*pow(sp.e_52,(2.+sp.p)/4.)*pow(sp.t_day,(2.-3.*sp.p)/4.)*pow(sp.d_L28,-2)*pow(nu_14,-sp.p/2.),
                   'Wind': lambda nu_14,sp: 0.0381*(7.11-sp.p)*np.exp(2.76*sp.p)*pow(1.+sp.z, (2.+sp.p)/4.)*pow(sp.epsilon_e_bar,sp.p-1.)*pow(sp.epsilon_B,(sp.p-2.)/4.)*pow(sp.e_52,(2.+sp.p)/4.)*pow(sp.t_day,(2.-3.*sp.p)/4.)*pow(sp.d_L28,-2)*pow(nu_14,-sp.p/2.)}
             }


DICT_BETA = {'D': lambda p: 1./3.,
             'E': lambda p: 1./3.,
             'F': lambda p: -1./2.,
             'G': lambda p: (1.-p)/2.,
             'H': lambda p: -p/2.}


def convert_eV_to_Hz(ev, unit='eV'):
    if unit is 'eV':
        ev = ev
    elif unit is 'keV':
        ev = ev * 1e3
    elif unit is 'MeV':
        ev = ev * 1e6
    elif unit is 'GeV':
        ev = ev * 1e9
    elif unit is 'TeV':
        ev = ev * 1e12
    return ev*e/h


def convert_Hz_to_eV(hz):
    return hz*h/e


def calc_flux(freq, srcpars, pls=None, cbm='ISM', unit='eV'):
    if unit=='GeV':
        freq = freq * 1E9
        unit = 'eV'
    elif unit=='MeV':
        freq = freq * 1E6
        unit = 'eV'
    elif unit=='keV':
        freq = freq * 1E3
        unit = 'eV'

    if unit=='eV':
        freq_eV = freq
        freq = convert_eV_to_Hz(freq_eV)
        logger.debug('{ev} eV equals to {hz} Hz'.format(ev=freq_eV, hz=freq))
    if pls is None:
        pls = find_PLS(freq=freq, srcpars=srcpars, cbm=cbm, unit='Hz')
    logger.debug('Power-law segment: {0}'.format(pls))
    flux_mJy = DICT_FLUX[pls][cbm](freq*1E-14, srcpars) #mJy
    flux = (flux_mJy * u.mJy).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
    return flux


def find_PLS(freq,srcpars,cbm='ISM', unit='eV'):
    logger.debug('CBM: {0}-like'.format(cbm))
    if unit=='GeV':
        freq = freq * 1E9
        unit = 'eV'
    elif unit=='MeV':
        freq = freq * 1E6
        unit = 'eV'
    elif unit=='keV':
        freq = freq * 1E3
        unit = 'eV'
    if unit=='eV':
        freq_eV = freq
        freq = convert_eV_to_Hz(freq_eV, 'eV')
        logger.debug('{ev} eV equals to {hz} Hz'.format(ev=freq_eV, hz=freq))
    dict_freq_break = {}
    dict_PLS = {}
    cool_consistent = None
    for cool in ('Slow', 'Fast'):
        logger.debug('{0}-cooling'.format(cool))
        dict_freq_break[cool] = {}
        dict_PLS[cool] = {}
        for brk in ('m', 'c'):
            dict_freq_break[cool][brk]= DICT_FREQ_BREAK[cool][brk][cbm](srcpars)
            logger.debug('Break frequency {b}: {v}'.format(b=brk, v=dict_freq_break[cool][brk]))
        if dict_freq_break[cool]['m']>=dict_freq_break[cool]['c']:
            if freq>=dict_freq_break[cool]['m']:
                dict_PLS[cool] = 'H'
            elif freq>=dict_freq_break[cool]['c']:
                dict_PLS[cool] = 'F'
            else:
                dict_PLS[cool] = 'E'
            if cool=='Slow':
                logger.debug('{coo} is inconsistent.'.format(coo=cool))
            if cool=='Fast':
                logger.debug('{coo} is consistent.'.format(coo=cool))
                cool_consistent = 'Fast'
        else:
            if freq>=dict_freq_break[cool]['c']:
                dict_PLS[cool] = 'H'
            elif freq>=dict_freq_break[cool]['m']:
                dict_PLS[cool] = 'G'
            else:
                dict_PLS[cool] = 'D'
            if cool=='Slow':
                logger.debug('{coo} is consistent.'.format(coo=cool))
                cool_consistent = 'Slow'
            if cool=='Fast':
                logger.debug('{coo} is inconsistent.'.format(coo=cool))
    if cool_consistent is None:
        logger.critical('No consistent cooling regime!!!')
        sys.exit(1)
    return dict_PLS[cool_consistent]


def integrate_flux(flux, phindex, e0, emin, emax, unit='eV'):
    if unit=='GeV':
        e0 = e0 * 1E9
        emin = emin * 1E9
        emax = emax * 1E9
        unit = 'eV'
    elif unit=='MeV':
        e0 = e0 * 1E6
        emin = emin * 1E6
        emax = emax * 1E6
        unit = 'eV'
    elif unit=='keV':
        e0 = e0 * 1E3
        emin = emin * 1E3
        emax = emax * 1E3
        unit = 'eV'

    if unit=='eV':
        freq0 = convert_eV_to_Hz(e0, 'eV')
    elif unit=='Hz':
        freq0 = e0
    else:
        logger.critical('Input the energies in eV or the frequencies in Hz!!!')
        sys.exit(1)
    vFv = freq0 * flux
    if phindex!=-2:
        return vFv/(phindex+2)*(pow(emax,phindex+2) - pow(emin, phindex+2))/pow(e0, phindex+2)
    else:
        return vFv * np.log(emax/emin)


def make_lightcurve(dict_par, tmin, tmax, emin=1E-1, emax=1E2, unit='GeV'):
    times = 10**np.linspace(np.log10(tmin), np.log10(tmax), int(np.log10(tmax/tmin)*10+1.5)) #10
    energies = 10**np.linspace(np.log10(emin), np.log10(emax), int(np.log10(emax/emin)*10+1.5)) #10
    erefs = np.array([np.sqrt(e0*e1) for e0,e1 in zip(energies[:-1], energies[1:])])
    integral_fluxes = []
    for t in times:
        sp = SRC_PARAMETERS(dict_par['pelec'], dict_par['zredshift'], dict_par['energy'], dict_par['epsilone'], dict_par['epsilonb'], dict_par['nism'], dict_par['awind'], t)
        integral_fluxes.append(0)
        list_pls = []
        for eref, e0, e1 in zip(erefs, energies[:-1], energies[1:]):
            pls = find_PLS(freq=eref, srcpars=sp, cbm='ISM', unit=unit)
            list_pls.append(pls)
            flux = calc_flux(freq=eref, srcpars=sp, pls=pls, cbm='ISM', unit=unit)
            integral_fluxes[-1] += integrate_flux(flux=flux, phindex=-(DICT_BETA[pls](sp.p)+1.), e0=eref, emin=e0, emax=e1, unit=unit)
        logger.info('Tobs = {t} s: Integral flux = {f} erg/cm2/s ({pls})'.format(t=t, f=integral_fluxes[-1], pls=set(list_pls)))
    return (times, np.array(integral_fluxes))


# def make_lightcurve_Sari(dict_par, tmin, tmax, emin=1E-1, emax=1E2, unit='GeV'):
#     times = 10**np.linspace(np.log10(tmin), np.log10(tmax), int(np.log10(tmax/tmin)*10+1.5))
#     energies = 10**np.linspace(np.log10(emin), np.log10(emax), int(np.log10(emax/emin)*10+1.5))
#     erefs = np.array([np.sqrt(e0*e1) for e0,e1 in zip(energies[:-1], energies[1:])])

#     evolsp = EvolvingFluid(dict_par)
        
#     integral_fluxes = []
#     integral_fluxes_IC = []
#     for it, t in enumerate(times):
#         integral_fluxes.append(0)
#         integral_fluxes_IC.append(0)
#         for eref, e0, e1 in zip(erefs, energies[:-1], energies[1:]):
#             flux = evolsp.sp_current(t).f_nu_redshift(t_day=(t*u.s).to(u.d), nu=convert_eV_to_Hz(eref, unit)*u.Hz).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
#             integral_fluxes[-1] += integrate_flux(flux=flux, phindex=-(evolsp.sp_current(t).beta_redshift(t_day=t/86400.*u.day, nu=convert_eV_to_Hz(eref, unit)*u.Hz) + 1.), e0=eref, emin=e0, emax=e1, unit=unit)

#             if dict_par['iceffect'] is True:
#                 flux_IC = evolsp.sp_current(t).f_nu_IC_redshift(t_day=(t*u.s).to(u.d), nu=convert_eV_to_Hz(eref, unit)*u.Hz).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
#                 integral_fluxes_IC[-1] += integrate_flux(flux=flux_IC, phindex=-(evolsp.sp_current(t).beta_IC_redshift(t_day=t/86400.*u.day, nu=convert_eV_to_Hz(eref, unit)*u.Hz) + 1.), e0=eref, emin=e0, emax=e1, unit=unit)

#         logger.debug('Tobs = {t:.2} s: Integral flux = {f:1.2E} erg/cm2/s'.format(t=t, f=integral_fluxes[-1]))
        
#     if dict_par['iceffect'] is True:
#         return (times, dict_report, {'Sync.': np.array(integral_fluxes), 'SSC':np.array(integral_fluxes_IC)})
#     else:
#         return (times, dict_report, {'Sync.': np.array(integral_fluxes)})


def make_snapshot(dict_par, t, emin=1E15, emax=1E26):
    sp = SRC_PARAMETERS(dict_par['pelec'], dict_par['zredshift'], dict_par['energy'], dict_par['epsilone'], dict_par['epsilonb'], dict_par['nism'], dict_par['awind'], t)    
    energies = 10**np.linspace(np.log10(emin), np.log10(emax), int(np.log10(emax/emin)*10+1.5))
    erefs = np.array([np.sqrt(e0*e1) for e0,e1 in zip(energies[:-1], energies[1:])])
    fluxes = np.zeros_like(erefs)
    for ie, (eref, e0, e1) in enumerate(zip(erefs, energies[:-1], energies[1:])):
        pls = find_PLS(freq=eref, srcpars=sp, cbm='ISM', unit='Hz')
        fluxes[ie] = calc_flux(freq=eref, srcpars=sp, pls=pls, cbm='ISM', unit='Hz')
    return (erefs, fluxes)


class EvolvingFluid:
    def __init__(self,dict_par):
        self.dict_par = dict_par
        if self.dict_par['evolution'] == 'adiabatic':
            self.sp = CalcCharaFreq.ShockedFluid_Sari(e52=self.dict_par['energy']/1e52, n0=self.dict_par['nism'], p=self.dict_par['pelec'], epsilonB=self.dict_par['epsilonb'], epsilonE=self.dict_par['epsilone'], dL28=None, z=self.dict_par['zredshift'], iceffect=self.dict_par['iceffect'], eblmodel=self.dict_par['eblmodel'])
        elif self.dict_par['evolution'] == 'radiative':
            self.sp = CalcCharaFreq.ShockedFluid_Sari_radiative(e52=self.dict_par['energy']/1e52, n0=self.dict_par['nism'], p=self.dict_par['pelec'], epsilonB=self.dict_par['epsilonb'], epsilonE=self.dict_par['epsilone'], dL28=None, z=self.dict_par['zredshift'], gamma0=self.dict_par['gamma0'], iceffect=self.dict_par['iceffect'], eblmodel=self.dict_par['eblmodel'])
            # After transition
            #tsols_transit = optimize.fsolve(self.sp.diff_nu_breaks, self.sp.time_transition_F2S)[0]*u.s
            e_residual = self.sp.residual_energy(self.sp.time_transition_F2S * u.s).to(u.erg)
            logger.info('Residual energy at transition time: {0}'.format(e_residual))
            r_transition = self.sp.radius(self.sp.time_transition_F2S*u.s)
            logger.info('Radius at transition time: {0}'.format(r_transition))
            
            logger.info('==== After switching to adiabatic evolution ====')
            self.sp_afterRadiative = CalcCharaFreq.ShockedFluid_Sari_afterRadiative(e52=e_residual.to(u.erg).value/1e52, n0=self.dict_par['nism'], p=self.dict_par['pelec'], epsilonB=self.dict_par['epsilonb'], epsilonE=self.dict_par['epsilone'], dL28=None, z=self.dict_par['zredshift'], gamma0=None, r_residual=r_transition, t_transition=self.sp.time_transition_F2S*u.s, iceffect=self.dict_par['iceffect'], eblmodel=self.dict_par['eblmodel'])

        else:
            logger.critical('Specify "adiabatic" or "radiative" as "evolution"!!!')
            sys.exit(1)

    def sp_current(self, t):
        #if self.dict_par['evolution'] == 'radiative' and self.sp.gamma_c(t)>self.sp.gamma_m(t): #t.to(u.s).value > self.sp.time_transition_redshift_F2S:
            #return 1 #
            #return self.sp_afterRadiative
        #else:
         #   return self.sp
        return self.sp
        
    def report_status(self, t_obs):
        t_burst = t_obs / (1.+self.dict_par['zredshift'])
        dict_status = {}
        sp_current = self.sp_current(t_obs)
        dict_status['gamma2'] = sp_current.gamma2(t_burst.to(u.day)).decompose()
        dict_status['mag'] = sp_current.mag(t_burst.to(u.day)).cgs.value * u.G
        dict_status['gamma_c'] = sp_current.gamma_c(t_burst.to(u.day)).decompose()
        dict_status['gamma_m'] = sp_current.gamma_m(t_burst.to(u.day)).decompose()
        dict_status['nu_obs_c'] = sp_current.nu_c_redshift(t_obs.to(u.day)).to(u.Hz)
        dict_status['nu_obs_m'] = sp_current.nu_m_redshift(t_obs.to(u.day)).to(u.Hz)
        dict_status['nu_obs_IC_c'] = sp_current.nu_IC_c_redshift(t_obs.to(u.day)).to(u.Hz)
        dict_status['nu_obs_IC_m'] = sp_current.nu_IC_m_redshift(t_obs.to(u.day)).to(u.Hz)
        return dict_status

    def report_evol(self, tmin, tmax):
        """Time is in the observer frame.
"""
        times = 10**np.linspace(np.log10(tmin), np.log10(tmax), int(np.log10(tmax/tmin)*10+1.5))

        dict_report = {'gamma2':{'val':np.zeros_like(times), 'unit':u.dimensionless_unscaled}, 
                       'mag':{'val':np.zeros_like(times), 'unit':u.G}, 
                       'gamma_c':{'val':np.zeros_like(times), 'unit':u.dimensionless_unscaled}, 
                       'gamma_m':{'val':np.zeros_like(times), 'unit':u.dimensionless_unscaled}, 
                       'nu_obs_c':{'val':np.zeros_like(times), 'unit':u.Hz}, 
                       'nu_obs_m':{'val':np.zeros_like(times), 'unit':u.Hz}, 
                       'nu_obs_IC_c':{'val':np.zeros_like(times), 'unit':u.Hz}, 
                       'nu_obs_IC_m':{'val':np.zeros_like(times), 'unit':u.Hz}}

        for it, t in enumerate(times):
            report = self.report_status(t*u.s)
            for krepo in dict_report.keys():
                dict_report[krepo]['val'][it] = report[krepo].to(dict_report[krepo]['unit']).value
            if dict_report['gamma_c']['val'][it]>dict_report['gamma_m']['val'][it]: #self.sp_current(t*u.s)==1:
                times = times[:it+1]
                for krepo in dict_report.keys():
                    dict_report[krepo]['val'] = dict_report[krepo]['val'][:it+1]
                logger.warning('Radiative phase finished.')
                break
        return times, dict_report


    def make_lightcurve_Sari(self, tmin, tmax, emin=1E-1, emax=1E2, unit='GeV'):
        times = 10**np.linspace(np.log10(tmin), np.log10(tmax), int(np.log10(tmax/tmin)*10+1.5))#Back to 10
        energies = 10**np.linspace(np.log10(emin), np.log10(emax), int(np.log10(emax/emin)*5+1.5)) #Back to 10
        erefs = np.array([np.sqrt(e0*e1) for e0,e1 in zip(energies[:-1], energies[1:])])

        integral_fluxes = []
        integral_fluxes_IC = []
        for it, t in enumerate(times):
            sp_current = self.sp_current(t*u.s)
            # if sp_current.gamma_c(t)>sp_current.gamma_m(t):
            #     times = times[:it]
            #     logger.warning('Radiative phase finished.')
            #     break
            integral_fluxes.append(0)
            integral_fluxes_IC.append(0)
            for eref, e0, e1 in zip(erefs, energies[:-1], energies[1:]):
                flux = sp_current.f_nu_redshift(t_day=(t*u.s).to(u.d), nu=convert_eV_to_Hz(eref, unit)*u.Hz).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
                integral_fluxes[-1] += integrate_flux(flux=flux, phindex=-(sp_current.beta_redshift(t_day=t/86400.*u.day, nu=convert_eV_to_Hz(eref, unit)*u.Hz) + 1.), e0=eref, emin=e0, emax=e1, unit=unit)

                if self.dict_par['iceffect'] is True:
                    flux_IC = sp_current.f_nu_IC_redshift(t_day=(t*u.s).to(u.d), nu=convert_eV_to_Hz(eref, unit)*u.Hz).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
                    integral_fluxes_IC[-1] += integrate_flux(flux=flux_IC, phindex=-(sp_current.beta_IC_redshift(t_day=t/86400.*u.day, nu=convert_eV_to_Hz(eref, unit)*u.Hz) + 1.), e0=eref, emin=e0, emax=e1, unit=unit)
            logger.debug('Tobs = {t:.2} s: Integral flux = {f:1.2E} erg/cm2/s'.format(t=t, f=integral_fluxes[-1]))
            if sp_current.nu_c_redshift(t*u.s)>sp_current.nu_m_redshift(t*u.s):
                times = times[:it+1]
                logger.warning('Radiative phase finished.')
                break
        
        if self.dict_par['iceffect'] is True:
            return (times, {'Sync.': np.array(integral_fluxes), 'SSC':np.array(integral_fluxes_IC)})
        else:
            return (times, {'Sync.': np.array(integral_fluxes)})


    def calc_spectrum_Sari(self, t, emin=1E5, emax=1E26):
        return self.sp_current(t).calc_spectrum_Sari(t, emin, emax)
        #if self.dict_par['evolution'] == 'radiative' and t > self.sp.time_transition_redshift:
        #    return self.sp_afterRadiative.calc_spectrum_Sari(t, emin, emax)
        #else:
        #    return self.sp.calc_spectrum_Sari(t, emin, emax)


    def calc_spectrum_IC_Sari(self, t, emin=1E5, emax=1E26):
        return self.sp_current(t).calc_spectrum_IC_Sari(t, emin, emax)
        #if self.dict_par['evolution'] == 'radiative' and t > self.sp.time_transition_redshift:
        #    return self.sp_afterRadiative.calc_spectrum_IC_Sari(t, emin, emax)
        #else:
        #    return self.sp.calc_spectrum_IC_Sari(t, emin, emax)
        

#def make_snapshot_Sari(dict_par, t, emin=1E5, emax=1E26): #emin=1E15, emax=1E26):

 #   return sp.calc_spectrum_Sari(t, emin, emax)


def plot_snapshot(ax, energies, fluxes, label=None):
    xlim = [min(energies[0]), max(energies[0])]
    ax.set_xlim(xlim)
    ylim = [1e-37, 1e-24]
    ax.set_ylim(ylim)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'$F_{\nu} \mathrm{[erg \cdot cm^{-2} \cdot s^{-1} \cdot Hz^{-1}]}$')
    for i, (e,f) in enumerate(zip(energies, fluxes)):
        ax.plot(e, f, label=label, c=pMatplot.TPL_COLOR[i], scalex=False, scaley=False)
    ax.set_xscale('log', nonposx='clip')
    ax.set_yscale('log', nonposx='clip')
    ax.grid()
    #ax.plot(energies, energies*fluxes, label=label, c='k', scalex=False, scaley=False)
    ax.fill_between(x=[2.4180e22, 2.4180e25], y1=[ylim[0],ylim[0]], y2=[ylim[1],ylim[1]],facecolor='y',alpha=0.2)
    ax.fill_between(x=[7.2540e16, 2.4180e18], y1=[ylim[0],ylim[0]], y2=[ylim[1],ylim[1]],facecolor='orange',alpha=0.2)
    ax_eV = ax.twiny()
    ax_eV.set_xscale("log", nonposx='clip')
    #ax_eV.set_yscale("log", nonposx='clip')
    ax_eV.set_xlim(convert_Hz_to_eV(np.array(ax.get_xlim())))
    ax_eV.set_xlabel('Energy [eV]', color='g')
    ax_eV.grid(ls='-', lw=0.5, alpha=0.5, axis='x', c='g')
    ax_eV.tick_params(axis='x', colors='g')
    ax.set_ylim(ylim)


@click.command()
@click.argument('grb', type=str)
@click.argument('pelec', type=float)
@click.argument('zredshift', type=float)
@click.argument('energy', type=float)
@click.argument('epsilone', type=float)
@click.argument('epsilonb', type=float)
@click.argument('nism', type=float)
@click.argument('awind', type=float)
@click.option('--tmin', type=float, default=100.)
@click.option('--tmax', type=float, default=100000.)
@click.option('--mode', type=click.Choice(['lightcurve', 'snapshot', 'segment', 'movie']))
@click.option('--especific', type=float, default=1E8, help='In eV. For mode segment.')
@click.option('--latemax', type=float, default=100., help='In GeV.')
@click.option('--xrtdata', type=str, default=None)
@click.option('--latdata', type=str, default=None)
@click.option('--evolution', type=click.Choice(['adiabatic', 'radiative']))
@click.option('--iceffect', is_flag=True)
@click.option('--gamma0', type=float, default=100.)
#@click.argument('dst', nargs=-1)
@click.option('--ebl', default=None, type=click.Choice([None, 'Franceschini08']))
#@click.option('--values', type=(str, int))
#@click.option('--values', multiple=True)
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
#@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(grb,pelec,zredshift,energy,epsilone,epsilonb,nism,awind, tmin, tmax, mode, especific, latemax, xrtdata, latdata, evolution, gamma0, iceffect, ebl, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)
    
    dict_par = {'pelec':pelec, 
                'zredshift':zredshift, 
                'energy': energy, 
                'epsilone': epsilone, 
                'epsilonb': epsilonb, 
                'nism': nism, 
                'awind': awind,
                'evolution': evolution,
                'gamma0': gamma0,
                'iceffect': iceffect,
                'eblmodel': InterpolateEBLmodels.DICT_PATH_MODEL[ebl] if ebl is not None else None}

    # dict_par_name = {'pelec':r'$p$', 
    #                  'zredshift':r'$z$', 
    #                  'energy': r'$E$', 
    #                  'epsilone': , 
    #                  'epsilonb': epsilonb, 
    #                  'nism': nism, 
    #                  'awind': awind,
    #                  'evolution': evolution,
    #                  'gamma0': gamma0,
    #                  'iceffect': iceffect
    #                  }
    evolsp = EvolvingFluid(dict_par)
    if mode == 'lightcurve':
        logger.info('Details of time evolution')
        fig_repo, ax_repo = plt.subplots(2, 2, figsize=(14, 10))

        time_repo, dict_repo = evolsp.report_evol(tmin, tmax)

        nax_repo = 0
        ax_repo[nax_repo/2][nax_repo%2].set_xlabel('Time [s]')
        ax_repo[nax_repo/2][nax_repo%2].set_ylabel(r'$\gamma_{2}$')
        #ax_repo[nax_repo/2][nax_repo%2].set_ylim((1, 1e5))
        ax_repo[nax_repo/2][nax_repo%2].set_xscale('log')
        ax_repo[nax_repo/2][nax_repo%2].set_yscale('log')
        ax_repo[nax_repo/2][nax_repo%2].grid()
        #ax_repo[nax_repo/2][nax_repo%2].grid(which='minor', axis='y')
        ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['gamma2']['val'], lw=2)
        
        nax_repo+=1
        ax_repo[nax_repo/2][nax_repo%2].set_xlabel('Time [s]')
        ax_repo[nax_repo/2][nax_repo%2].set_ylabel(r'$B \ \mathrm{[G]}$')
        ax_repo[nax_repo/2][nax_repo%2].set_xscale('log')
        ax_repo[nax_repo/2][nax_repo%2].set_yscale('log')
        ax_repo[nax_repo/2][nax_repo%2].grid()
        #ax_repo[nax_repo/2][nax_repo%2].grid(which='minor', axis='y')
        ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['mag']['val'], lw=2)
        
        nax_repo+=1
        ax_repo[nax_repo/2][nax_repo%2].set_xlabel('Time [s]')
        ax_repo[nax_repo/2][nax_repo%2].set_ylabel(r'$\gamma_{e}$')
        #ax_repo[nax_repo/2][nax_repo%2].set_ylim((1, 1e5))
        ax_repo[nax_repo/2][nax_repo%2].set_xscale('log')
        ax_repo[nax_repo/2][nax_repo%2].set_yscale('log')
        ax_repo[nax_repo/2][nax_repo%2].grid()
        #ax_repo[nax_repo/2][nax_repo%2].grid(which='minor', axis='y')
        ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['gamma_c']['val'], label=r'$\gamma_{c}$', lw=2, c=DICT_COLOR['c'])
        ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['gamma_m']['val'], label=r'$\gamma_{m}$', lw=2, c=DICT_COLOR['m'])
        ax_repo[nax_repo/2][nax_repo%2].legend(loc=0)

        nax_repo+=1
        ax_repo[nax_repo/2][nax_repo%2].set_xlabel('Time [s]')
        ax_repo[nax_repo/2][nax_repo%2].set_ylabel(r'$\nu \ \mathrm{[Hz]}$')
        #ax_repo[nax_repo/2][nax_repo%2].set_ylim((1, 1e5))
        ax_repo[nax_repo/2][nax_repo%2].set_xscale('log')
        ax_repo[nax_repo/2][nax_repo%2].set_yscale('log')
        ax_repo[nax_repo/2][nax_repo%2].grid()
        ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['nu_obs_c']['val'], label=r'$\nu_{c}$', lw=2, c=DICT_COLOR['c'])
        ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['nu_obs_m']['val'], label=r'$\nu_{m}$', lw=2, c=DICT_COLOR['m'])
        if iceffect is True:
            ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['nu_obs_IC_c']['val'], label=r'$\nu^{IC}_{c}$', lw=2, ls=DICT_LINE['SSC'], c=DICT_COLOR['c'])
            ax_repo[nax_repo/2][nax_repo%2].plot(time_repo, dict_repo['nu_obs_IC_m']['val'], label=r'$\nu^{IC}_{m}$', lw=2, ls=DICT_LINE['SSC'], c=DICT_COLOR['m'])

        ax_repo_energy = ax_repo[nax_repo/2][nax_repo%2].twinx()
        ax_repo_energy.set_yscale("log", nonposx='clip')
        ax_repo_energy.set_ylim((np.array(ax_repo[nax_repo/2][nax_repo%2].get_ylim())*u.Hz * constants.h).to(u.eV).value)
        ax_repo_energy.set_ylabel('Energy [eV]', color='magenta')
        ax_repo_energy.grid(ls='-', lw=0.5, alpha=0.5, axis='y', c='magenta')
        ax_repo_energy.tick_params(axis='y', colors='magenta')

        ax_repo[nax_repo/2][nax_repo%2].legend(loc=0)

        fig_repo.tight_layout()
        #fig_repo.subplots_adjust(hspace=0)
        for ff in ('png', 'pdf'):
            fig_repo.savefig('./SynchrotronModelEvolReport_Sari98{ic}{ebl}.{fmt}'.format(fmt=ff,ic='_wIC' if iceffect is True else '', ebl='_'+ebl if ebl is not None else ""))

        logger.info('LAT light curve')
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_axes((0.15, 0.15, 0.75, 0.75))
        ax.grid()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim((1e-14, 1e-4))

        emin_lat, emax_lat = 1E-1, latemax # GeV
        # Observed data table
        tb_lat = read_observation_table(latdata)
        obslat = ObservedLightCurve(name='LAT', title='LAT data ({0:.1f} - {1:.0f} GeV)'.format(emin_lat, emax_lat))
        obslat.load(qtable=tb_lat, colx='time', coly='eflux', colxerr=['time_err_lo', 'time_err_hi'], colyerr=['eflux_err_lo', 'eflux_err_hi'], selector=None)
        obslat.plot(ax)
        # Model
        times_lat, fluxes_lat = evolsp.make_lightcurve_Sari(tmin, tmax, emin=emin_lat, emax=emax_lat, unit='GeV')
        curve_lat_sum = np.zeros_like(times_lat)
        for key, curve in fluxes_lat.iteritems():
            ax.plot(times_lat, curve, ls=DICT_LINE[key], label='LAT model ({0})'.format(key), c=DICT_COLOR['LAT'], lw=1, scaley=False)
            curve_lat_sum += curve
        if len(fluxes_lat.keys())>1:
            ax.plot(times_lat, curve_lat_sum, label='LAT model (sum)', c=DICT_COLOR['LAT'], lw=2, ls='-') 
        #plot_lightcurve(ax, times_lat, fluxes_lat, label=('LAT (Sync.)','LAT (SSC)'), color=DICT_COLOR['LAT'])

        logger.info('XRT light curve')
        emin_xrt, emax_xrt = 0.3, 10. # keV
        # Observed data table
        tb_xray = read_observation_table(xrtdata)
        obsxrt = ObservedLightCurve(name='XRT', title='XRT data ({0:.1f} - {1:.0f} keV)'.format(emin_xrt, emax_xrt))
        obsxrt.load(qtable=tb_xray, colx='Time[s]', coly='Flux[erg/cm^2/s]', colxerr=['T_-ve[s]', 'T_+ve[s]'], colyerr=['Fluxneg[erg/cm^2/s]', 'Fluxpos[erg/cm^2/s]'], selector=None)
        obsxrt.plot(ax)
        # Model
        times_xrt, fluxes_xrt = evolsp.make_lightcurve_Sari(tmin, tmax, emin=emin_xrt, emax=emax_xrt, unit='keV')
        curve_xrt_sum = np.zeros_like(times_xrt)
        for key, curve in fluxes_xrt.iteritems():
            ax.plot(times_xrt, curve, ls=DICT_LINE[key], label='XRT model ({0})'.format(key), c=DICT_COLOR['XRT'], lw=1, scaley=False)
            curve_xrt_sum += curve
        if len(fluxes_xrt.keys())>1:
            ax.plot(times_xrt, curve_xrt_sum, label='XRT model (sum)', c=DICT_COLOR['XRT'], lw=2, ls='-') 

        ax.legend(loc=0, fontsize=8)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(r'Integral energy flux $\mathrm{[erg/cm^2 s]}$')
        for ff in ('png', 'pdf'):
            fig.savefig('./SynchrotronModelLightCurves_Sari98{ic}{ebl}.{fmt}'.format(fmt=ff,ic='_wIC' if iceffect is True else '', ebl='_'+ebl if ebl is not None else ""))

    elif mode == 'snapshot':
        logger.info('Spectrum snapshot')
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_axes((0.15, 0.15, 0.75, 0.75))
        #evolsp = EvolvingFluid(dict_par)
        energies = []
        fluxes = []
        spec_Sari = evolsp.calc_spectrum_Sari(tmin) #make_snapshot_Sari(dict_par, tmin)
        energies.append(spec_Sari[0])
        fluxes.append(spec_Sari[1])
        if iceffect is True:
            spec_Sari_IC = evolsp.calc_spectrum_IC_Sari(tmin)
            energies.append(spec_Sari_IC[0])
            fluxes.append(spec_Sari_IC[1])
        plot_snapshot(ax, energies, fluxes)
        # ax_eV = ax.twiny()
        # ax_eV.set_xscale("log", nonposx='clip')
        # ax_eV.set_yscale("log", nonposx='clip')
        # ax_eV.set_xlim(convert_Hz_to_eV(np.array(ax.get_xlim())))
        # ax_eV.set_xlabel('Energy [eV]', color='g')
        # ax_eV.grid(ls='-', lw=0.5, alpha=0.5, axis='x', c='g')
        # ax_eV.tick_params(axis='x', colors='g')
        for ff in ('png', 'pdf'):
            fig.savefig('./SynchrotronModelSnapshot_Sari98_T{t:0>8.0f}{ic}{ebl}.{fmt}'.format(fmt=ff, t=tmin, ic='_wIC' if iceffect is True else '', ebl='_'+ebl if ebl is not None else ""))

    elif mode == 'movie':
        logger.info('Spectrum movie')
        fig = plt.figure(figsize=(6, 5))
        times = 10**np.linspace(np.log10(tmin), np.log10(tmax), 51)
        snapshots = []
        evolsp = EvolvingFluid(dict_par)
        def animate(iframe):
            logger.debug('Frame: No. {0}, Time: {1} s:'.format(iframe, times[iframe]))
            fig.clf() 
            ax = fig.add_axes((0.15, 0.25, 0.75, 0.65))
            ax_tbar = fig.add_axes((0.15, 0.075, 0.75, 0.05))
            ax_tbar.barh(bottom=0, width=[times[iframe]], left=0, height=0.025, log=True, tick_label='Time [s]')
            ax_tbar.set_xlim(min(times), max(times))
            energies = []
            fluxes = []
            spec_Sari = evolsp.calc_spectrum_Sari(times[iframe])
            energies.append(spec_Sari[0])
            fluxes.append(spec_Sari[1])
            if iceffect is True:
                spec_Sari_IC = evolsp.calc_spectrum_IC_Sari(times[iframe])
                energies.append(spec_Sari_IC[0])
                fluxes.append(spec_Sari_IC[1])
            #energies, fluxes = evolsp.calc_spectrum_Sari(tmin)
            #energies, fluxes = make_snapshot_Sari(dict_par, times[iframe])
            plot_snapshot(ax, energies, fluxes)

        anim = animation.FuncAnimation(fig, animate, frames=len(times), interval=200)
        anim.save('./SynchrotronModelMovie_Sari98_T{tmin:0>8.0f}-{tmax:0>8.0f}.{fmt}'.format(fmt='gif', tmin=tmin, tmax=tmax), writer='imagemagick', fps=4)


    elif mode == 'segment':
        sp = SRC_PARAMETERS(dict_par['pelec'], dict_par['zredshift'], dict_par['energy'], dict_par['epsilone'], dict_par['epsilonb'], dict_par['nism'], dict_par['awind'], tmin)
        pls = find_PLS(freq=especific,srcpars=sp,cbm='ISM', unit='eV')
        logger.info('Power-law segment: {0} at T = {1} s, E = {2} eV (nu = {3} Hz)'.format(pls, tmin, especific, convert_eV_to_Hz(especific)))


if __name__ == '__main__':
    main()
