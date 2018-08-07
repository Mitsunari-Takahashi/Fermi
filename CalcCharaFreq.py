#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import sqrt, cos, sin, tan, acos, asin, atan, radians, degrees, pi, log10
from astropy import units as u
from astropy.units import cds
import astropy.constants as constants
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import Distance
#from astropy.cosmology import WMAP5, WMAP7
from astropy.cosmology import WMAP9 as cosmo
from scipy import optimize
import click
from ROOT import TGraph2D
import InterpolateEBLmodels
#import pEBLatten
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL


##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


COSMO = FlatLambdaCDM(H0=70, Om0=0.3)


class ShockedFluid:
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., iceffect=False, eblmodel=None):
        self.energy = e52 * 1E52 * u.erg
        self.ndensity = n0 * u.cm**-3
        self.p = p
        self.epsilonB = epsilonB
        self.epsilonE = epsilonE
        self.epsilonE_bar = epsilonE * (self.p-2.) / (self.p-1.)
        if dL28==None and z==None:
            logger.error('Distance and redshift are neither provided!!!')
            dL28 = 5.
        if dL28 is not None and z==None:
            self.dL = dL28 * 1E28 * u.cm
            distance = Distance(value=self.dL)
            self.z = distance.z
        elif dL28 is not None and z is not None:
            self.dL = dL28 * 1E28 * u.cm
            self.z = z
        else:
            self.z = z
            self.dL = cosmo.luminosity_distance(self.z)
        logger.debug('EBL model: {0}'.format(eblmodel))
        self.ebl = InterpolateEBLmodels.read_model(eblmodel)[0] if eblmodel is not None else None
        self.gamma0 = gamma0
        self.iceffect = iceffect


    def ebl_atten(self, nu):
        e_TeV = (constants.h * nu).to(u.TeV).value
        if e_TeV>=0.01:
            tau = self.ebl.Interpolate(self.z, e_TeV)
            atten = np.exp(-tau)
        else:
            atten = 1.
        return atten


    def info(self):
        logger.info('Energy: {0}'.format(self.energy))
        logger.info('Number density: {0}'.format(self.ndensity))
        logger.info('Electron index: {0}'.format(self.p))
        logger.info('epsilon_B: {0}'.format(self.epsilonB))
        logger.info('epsilon_e: {0}'.format(self.epsilonE))
        logger.info('epsilon_e_bar: {0}'.format(self.epsilonE_bar))
        logger.info('Luminosity distance: {0}'.format(self.dL))
        logger.info('Redshift: {0}'.format(self.z))
        logger.info('Initial Lorentz factor of shock: {0}'.format(self.gamma0))


    def hnu_c(self, t_day):
        return constants.h * self.nu_c(t_day)


    def hnu_m(self, t_day):
        return constants.h * self.nu_m(t_day)


    def hnu_c_redshift(self, t_day):
        return constants.h * self.nu_c_redshift(t_day)


    def hnu_m_redshift(self, t_day):
        return constants.h * self.nu_m_redshift(t_day)


class ShockedFluid_Sari(ShockedFluid):
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., iceffect=False, eblmodel=None):
        ShockedFluid.__init__(self, e52, n0, p, epsilonB, epsilonE, dL28, z, gamma0, iceffect, eblmodel)
        if self.iceffect is True:
            time_transitions = self.get_time_transition()
            time_transitions_redshift = self.get_time_transition_redshift()
            logger.debug('Transition time: {0}'.format(time_transitions_redshift))
            self.time_transition_F2S = time_transitions[0]
            self.time_transition_redshift_F2S = time_transitions_redshift[0]
            self.time_transition_C2S = time_transitions[1]
            self.time_transition_redshift_C2S = time_transitions_redshift[1]
        else:
            self.time_transition_F2S = self.get_time_transition()
            self.time_transition_redshift_F2S = self.get_time_transition_redshift()

   # Common for adiabatic and radiative

    def gamma2(self, t_day):
        gamma2 = ((self.radius(t_day) / 4. / constants.c / t_day)**0.5).decompose() #( 17. * self.energy / (16. * pi * self.ndensity * u.M_p.cgs * cds.c.cgs**2 * self.radius(t_day)**3) )**0.5

        if gamma2.unit != u.dimensionless_unscaled and not isinstance(gamma2, float):
            logger.error('gamma2 is not in dimensionless!')
            logger.error('Radius: {0}'.format(self.radius(t_day)))
            logger.error('Radius [cm]: {0}'.format(self.radius(t_day).to(u.cm).value))
            logger.error('t_day: {0}'.format(t_day))
            logger.error('gamma2: {0}'.format(gamma2))
            logger.error('Unit of gamma2: {0}'.format(gamma2.unit))
            #sys.exit(1)
        return gamma2


    def convert_gamma_to_nu(self, gamma, t_day):
        # if gamma.decompose().unit == u.dimensionless_unscaled:
        #     logger.error('Input gamma is not dimensionless.')
        #     logger.error(gamma)
        #     sys.exit(1)
        nu = self.gamma2(t_day) * gamma**2 * constants.e.esu.value * self.mag(t_day).cgs.value / (2.*pi*constants.m_e.cgs.value*constants.c.cgs.value) * u.Hz
        if nu.decompose().unit != u.Hz:
            logger.error('nu is not in Hz!!')
            logger.error('gamma2: {0}'.format(self.gamma2(t_day)))
            logger.error('gamma_c: {0}'.format(self.gamma_c(t_day)))
            logger.error('gamma_m: {0}'.format(self.gamma_m(t_day)))
            logger.error('Input gamma: {0}'.format(gamma))
            logger.error('Frequency: {0}'.format(nu))
            #sys.exit(1)
        return nu


    def mag(self, t_day):
        return (32. * pi * u.M_p.cgs * self.epsilonB * self.ndensity)**0.5 * self.gamma2(t_day) * cds.c


    def gamma_c_reduce_factot(self, t_day):
        if self.iceffect == True:
            logger.debug('Transition time from fast-cooling to slow-cooling: {0}'.format(self.time_transition_redshift_F2S))
            logger.debug('Transition time from IC-dominant to synchrotron-dominant: {0}'.format(self.time_transition_redshift_C2S))
            DICT_RED_FACTOR = {'Fast': (1.+(1.+4.*self.epsilonE/self.epsilonB)**0.5)/2.,
                               'IC-dominant Slow': (self.epsilonE/self.epsilonB)**0.5 * (t_day.to(u.s).value/self.time_transition_redshift_F2S)**(-(self.p-2.)/2./(4.-self.p)),
                               'Synchrotron-dominant Slow': 1.
                               }
            red_factor = DICT_RED_FACTOR[self.get_cooling_phase_redshift(t_sec=t_day.to(u.s).value)]
            if not isinstance(red_factor, float):
                logger.warning('Reduce factor: {0}'.format(red_factor))
                logger.warning('Transit time from fast-cooling to slow-cooling: {0}'.format(self.time_transition_redshift_F2S))
                logger.warning('Transit time from IC-dominant to synchrotron-dominant: {0}'.format(self.time_transition_redshift_C2S))
                logger.warning('Cooling phase: {0}'.format(self.get_cooling_phase_redshift(t_sec=t_day.to(u.s).value)))
        else:
            red_factor = 1.
        return red_factor


    def gamma_c(self, t_day):
        gamma_c = (3.*constants.m_e.cgs / (16.*self.epsilonB*constants.sigma_T.cgs*constants.m_p.cgs*constants.c * self.gamma2(t_day)**3 * self.ndensity * t_day.to(u.s)) / self.gamma_c_reduce_factot(t_day)).decompose()
        if gamma_c.unit != u.dimensionless_unscaled:
            logger.error('gamma_c is not dimensionless.')
            logger.error('gamma_c: {0}'.format(gamma_c))
            logger.error('gamma2: {0}'.format(self.gamma2(t_day)))
            logger.error('gamma_c reduced by: {0}'.format(self.gamma_c_reduce_factot(t_day)))

        #     sys.exit(1)
        return gamma_c


    def gamma_m(self, t_day):
        return self.epsilonE_bar * u.M_p.cgs/u.M_e.cgs * self.gamma2(t_day)


    def diff_nu_breaks(self, t_sec):
        """nu_m - nu_c"""
        return log10(self.nu_m(t_sec*u.s).to(u.Hz).value)-log10(self.nu_c(t_sec*u.s).to(u.Hz).value)


    def diff_nu_redshift_breaks(self, t_sec):
        """nu_m - nu_c"""
        return log10(self.nu_m_redshift(t_sec*u.s).to(u.Hz).value)-log10(self.nu_c_redshift(t_sec*u.s).to(u.Hz).value)


    # def radiative_energy_fraction(self, t_day):
    #     ratio = self.gamma_c(t_day) / self.gamma_m(t_day)
    #     if ratio < 1: # Fast-cooling
    #         return 1.
    #     else: # Slow-cooling
    #         return ratio**(-self.p+2.)


    # def luminosity_ratio(self, t_day):
    #     eta = self.radiative_energy_fraction(t_day)
    #     return (-1.+(1.+4.*eta*self.epsilonE/self.epsilonB)**0.5) / 2.



    def nu_c(self, t_day):
        return self.convert_gamma_to_nu(self.gamma_c(t_day), t_day).to(u.Hz)


    def nu_m(self, t_day):
        return self.convert_gamma_to_nu(self.gamma_m(t_day), t_day).to(u.Hz)


    # SSC
    def nu_IC_c(self, t_day):
        return 2.*self.gamma_c(t_day)**2 * self.nu_c(t_day)


    def nu_IC_m(self, t_day):
        return 2.*self.gamma_m(t_day)**2 * self.nu_m(t_day)


    def f_nu_IC_max(self, t_day):
        return constants.sigma_T.cgs * self.ndensity * self.radius(t_day) / 3. * self.f_nu_max(t_day)


    # def nu_c(self, t_day):
    #     return 9./4.*sqrt(2./17.) * constants.e.esu * u.M_e.cgs/u.M_p.cgs * cds.c**0.5 / constants.sigma_T.cgs**2 * self.epsilonB**-1.5 * self.ndensity**-1 * (t_day.to(u.s))**-0.5 * self.energy**-0.5 # * (1.+self.z)**-0.5

    
    # def nu_m(self, t_day):
    #     return ( 17. / 128. * self.energy * self.epsilonB / cds.c**5 / (t_day.to(u.s))**3. )**0.5 * self.epsilonE**2 * ((self.p-2.)/(self.p-1.))**2 * (u.M_p.cgs/u.M_e.cgs)**2 * constants.e.esu / u.M_e.cgs / pi # * (1.+self.z)**0.5


    # def f_nu_max(self, t_day):
    #     return 17./18. * (self.ndensity*self.epsilonB/2./pi/u.M_p.cgs)**0.5 * self.energy * u.M_e.cgs * cds.c * constants.sigma_T.cgs / constants.e.esu / self.dL**2 # * (1.+self.z)


    def f_nu(self, t_day, nu):
        if self.nu_c(t_day) <= self.nu_m(t_day):
            # Fast cooling
            if nu<self.nu_c(t_day):
                f = (nu/self.nu_c(t_day))**(1./3.) * self.f_nu_max(t_day)
            elif nu<self.nu_m(t_day):
                f = (nu/self.nu_c(t_day))**(-1./2.) * self.f_nu_max(t_day)
            else:
                f = (nu/self.nu_m(t_day))**(-self.p/2.) * (self.nu_m(t_day)/self.nu_c(t_day))**(-1./2.) * self.f_nu_max(t_day)
        else:
            # Slow cooling
            if nu<self.nu_m(t_day):
                f = (nu/self.nu_m(t_day))**(1./3.) * self.f_nu_max(t_day)
            elif nu<self.nu_c(t_day):
                f = (nu/self.nu_m(t_day))**(-(self.p-1.)/2.) * self.f_nu_max(t_day)
            else:
                f = (nu/self.nu_c(t_day))**(-self.p/2.) * (self.nu_c(t_day)/self.nu_m(t_day))**(-(self.p-1.)/2.) * self.f_nu_max(t_day)
        return f * self.ebl_atten(nu)


    def nu_c_redshift(self, t_day):
        return self.nu_c(t_day/(1.+self.z)) / (1.+self.z)


    def nu_m_redshift(self, t_day):
        return self.nu_m(t_day/(1.+self.z)) / (1.+self.z)


    def f_nu_max_redshift(self, t_day):
        return self.f_nu_max(t_day/(1.+self.z)) * (1.+self.z)


    def f_nu_redshift(self, t_day, nu):
        if self.nu_c_redshift(t_day) <= self.nu_m_redshift(t_day):
            # Fast cooling
            if nu<self.nu_c_redshift(t_day):
                f = (nu/self.nu_c_redshift(t_day))**(1./3.) * self.f_nu_max_redshift(t_day)
            elif nu<self.nu_m(t_day):
                f = (nu/self.nu_c_redshift(t_day))**(-1./2.) * self.f_nu_max_redshift(t_day)
            else:
                f = (nu/self.nu_m_redshift(t_day))**(-self.p/2.) * (self.nu_m_redshift(t_day)/self.nu_c_redshift(t_day))**(-1./2.) * self.f_nu_max_redshift(t_day)
        else:
            # Slow cooling
            if nu<self.nu_m_redshift(t_day):
                f = (nu/self.nu_m_redshift(t_day))**(1./3.) * self.f_nu_max_redshift(t_day)
            elif nu<self.nu_c(t_day):
                f = (nu/self.nu_m_redshift(t_day))**(-(self.p-1.)/2.) * self.f_nu_max_redshift(t_day)
            else:
                f = (nu/self.nu_c_redshift(t_day))**(-self.p/2.) * (self.nu_c_redshift(t_day)/self.nu_m_redshift(t_day))**(-(self.p-1.)/2.) * self.f_nu_max_redshift(t_day)
        return f * self.ebl_atten(nu)


    # SSC
    def f_nu_IC(self, t_day, nu):
        if self.nu_c(t_day) <= self.nu_m(t_day):
            # Fast cooling
            if nu<self.nu_IC_c(t_day):
                f = (nu/self.nu_IC_c(t_day))**(1./3.) * self.f_nu_IC_max(t_day)
            elif nu<self.nu_IC_m(t_day):
                f = (nu/self.nu_IC_c(t_day))**(-1./2.) * self.f_nu_IC_max(t_day)
            else:
                f = (nu/self.nu_IC_m(t_day))**(-self.p/2.) * (self.nu_IC_m(t_day)/self.nu_IC_c(t_day))**(-1./2.) * self.f_nu_IC_max(t_day)
        else:
            # Slow cooling
            if nu<self.nu_IC_m(t_day):
                f = (nu/self.nu_IC_m(t_day))**(1./3.) * self.f_nu_IC_max(t_day)
            elif nu<self.nu_IC_c(t_day):
                f = (nu/self.nu_IC_m(t_day))**(-(self.p-1.)/2.) * self.f_nu_IC_max(t_day)
            else:
                f = (nu/self.nu_IC_c(t_day))**(-self.p/2.) * (self.nu_IC_c(t_day)/self.nu_IC_m(t_day))**(-(self.p-1.)/2.) * self.f_nu_IC_max(t_day)
        return f * self.ebl_atten(nu)


    def nu_IC_c_redshift(self, t_day):
        return self.nu_IC_c(t_day/(1.+self.z)) / (1.+self.z)


    def nu_IC_m_redshift(self, t_day):
        return self.nu_IC_m(t_day/(1.+self.z)) / (1.+self.z)


    def f_nu_IC_max_redshift(self, t_day):
        return self.f_nu_IC_max(t_day/(1.+self.z)) * (1.+self.z)


    def f_nu_IC_redshift(self, t_day, nu):
        if self.nu_c_redshift(t_day) <= self.nu_m_redshift(t_day):
            # Fast cooling
            if nu<self.nu_IC_c_redshift(t_day):
                f = (nu/self.nu_IC_c_redshift(t_day))**(1./3.) * self.f_nu_IC_max_redshift(t_day)
            elif nu<self.nu_IC_m(t_day):
                f = (nu/self.nu_IC_c_redshift(t_day))**(-1./2.) * self.f_nu_IC_max_redshift(t_day)
            else:
                f = (nu/self.nu_IC_m_redshift(t_day))**(-self.p/2.) * (self.nu_IC_m_redshift(t_day)/self.nu_IC_c_redshift(t_day))**(-1./2.) * self.f_nu_IC_max_redshift(t_day)
        else:
            # Slow cooling
            if nu<self.nu_IC_m_redshift(t_day):
                f = (nu/self.nu_IC_m_redshift(t_day))**(1./3.) * self.f_nu_IC_max_redshift(t_day)
            elif nu<self.nu_IC_c(t_day):
                f = (nu/self.nu_IC_m_redshift(t_day))**(-(self.p-1.)/2.) * self.f_nu_IC_max_redshift(t_day)
            else:
                f = (nu/self.nu_IC_c_redshift(t_day))**(-self.p/2.) * (self.nu_IC_c_redshift(t_day)/self.nu_IC_m_redshift(t_day))**(-(self.p-1.)/2.) * self.f_nu_IC_max_redshift(t_day)
        return f * self.ebl_atten(nu)

    # SSC end

    def beta(self, t_day, nu):
        if self.nu_c(t_day) <= self.nu_m(t_day):
            # Fast cooling
            if nu<self.nu_c(t_day):
                return -1./3.
            elif nu<self.nu_m(t_day):
                return 1./2.
            else:
                return self.p/2.
        else:
            # Slow cooling
            if nu<self.nu_m(t_day):
                return -1./3.
            elif nu<self.nu_c(t_day):
                return (self.p-1.)/2.
            else:
                return self.p/2.


    def beta_redshift(self, t_day, nu):
        if self.nu_c_redshift(t_day) <= self.nu_m_redshift(t_day):
            # Fast cooling
            if nu<self.nu_c_redshift(t_day):
                return -1./3.
            elif nu<self.nu_m_redshift(t_day):
                return 1./2.
            else:
                return self.p/2.
        else:
            # Slow cooling
            if nu<self.nu_m_redshift(t_day):
                return -1./3.
            elif nu<self.nu_c_redshift(t_day):
                return (self.p-1.)/2.
            else:
                return self.p/2.


    #SSC
    def beta_IC(self, t_day, nu):
        if self.nu_c(t_day) <= self.nu_m(t_day):
            # Fast cooling
            if nu<self.nu_IC_c(t_day):
                return -1./3.
            elif nu<self.nu_IC_m(t_day):
                return 1./2.
            else:
                return self.p/2.
        else:
            # Slow cooling
            if nu<self.nu_IC_m(t_day):
                return -1./3.
            elif nu<self.nu_IC_c(t_day):
                return (self.p-1.)/2.
            else:
                return self.p/2.


    def beta_IC_redshift(self, t_day, nu):
        if self.nu_c_redshift(t_day) <= self.nu_m_redshift(t_day):
            # Fast cooling
            if nu<self.nu_IC_c_redshift(t_day):
                return -1./3.
            elif nu<self.nu_IC_m_redshift(t_day):
                return 1./2.
            else:
                return self.p/2.
        else:
            # Slow cooling
            if nu<self.nu_IC_m_redshift(t_day):
                return -1./3.
            elif nu<self.nu_IC_c_redshift(t_day):
                return (self.p-1.)/2.
            else:
                return self.p/2.

    
    def calc_spectrum_Sari(self, t, emin, emax):
        energies = 10**np.linspace(np.log10(emin), np.log10(emax), int(np.log10(emax/emin)*10+1.5))
        erefs = np.array([np.sqrt(e0*e1) for e0,e1 in zip(energies[:-1], energies[1:])])
        fluxes = np.zeros_like(erefs)
        for ie, (eref, e0, e1) in enumerate(zip(erefs, energies[:-1], energies[1:])):
            fluxes[ie] = self.f_nu_redshift(t_day=(t*u.s).to(u.d), nu=eref*u.Hz).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
        return (erefs, fluxes)    


    def calc_spectrum_IC_Sari(self, t, emin, emax):
        energies = 10**np.linspace(np.log10(emin), np.log10(emax), int(np.log10(emax/emin)*10+1.5))
        erefs = np.array([np.sqrt(e0*e1) for e0,e1 in zip(energies[:-1], energies[1:])])
        fluxes = np.zeros_like(erefs)
        for ie, (eref, e0, e1) in enumerate(zip(erefs, energies[:-1], energies[1:])):
            fluxes[ie] = self.f_nu_IC_redshift(t_day=(t*u.s).to(u.d), nu=eref*u.Hz).to(u.erg / u.s / u.cm / u.cm / u.Hz).value
        return (erefs, fluxes)    


    def get_cooling_phase_redshift(self, t_sec):
        if self.iceffect==False:
            if t_sec<self.time_transition_redshift_F2S:
                return 'Fast'
            else:
                return 'Slow'
        else:
            if t_sec<self.time_transition_redshift_F2S:
                return 'Fast'
            elif t_sec<self.time_transition_redshift_C2S:
                return 'IC-dominant Slow'
            else:
                return 'Synchrotron-dominant Slow'


    # Adiabatic evolution
    def f_nu_max(self, t_day):
        return u.M_e.cgs * cds.c.cgs**2 * constants.sigma_T.cgs * self.ndensity * self.radius(t_day)**3 * self.mag(t_day) * self.gamma2(t_day) / 9. / constants.e.esu / self.dL**2


    def get_time_transition(self):
        """Transition time from fast-cooling to slow-cooling in sec in the GRB frame"""
        t_f2s = 210. * self.epsilonB**2 * (self.epsilonE_bar*3.)**2 * (self.energy.to(u.erg).value/1e52) * (self.ndensity.to(u.cm**-3).value) * 86400. # No IC effect
        if self.iceffect is not True:
            return t_f2s
        else:
            t_f2sIC = (self.epsilonE/self.epsilonB) * t_f2s #(self.epsilonE/self.epsilonB)**2 * t_f2s
            t_sIC2sSync = (self.epsilonE/self.epsilonB)**((4.-self.p)/(self.p-2.)) * t_f2sIC
            return  np.array((t_f2sIC, t_sIC2sSync))


    def get_time_transition_redshift(self):
        return self.get_time_transition() * (1.+self.z)


    def radius(self, t_day):
        return ( (17. * self.energy * t_day.to(u.s) / 4. / pi / u.M_p.cgs / self.ndensity / cds.c.cgs)**0.25 ).to(u.cm)


class ShockedFluid_Sari_afterRadiative(ShockedFluid_Sari):
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., r_residual = 0, t_transition = 0, iceffect=False, eblmodel=None):
        ShockedFluid_Sari.__init__(self, e52, n0, p, epsilonB, epsilonE, dL28, z, gamma0, iceffect, eblmodel)
        self.r_residual = r_residual
        self.t_transition = t_transition


    def radius(self, t_day):
        def func(r):
            return r**4 - (self.r_residual.to(u.cm).value)**3 * r - 17.*self.energy.to(u.erg).value*(t_day.to(u.s).value)/(4.*pi*(self.ndensity.to(u.cm**-3).value)*u.M_p.cgs.to(u.g)*(cds.c.to(u.cm/u.s)))

        solutions = optimize.fsolve(func, self.r_residual.to(u.cm).value * ((t_day.to(u.s)/self.t_transition.to(u.s)))**0.25)#, xtol=1e-3)
        #print solutions
        return solutions[0]*u.cm

        #r = ( (17. * self.energy * (t_day.to(u.s)) / 4. / pi / u.M_p.cgs / self.ndensity / cds.c.cgs)**0.25 ).to(u.cm)




class ShockedFluid_Sari_radiative(ShockedFluid_Sari):
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., r_residual = 0, t_transition = 0, iceffect=False, eblmodel=None):
        ShockedFluid_Sari.__init__(self, e52, n0, p, epsilonB, epsilonE, dL28, z, gamma0, iceffect, eblmodel)

    def length(self):
        return ((17.*self.energy/self.gamma0/(cds.c.cgs)**2 / 16./pi/u.M_p.cgs/self.ndensity).to(u.cm**3))**(1./3.)


    def radius(self, t_day):
        return (4.*cds.c.cgs*t_day.to(u.s)/self.length())**(1./7.) * self.length()


    #def gamma2(self, t_day):
    #    return (4.*cds.c.cgs*t_day.to(u.s)/self.length())**(-3./7.)


    #def nu_c(self, t_day):
    #    return self.convert_gamma_to_nu(self.gamma_c(t_day), t_day).to(u.Hz)


    #def nu_m(self, t_day):
    #    return self.convert_gamma_to_nu(self.gamma_m(t_day), t_day).to(u.Hz)


    #def f_nu_max(self, t_day):
    #    return u.M_e.cgs * cds.c.cgs**2 * constants.sigma_T.cgs * self.ndensity * self.radius(t_day)**3 * self.mag(t_day) * self.gamma2(t_day) / 9. / constants.e.esu / self.dL**2


    def residual_energy(self, t_day):
        print 'Calculating residual energy...'
        print 'Initial energy', self.energy.to(u.erg)
        print 'Sweep length:', self.length().to(u.cm)
        print 'Radius at {0}: {1}'.format(t_day, self.radius(t_day).to(u.cm))
        print 'Lorentz factor at {0}: {1}'.format(t_day, self.gamma2(t_day).decompose())
        return self.energy * (4.*cds.c.cgs*t_day.to(u.s)/self.length())**(-3./7.) / self.gamma0 * sqrt(2.)


    def residual_energy_fluid(self, t_day):
        return 16. * pi * self.gamma2(t_day)**2 * self.radius(t_day)**3 * self.ndensity * u.M_p.cgs * cds.c.cgs**2


    def get_time_transition(self):
        t_f2s = 4.6 * self.epsilonB**(7./5.) * (self.epsilonE_bar*3.)**(7./5.) * (self.energy.to(u.erg).value/1e52)**(4./5.) * (self.gamma0/100.)**(-4./5.) * (self.ndensity.to(u.cm**-3).value)**(3./5.) * 86400.
        if self.iceffect is False:
            return t_f2s
        else:
            t_f2sIC = (self.epsilonE/self.epsilonB) * t_f2s #(self.epsilonE/self.epsilonB)**2 * t_f2s
            t_sIC2sSync = (self.epsilonE/self.epsilonB)**((4.-self.p)/(self.p-2.)) * t_f2sIC
            return  np.array((t_f2sIC, t_sIC2sSync))


#    def get_time_transition_redshift(self):
#        return self.get_time_transition() * (1.+self.z)


    def get_cooling_phase_redshift(self, t_sec):
        if t_sec<self.time_transition_redshift_F2S:
            return 'Fast'
        else:
            logger.error('Slow cooling!! No longer radiative evolution!!')
            return 'Slow'


class ShockedFluid_Mitsunari(ShockedFluid):
    def nu_c(self, t_day):
        return ( 3.*sqrt(3.) * u.M_e.cgs * constants.e.esu * (cds.c.cgs)**0.5 ) / ( 4. * constants.sigma_T.cgs**2 * (self.energy)**0.5 * self.ndensity * u.M_p.cgs * self.epsilonB**(3./2.) * (t_day.to(u.s))**0.5 )

    
    def nu_m(self, t_day):
        return ( 3. * self.energy * self.epsilonB / cds.c**5 / (t_day.to(u.s))**3. )**0.5 * self.epsilonE**2 * ((self.p-2.)/(self.p-1.))**2 * (u.M_p.cgs/u.M_e.cgs)**2 * constants.e.esu / u.M_e.cgs / 8. / pi


    def f_nu_max(self, t_day):
        return 4./3. * (self.ndensity*self.epsilonB/pi/u.M_p.cgs)**0.5 * self.energy * u.M_e.cgs * cds.c * constants.sigma_T.cgs / constants.e.esu / self.dL**2
        

@click.command()
@click.argument('tday', type=float)
@click.argument('nu', type=float)
@click.option('--e52', type=float, default=1)
@click.option('--n0', type=float, default=1)
@click.option('--pelec', '-p', type=float, default=2.5)
@click.option('--epsilonb', '-b', type=float, default=1)
@click.option('--epsilone', '-e', type=float, default=1)
@click.option('--dl28', '-d', type=float, default=None)
@click.option('--redshift', '-z', type=float, default=None)
#@click.option('--values', type=(str, int))
#@click.option('--values', multiple=True)
@click.option('--formula', '-f', type=click.Choice(['Sari', 'Mitsunari']), default='Sari')
@click.option('--evolution', type=click.Choice(['adiabatic', 'radiative']), default='adiabatic')
@click.option('--gamma0', type=float, default=100.)
@click.option('--tunit', type=click.Choice(['s', 'm', 'h', 'd']), default='d')
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(tday, nu, e52, n0, pelec, epsilonb, epsilone, dl28, redshift, formula, evolution, gamma0, tunit, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    DICT_TUNIT = {'s': u.s,
                  'm': u.min,
                  'h': u.h,
                  'd': u.d}
    tday = tday * DICT_TUNIT[tunit]

    if formula == 'Sari':
        if evolution == 'adiabatic':
            shocked_fluid = ShockedFluid_Sari(e52, n0, pelec, epsilonb, epsilone, dl28, redshift)
        elif evolution == 'radiative':
            shocked_fluid = ShockedFluid_Sari_radiative(e52, n0, pelec, epsilonb, epsilone, dl28, redshift, gamma0)
    elif formula == 'Mitsunari':
        shocked_fluid = ShockedFluid_Mitsunari(e52, n0, pelec, epsilonb, epsilone, dl28, redshift)

    shocked_fluid.info()
    logger.info('v_c = {nuc:1.2E} at T = {td}'.format(nuc=shocked_fluid.nu_c_redshift(tday).to(u.Hz), td=tday))
    logger.info('hv_c = {hnuc:1.2E} at T = {td}'.format(hnuc=shocked_fluid.hnu_c_redshift(tday).to(u.eV), td=tday))
    logger.info('v_m = {num:1.2E} at T = {td}'.format(num=shocked_fluid.nu_m_redshift(tday).to(u.Hz), td=tday))
    logger.info('hv_m = {hnum:1.2E} at T = {td}'.format(hnum=shocked_fluid.hnu_m_redshift(tday).to(u.eV), td=tday))
    logger.info('F_max = {fmax:1.2E} at T = {td}'.format(fmax=shocked_fluid.f_nu_max_redshift(tday).to(u.uJy), td=tday))
    logger.info('F = {fval:1.2E} at T = {td}, nu = {nu:1.2E}'.format(fval=shocked_fluid.f_nu_redshift(tday, nu*u.Hz).to(u.uJy), td=tday, nu=nu))

if __name__ == '__main__':
    main()
