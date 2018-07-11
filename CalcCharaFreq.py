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
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., iceffect=False):
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
        self.gamma0 = gamma0
        self.iceffect = iceffect


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
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., iceffect=False):
        ShockedFluid.__init__(self, e52, n0, p, epsilonB, epsilonE, dL28, z, gamma0, iceffect)
        self.time_transition = self.get_time_transition()
        self.time_transition_redshift = self.get_time_transition_redshift()


    # Common for adiabatic and radiative
    def convert_gamma_to_nu(self, gamma, t_day):
        return self.gamma2(t_day) * gamma**2 * constants.e.esu * self.mag(t_day) / (2.*pi*u.M_e.cgs*cds.c)


    def mag(self, t_day):
        return (32. * pi * u.M_p.cgs * self.epsilonB * self.ndensity)**0.5 * self.gamma2(t_day) * cds.c


    def gamma_c(self, t_day):
        return 3.*u.M_e.cgs / (16.*self.epsilonB*constants.sigma_T.cgs*u.M_p.cgs*cds.c * self.gamma2(t_day)**3 * self.ndensity * t_day.to(u.s))


    def gamma_m(self, t_day):
        return self.epsilonE_bar * u.M_p.cgs/u.M_e.cgs * self.gamma2(t_day)


    def diff_nu_breaks(self, t_sec):
        """nu_m - nu_c"""
        return log10(self.nu_m(t_sec*u.s).to(u.Hz).value)-log10(self.nu_c(t_sec*u.s).to(u.Hz).value)


    def diff_nu_redshift_breaks(self, t_sec):
        """nu_m - nu_c"""
        return log10(self.nu_m_redshift(t_sec*u.s).to(u.Hz).value)-log10(self.nu_c_redshift(t_sec*u.s).to(u.Hz).value)


    def get_time_transition(self):
        """Transition time from fast-cooling to slow-cooling in sec in the GRB frame"""
        return 210. * self.epsilonB**2 * (self.epsilonE_bar*3.)**2 * (self.energy.to(u.erg).value/1e52) * (self.ndensity.to(u.cm**-3).value) * 86400.


    def get_time_transition_redshift(self):
        return self.get_time_transition() * (1.+self.z)


    def radiative_energy_fraction(self, t_day):
        ratio = self.gamma_c(t_day) / self.gamma_m(t_day)
        if ratio < 1: # Fast-cooling
            return 1.
        else: # Slow-cooling
            return ratio**(-self.p+2.)


    def luminosity_ratio(self, t_day):
        eta = self.radiative_energy_fraction(t_day)
        return (-1.+(1.+4.*eta*self.epsilonE/self.epsilonB)**0.5) / 2.


    # Adiabatic evolution
    def gamma2(self, t_day):
        return (self.radius(t_day) / 4. / cds.c.cgs / t_day)**0.5 #( 17. * self.energy / (16. * pi * self.ndensity * u.M_p.cgs * cds.c.cgs**2 * self.radius(t_day)**3) )**0.5


    def radius(self, t_day):
        return (17. * self.energy * t_day.to(u.s) / 4. / pi / u.M_p.cgs / self.ndensity / cds.c.cgs)**0.25


    def nu_c(self, t_day):
        return self.convert_gamma_to_nu(self.gamma_c(t_day), t_day).to(u.Hz)


    def nu_m(self, t_day):
        return self.convert_gamma_to_nu(self.gamma_m(t_day), t_day).to(u.Hz)


    def f_nu_max(self, t_day):
        return u.M_e.cgs * cds.c.cgs**2 * constants.sigma_T.cgs * self.ndensity * self.radius(t_day)**3 * self.mag(t_day) * self.gamma2(t_day) / 9. / constants.e.esu / self.dL**2


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
                return (nu/self.nu_c(t_day))**(1./3.) * self.f_nu_max(t_day)
            elif nu<self.nu_m(t_day):
                return (nu/self.nu_c(t_day))**(-1./2.) * self.f_nu_max(t_day)
            else:
                return (nu/self.nu_m(t_day))**(-self.p/2.) * (self.nu_m(t_day)/self.nu_c(t_day))**(-1./2.) * self.f_nu_max(t_day)
        else:
            # Slow cooling
            if nu<self.nu_m(t_day):
                return (nu/self.nu_m(t_day))**(1./3.) * self.f_nu_max(t_day)
            elif nu<self.nu_c(t_day):
                return (nu/self.nu_m(t_day))**(-(self.p-1.)/2.) * self.f_nu_max(t_day)
            else:
                return (nu/self.nu_c(t_day))**(-self.p/2.) * (self.nu_c(t_day)/self.nu_m(t_day))**(-(self.p-1.)/2.) * self.f_nu_max(t_day)


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
                return (nu/self.nu_c_redshift(t_day))**(1./3.) * self.f_nu_max_redshift(t_day)
            elif nu<self.nu_m(t_day):
                return (nu/self.nu_c_redshift(t_day))**(-1./2.) * self.f_nu_max_redshift(t_day)
            else:
                return (nu/self.nu_m_redshift(t_day))**(-self.p/2.) * (self.nu_m_redshift(t_day)/self.nu_c_redshift(t_day))**(-1./2.) * self.f_nu_max_redshift(t_day)
        else:
            # Slow cooling
            if nu<self.nu_m_redshift(t_day):
                return (nu/self.nu_m_redshift(t_day))**(1./3.) * self.f_nu_max_redshift(t_day)
            elif nu<self.nu_c(t_day):
                return (nu/self.nu_m_redshift(t_day))**(-(self.p-1.)/2.) * self.f_nu_max_redshift(t_day)
            else:
                return (nu/self.nu_c_redshift(t_day))**(-self.p/2.) * (self.nu_c_redshift(t_day)/self.nu_m_redshift(t_day))**(-(self.p-1.)/2.) * self.f_nu_max_redshift(t_day)


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


class ShockedFluid_Sari_afterRadiative(ShockedFluid_Sari):
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100., r_transition = 0, t_transition = 0):
        ShockedFluid_Sari.__init__(self, e52, n0, p, epsilonB, epsilonE, dL28, z, gamma0)
        self.r_transition = r_transition
        self.t_transition = t_transition


    def radius(self, t_day):
        def func(r):
            return r**4 - (self.r_transition.to(u.cm).value)**3 * r - 17.*self.energy.to(u.erg).value*(t_day.to(u.s).value)/(4.*pi*(self.ndensity.to(u.cm**-3).value)*u.M_p.cgs.to(u.g)*(cds.c.to(u.cm/u.s)))

        solutions = optimize.fsolve(func, self.r_transition.to(u.cm).value * ((t_day.to(u.s)/self.t_transition.to(u.s)))**0.25)#, xtol=1e-3)
        #print solutions
        return solutions[0]*u.cm



class ShockedFluid_Sari_radiative(ShockedFluid_Sari):

    def length(self):
        return ((17.*self.energy/self.gamma0/(cds.c.cgs)**2 / 16./pi/u.M_p.cgs/self.ndensity).to(u.cm**3))**(1./3.)


    def radius(self, t_day):
        return (4.*cds.c.cgs*t_day.to(u.s)/self.length())**(1./7.) * self.length()


    #def gamma2(self, t_day):
    #    return (4.*cds.c.cgs*t_day.to(u.s)/self.length())**(-3./7.)


    def nu_c(self, t_day):
        return self.convert_gamma_to_nu(self.gamma_c(t_day), t_day).to(u.Hz)


    def nu_m(self, t_day):
        return self.convert_gamma_to_nu(self.gamma_m(t_day), t_day).to(u.Hz)


    def f_nu_max(self, t_day):
        return u.M_e.cgs * cds.c.cgs**2 * constants.sigma_T.cgs * self.ndensity * self.radius(t_day)**3 * self.mag(t_day) * self.gamma2(t_day) / 9. / constants.e.esu / self.dL**2


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
        return 4.6 * self.epsilonB**(7./5.) * (self.epsilonE_bar*3.)**(7./5.) * (self.energy.to(u.erg).value/1e52)**(4./5.) * (self.gamma0/100.)**(-4./5.) * (self.ndensity.to(u.cm**-3).value)**(3./5.) * 86400.


    def get_time_transition_redshift(self):
        return self.get_time_transition() * (1.+self.z)


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
#@click.option('--shout', is_flag=True)
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
