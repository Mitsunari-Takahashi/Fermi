#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import sqrt, cos, sin, tan, acos, asin, atan, radians, degrees, pi
from astropy import units as u
from astropy.units import cds
import astropy.constants as constants
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import Distance
#from astropy.cosmology import WMAP5, WMAP7
from astropy.cosmology import WMAP9 as cosmo
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
    def __init__(self, e52, n0, p, epsilonB, epsilonE, dL28=None, z=None, gamma0=100.):
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


class ShockedFluid_Sari(ShockedFluid):
    def nu_c(self, t_day):
        return 9./4.*sqrt(2./17.) * (1.+self.z)**-0.5 * constants.e.esu * u.M_e.cgs/u.M_p.cgs * cds.c**0.5 / constants.sigma_T.cgs**2 * self.epsilonB**-1.5 * self.ndensity**-1 * (t_day.to(u.s))**-0.5 * self.energy**-0.5

    
    def nu_m(self, t_day):
        return ( 17. / 128. * self.energy * self.epsilonB / cds.c**5 / (t_day.to(u.s))**3. )**0.5 * (1.+self.z)**0.5 * self.epsilonE**2 * ((self.p-2.)/(self.p-1.))**2 * (u.M_p.cgs/u.M_e.cgs)**2 * constants.e.esu / u.M_e.cgs / pi


    def f_nu_max(self, t_day):
        return 17./18. * (1.+self.z) * (self.ndensity*self.epsilonB/2./pi/u.M_p.cgs)**0.5 * self.energy * u.M_e.cgs * cds.c * constants.sigma_T.cgs / constants.e.esu / self.dL**2


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
@click.option('--tunit', type=click.Choice(['s', 'm', 'h', 'd']), default='d')
#@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(tday, nu, e52, n0, pelec, epsilonb, epsilone, dl28, redshift, formula, tunit, loglevel):
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
        shocked_fluid = ShockedFluid_Sari(e52, n0, pelec, epsilonb, epsilone, dl28, redshift)
    elif formula == 'Mitsunari':
        shocked_fluid = ShockedFluid_Mitsunari(e52, n0, pelec, epsilonb, epsilone, dl28, redshift)

    shocked_fluid.info()
    logger.info('v_c = {nuc:1.2E} at T = {td}'.format(nuc=shocked_fluid.nu_c(tday).to(u.Hz), td=tday))
    logger.info('hv_c = {hnuc:1.2E} at T = {td}'.format(hnuc=shocked_fluid.hnu_c(tday).to(u.eV), td=tday))
    logger.info('v_m = {num:1.2E} at T = {td}'.format(num=shocked_fluid.nu_m(tday).to(u.Hz), td=tday))
    logger.info('hv_m = {hnum:1.2E} at T = {td}'.format(hnum=shocked_fluid.hnu_m(tday).to(u.eV), td=tday))
    logger.info('F_max = {fmax:1.2E} at T = {td}'.format(fmax=shocked_fluid.f_nu_max(tday).to(u.uJy), td=tday))
    logger.info('F = {fval:1.2E} at T = {td}, nu = {nu:1.2E}'.format(fval=shocked_fluid.f_nu(tday, nu*u.Hz).to(u.uJy), td=tday, nu=nu))

if __name__ == '__main__':
    main()
