#!/usr/bin/env python

from astropy.io import fits
import click


def open_table(nextention=1, path_catalogue='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'):
    """Open your FITS file and return the extention table you assigned.
"""
    f = fits.open(path_catalogue)
    return f[nextention].data


def select_one_by_name(tb, grbname):
    """Return a masked table consists of one paticular GRB.
"""
    tb_masked = tb[tb['GRBNAME'] == grbname]
    if len(tb_masked)==1:
        return tb_masked[0]
    elif len(tb_masked)<1:
        print 'No data of', grbname
        return 1
    elif len(tb_masked)>1:
        print 'No unique candidate of', grbname
        return 1


def select_by_name(tb, name_min='0', name_max='200000000'):
    """Return a masked table consists of GRBs from name_min to name_max. The names should be string.
"""
    return tb[(tb['GRBNAME'] >= name_min) * (tb['GRBNAME'] < name_max)]


def select_gbm_exist(tb):
    """Return a masked table for GRBs whose GBM data exist.
"""
    return tb[tb['TRIGGER_TIME'] == tb['TRIGGER_TIME']]


def select_by_mjd(tb, mjd_min=0., mjd_max=58000.):
    tb_gbm = select_gbm_exist(tb)
    return tb_gbm[(tb_gbm['TRIGGER_TIME']>=mjd_min) * (tb_gbm['TRIGGER_TIME']<mjd_max)]


def select_by_fluence(tb, flu_min=0., flu_max=1.,):
    tb_gbm = select_gbm_exist(tb)
    return tb_gbm[(tb_gbm['TRIGGER_TIME']>=mjd_min) * (tb_gbm['TRIGGER_TIME']<mjd_max)]    


def select_long(tb):
    tb_gbm = select_gbm_exist(tb)
    return tb_gbm[tb_gbm['T90']>=2.0]    


def select_short(tb):
    tb_gbm = select_gbm_exist(tb)
    return tb_gbm[tb_gbm['T90']<2.0] 


def select_redshift_known(tb):
    """Return a masked table for GRBs whose redshift are known.
"""
    return tb[tb['REDSHIFT'] > 0]


def select_by_redshift(tb, z_min=0., z_max=10.,):
    return tb[(tb['REDSHIFT']>z_min) * (tb['REDSHIFT']<=z_max)]    


@click.command()
@click.argument('grb', type=str)
@click.option('--table', '-t',default="/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits", help="Path of FITS table.")
@click.option('--items', '-i', multiple=True, default=None)
def main(grb, table, items):
    tb = open_table(1, table)
    tb_masded = select_one_by_name(tb, grb)
    if items is None:
        print tb_masded
    elif len(items)>0:
        for item in items:
            print item, tb_masded[0][item]
    else:
        print tb_masded

if __name__ == '__main__':
    main()
