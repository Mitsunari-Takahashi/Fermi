#!/usr/bin/env python

import sys
import os
import os.path
import datetime
import click


def proc(name, model, time, suffix='', nfit=1000):
    name_plot = 'FIG_GRB{0}_{1}_{2}{3}'.format(name, time, model, suffix)
    str_proc = """fit {0}
cpd /xs
pl eeufspec ratio
cpd {1}.ps/cps color
pl eeufspec ratio
cpd {1}.gif/gif
pl eeufspec ratio
mv {1}.gif_2 {1}.gif
""".format(nfit, name_plot)
    return str_proc


def ignore(detector):
    if detector=='LAT':
        return '**-1e5'
    elif detector=='LLE':
        return '**-3e4'
    elif detector[:3]=='NAI':
        return '**-8.0 1e3-**'
    elif detector[:3]=='BGO':
        return '**-150.0 3e4-**'
    else:
        return ''

def rebin(detector):
    if detector=='LAT':
        return '2 20'
    elif detector=='LLE':
        return '2 20'
    elif detector[:3]=='NAI':
        return '5 20'
    elif detector[:3]=='BGO':
        return '5 20'
    else:
        return ''


@click.command()
@click.option('--name', type=str, default='160509374')
@click.option('--suffix', '-s', type=str, default='LATcstat')
def main(name, suffix):
    if suffix!='':
        suffix = '_' + suffix

    LST_DATA = ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'BGO_01', 'LLE']

    LST_TIME_EDGE = [0.0, 10.0, 14.0, 16.0, 18.0, 20.0, 30.0, 64.0, 82.0, 285.0, 358.0, 385.0]

    NDIV_TIME = len(LST_TIME_EDGE)-1

    for idiv in range(NDIV_TIME):
        str_time = '{0:0>5d}-{1:0>5d}sec'.format(int(LST_TIME_EDGE[idiv]+0.5), int(LST_TIME_EDGE[idiv+1]+0.5))
        print '=====', str_time, '====='
        str_data = """log Fitting_GRB{0}_{1}.log
""".format(name, str_time)
        for (idata, dataset) in enumerate(LST_DATA):
            if dataset=='LAT':
                str_data = str_data + """data {0:d}:{1:d} {2}_{3}_ROI_E.pha2{{4}} #{3}
""".format(idiv+1, idata+1, name, dataset, idiv+1)
            else:
                str_data = str_data + """data {0:d}:{1:d} GRB{2}_{3}.pha{{4}} #{3}
""".format(idiv+1, idata+1, name, dataset, idiv+1)
        for (idata, dataset) in enumerate(LST_DATA):
            str_data = str_data + """ignore {0:d}:{1} #{2}
""".format(idata+1, ignore(dataset), dataset)
        str_data = str_data + """ignore bad
setplot energy
"""
        for (idata, dataset) in enumerate(LST_DATA):
            str_data = str_data + """setplot rebin {0} {1:d} #{2}
""".format(rebin(dataset), idata+1, dataset)

        str_data = str_data + """weight churazov
"""
        for (idata, dataset) in enumerate(LST_DATA):
            if dataset=='LAT' or dataset=='LLE':
                str_data = str_data + """statistic cstat {0:d} #{1}
""".format(idata+1, dataset)

        # MODEL
        # Black Body
        str_data = str_data + """mod bbody
500 1 10 100 1000 2000
1E5
{0}
""".format(proc(name, 'bbody', str_time, suffix))

        # Band
        str_data = str_data + """mod grbm
-1
-2
300
0.1
{0}
""".format(proc(name, 'grbm', str_time, suffix))

        # ExpBand
        str_data = str_data + """editmod grbm*highect
0 0 0 0 0 0
1E5 1E3 1E3 1E4 1E6 2E6
{0}
""".format(proc(name, 'grbm-highect', str_time, suffix))

        #CLOSING
        str_data = str_data + """cpd none
data none
log none
"""
    
        with open("cmd_fit_{0}{1}.xcm".format(str_time, suffix), 'w') as xcm:
            xcm.write(str_data)


if __name__ == '__main__':
    main()
