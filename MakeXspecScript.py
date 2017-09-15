#!/usr/bin/env python

import sys
import os
import os.path
import datetime
import click


def proc(name, model, time, suffix='', lst_plot=['ldata'], ebounds=None, nfit=1000):
    str_proc = """fit {0}
""".format(nfit)
    if not os.path.exists('./plots'):
        os.mkdir('plots')
    str_ebounds = 'setplot command Rescale X '
    if ebounds!=None and len(ebounds)==2:
        str_ebounds += """{0:e} {1:e}
""".format(ebounds[0], ebounds[1])
    for plot in lst_plot:
        name_plot = 'plots/FIG_GRB{0}_{1}_{2}_{3}{4}'.format(name, plot, time, model, suffix)
        str_proc = str_proc + """cpd /xs
pl {0} ratio
{2}
setplot command Rescale Y2 -3.0 5.0
hardcopy {1}.ps color
setplot device {1}.gif/vgif
pl {0} ratio
{2}
setplot command Rescale Y2 -1.0 3.0
mv {1}.gif_2 {1}.gif
""".format(plot, name_plot, str_ebounds)
    return str_proc


def ignore(detector):
    if detector=='LAT':
        return '**-1e5'
    elif detector=='LLE':
        return '**-3e4 3e5-**'
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
        return '3 1000'
    elif detector[:3]=='NAI':
        return '5 50'
    elif detector[:3]=='BGO':
        return '5 50'
    else:
        return ''


@click.command()
@click.option('--name', type=str, default='160509374')
@click.option('--suffix', '-s', type=str, default='LAT-LLEcstat')
@click.option('--plot', '-p', type=str, default='')
def main(name, suffix, plot):
    if plot!='':
        suffix += '_' + plot
    if suffix!='':
        suffix = '_' + suffix

    lst_plot = ['eeufspec', 'icounts']
    if plot!='':
        lst_plot = [plot]

    LST_TIME_EDGE = [0.0, 10.0, 14.0, 16.0, 18.0, 20.0, 30.0, 64.0, 82.0, 285.0, 358.0, 385.0]
    LST_DATA_MOTHER = [['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'LLE'], #0-10s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'LLE'], #10-14s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'LLE'], # 14-16s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'LLE'], #16-18s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'LLE'], #18-20s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'BGO_00', 'LLE'], #20-30s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'NAI_07', 'NAI_09', 'BGO_00'], #30-64s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'NAI_07', 'NAI_09', 'BGO_00'], #64-82s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'NAI_07', 'NAI_09', 'BGO_00', 'BGO_01'], #82-285s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'NAI_07', 'NAI_09', 'BGO_00', 'BGO_01'], #285-358s
                       ['LAT', 'NAI_00', 'NAI_01', 'NAI_03', 'NAI_06', 'NAI_07', 'NAI_09', 'BGO_00', 'BGO_01']] #358-385s
    #NAI<50 deg, BGO<100 deg
    NDIV_TIME = len(LST_TIME_EDGE)-1

    if not os.path.exists('./logs'):
        os.mkdir('logs')

    for idiv in range(NDIV_TIME):
        str_time = '{0:0>5d}-{1:0>5d}sec'.format(int(LST_TIME_EDGE[idiv]+0.5), int(LST_TIME_EDGE[idiv+1]+0.5))
        print '=====', str_time, '====='
        LST_DATA = LST_DATA_MOTHER[idiv]
        print LST_DATA
        str_data = """log logs/Fitting_GRB{0}_{1}.log
""".format(name, str_time)
        for (idata, dataset) in enumerate(LST_DATA):
            if dataset=='LAT':
                str_data = str_data + """data 1:{0:d} {1}_{2}_ROI_E.pha2{{{3}}}
""".format(idata+1, name, dataset, idiv+1)
            else:
                str_data = str_data + """data 1:{0:d} GRB{1}_{2}.pha{{{3}}}
""".format(idata+1, name, dataset, idiv+1)
        for (idata, dataset) in enumerate(LST_DATA):
            str_data = str_data + """ignore {0:d}:{1}
""".format(idata+1, ignore(dataset), dataset)
        str_data = str_data + """ignore bad
setplot energy
"""
        for (idata, dataset) in enumerate(LST_DATA):
            str_data = str_data + """setplot rebin {0} {1:d}
""".format(rebin(dataset), idata+1, dataset)

        str_data = str_data + """weight churazov
"""
        for (idata, dataset) in enumerate(LST_DATA):
            if dataset=='LAT' or dataset=='LLE':
                str_data = str_data + """statistic cstat {0:d}
""".format(idata+1, dataset)
            else:
                str_data = str_data + """statistic pgstat {0:d}
""".format(idata+1, dataset)

        # MODEL
        #Black Body
        str_data = str_data + """mod bbody
100 1 1 50 1000 2000
1E5
{0}
""".format(proc(name, 'bbody', str_time, suffix, lst_plot, [8.0, 3e4]))

        # Band
        str_data = str_data + """mod grbm
-1
-2
300 10 10 50 5000 10000
0.1
{0}
""".format(proc(name, 'grbm', str_time, suffix, lst_plot))

        # ExpBand
        str_data = str_data + """editmod grbm*highecut
0 0 0 0 0 0
1E5 1E2 1E2 1E3 5E5 2E6
{0}
""".format(proc(name, 'grbmhighect', str_time, suffix, lst_plot))

        # ExpPwl
        str_data = str_data + """editmod powerlaw*highecut
1.5 0.01 -4.0 -3.0 3.0 4.0
10.0
{0}
""".format(proc(name, 'powerlawhighecut', str_time, suffix, lst_plot))


        # ExpBand + Black body
#         str_data = str_data + """editmod grbm*highecut+bbody
# 50 1 10 100 1000 2000
# 25
# {0}
# """.format(proc(name, 'grbmhighecut-cutoffpwl', str_time, suffix))

#         # ExpBand+cutoffpwl
#         str_data = str_data + """editmod grbm*highecut+cutoffpl
# 2
# 1E8 1E5 1E5 1E6 1E9 2E9
# 1E-8
# {0}
# """.format(proc(name, 'grbmhighecut-cutoffpwl', str_time, suffix))

        #CLOSING
        str_data = str_data + """cpd none
data none
log none
"""
    
        with open("cmd_fit_{0}{1}.xcm".format(str_time, suffix), 'w') as xcm:
            xcm.write(str_data)


if __name__ == '__main__':
    main()
