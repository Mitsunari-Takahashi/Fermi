#!/usr/bin/env python
import sys
import click
import pandas as pd
import itertools
from scipy.stats import poisson


DCT_TINTERVAL = {'T1':'T0-100s', 'T2':'T100-1000s', 'T3':'T1000-10000s'}

LST_NCOINCI = [{'T1':0, 'T2':1, 'T3':3}]
# LST_NCOINCI = [{'T1':0, 'T2':0, 'T3':0},
#                {'T1':0, 'T2':0, 'T3':1},
#                {'T1':0, 'T2':0, 'T3':2},
#                {'T1':0, 'T2':0, 'T3':3},
#                {'T1':0, 'T2':1, 'T3':0},
#                {'T1':0, 'T2':1, 'T3':1},
#                {'T1':0, 'T2':1, 'T3':2}]


def Read_cvs(path_cvs, nskip=None):
    df = pd.read_csv(path_cvs, index_col='GRB', skipfooter=nskip)
    return df


def CalcMultivariateMultinomial(path_cvs, skip, noreplace):
    DF_PROB = Read_cvs(path_cvs, skip)
    print DF_PROB
    NGRB = len(DF_PROB.index)
    dct_combis = {}
    prob_sum = 0.
    combi_time = {}

    for dct_coinci in LST_NCOINCI:
        print '===================='
        print dct_coinci
        for keyt in DCT_TINTERVAL.keys(): #preparation for next loop
            if noreplace==False:
                dct_combis[keyt] = list(itertools.combinations_with_replacement(DF_PROB.index, dct_coinci[keyt]))
            else:
                dct_combis[keyt] = list(itertools.combinations(DF_PROB.index, dct_coinci[keyt]))
            print keyt, dct_combis[keyt]
        #for combi_time['T1'], combi_time['T2'], combi_time['T3'] in itertools.product(dct_combis.values()):
        for combi_time['T1'] in dct_combis['T1']:
            for combi_time['T2'] in dct_combis['T2']:
                for combi_time['T3'] in dct_combis['T3']:
                    print combi_time['T1'], combi_time['T2'], combi_time['T3']
                    prob = 1.
                    for grb in DF_PROB.index:
                        for keytime, strtime in DCT_TINTERVAL.items():
                            prob = prob * poisson.pmf(combi_time[keytime].count(grb), DF_PROB.ix[grb][strtime])
                    print prob
                    sys.stdout.flush()
                    prob_sum += prob
    print 'Total probability:', prob_sum
    return prob_sum

                        
@click.command()
@click.argument('cvs', type=str)
#@click.argument('probs', nargs=-1, type=float)
@click.option('--skip', default=None, type=int)
@click.option('--noreplace', is_flag=True)
def main(cvs, skip, noreplace):
    CalcMultivariateMultinomial(cvs, skip, noreplace)


if __name__ == '__main__':
    main()
