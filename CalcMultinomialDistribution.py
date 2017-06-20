#!/usr/bin/env python
from math import factorial
import click


def CalcMultinomialDistribution3(ntrial, probs, count):
    xp = probs[0]
    yp = probs[1]
    zp = probs[2]
    rp = 1.-xp-yp-zp
    prob_zero = 0
    prob_total = 0
    print 'Input probabilities:', xp, yp, zp, rp
    for xi in range(ntrial+1):
        for yi in range(ntrial+1):
            for zi in range(ntrial+1):
                ri = ntrial-xi-yi-zi
                if ri>=0:
                    prob_coincide = factorial(ntrial)/factorial(xi)/factorial(yi)/factorial(zi)/factorial(ri)*(xp**xi)*(yp**yi)*(zp**zi)*(rp**ri)
                    print '({0}, {1}, {2}, {3}) : {4}'.format(xi, yi, zi, ri, prob_coincide)
                    prob_total += prob_coincide
                    if count=='any':
                        if xi==0 or yi==0 or zi==0:
                            prob_zero += prob_coincide
                    elif count=='all':
                        if xi==0 and yi==0 and zi==0:
                            prob_zero += prob_coincide

    print 'Total probability of cases whose', count, 'coincidence event numbers are zero:', prob_zero
    print 'Probability of the other cases:', 1.-prob_zero
    print 'Total probability check:', prob_total
    return prob_zero
    
                        
@click.command()
@click.argument('ntrial', nargs=1, type=int)
@click.argument('probs', nargs=-1, type=float)
@click.option('--zero', type=click.Choice(['any', 'all']))
def main(ntrial, probs, zero):
    if len(probs)==3:
        CalcMultinomialDistribution3(ntrial, probs, zero)
    else:
        print 'Please input three probabilities.'


if __name__ == '__main__':
    main()
