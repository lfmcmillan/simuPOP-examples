#!/usr/bin/env python
#
# Demonstrate the decay of linkage disequilibrium
#
# Author: Bo Peng (bpeng@mdanderson.org)
#

"""
This program demonstrate the decay of linkage disequilibrium due
to recombination.
"""

import simuOpt, os, sys, types, time
from simuPOP import *

try:
    from simuPOP.plotter import *
except:
    print("simuRPy import failed. Please check your rpy installation.")
    print("LD values will not be plotted")
    useRPy = False
else:
    useRPy = True


def simuLDDecay(popSize, gen, recRate, numRep, measure, saveFigure, useRPy):
    '''Simulate the decay of linkage disequilibrium as a result
    of recombination.
    '''
    # diploid population, one chromosome with 2 loci
    # random mating with sex
    simu = Simulator(
        Population(size=popSize, ploidy=2, loci=[2]),
        rep = numRep)

    # get method value used to plot and evolve
    if measure=="D'":
        methodplot = "LD_prime[0][1]"
        upperlim = 1
        methodeval = r"'%.4f\t' % LD_prime[0][1]"
    elif measure=='R2':
        methodplot = "R2[0][1]"
        upperlim = 1
        methodeval = r"'%.4f\t' % R2[0][1]"
    else:
        methodplot = "LD[0][1]"
        upperlim = 0.25
        methodeval = r"'%.4f\t' % LD[0][1]"

    if useRPy:
        print(saveFigure)
        plotter = VarPlotter(methodplot, 
            ylim = [0, upperlim], saveAs=saveFigure,
            update = gen - 1, ylab=method, leaveOpen=True,
            main="Decay of Linkage Disequilibrium r=%f" % recRate)
    else:
        plotter = NoneOp()

    simu.evolve(
        # everyone will have the same genotype: 01/10
        initOps = [
            InitSex(),
            InitGenotype(genotype=[0,1,1,0])
        ],
        matingScheme = RandomMating(ops=Recombinator(rates = recRate)), 
        postOps = [
            Stat(alleleFreq=[0], LD=[0, 1]),
            PyEval(methodeval),
            PyOutput('\n', reps=-1),
            plotter
        ],
        gen = gen
    )

def check_positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

def check_probability(value):
    fvalue = float(value)
    if fvalue <= 0 or fvalue >= 1:
        raise argparse.ArgumentTypeError("%s is an invalid frequency value" % value)
    return fvalue

if __name__ == '__main__':
    import argparse
    args = argparse.ArgumentParser(description="This program demonstrates changes of allele frequency on single locus due to genetic drift.")

    args.add_argument("--popSize",
        default=1000,
        type=check_positive_int,
        help="Population Size") 
    args.add_argument("--gen",
        default=50,
        type=check_positive_int,
        help="Generations to Evolve")
    args.add_argument("--recRate",
        default=0.01,
        type=check_probability,
        help="Recombination Rate")
    args.add_argument("--numRep",
        default=5,
        type=check_positive_int,
        help="Number of Replicates")
    args.add_argument("--measure",
        default='D',
        choices=['D', "D'", 'R2'],
        help="LD measure to be displayed.")
    args.add_argument("--saveFigure",
        default='',
        type=str,
        help="Filename template to save figures to. They will be saved as e.g. 'filename_10.eps' The format of the figures is determined by the file extension.")
    args = args.parse_args()

    simuLDDecay(**vars(args), useRPy = useRPy)

    # wait five seconds before exit
    if useRPy:
        print("Figure will be closed after five seconds.")
        time.sleep(5)
