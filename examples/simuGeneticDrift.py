#!/usr/bin/env python
#
# Demonstrate changes of allele frequency due to genetic drift. 

"""
This program demonstrates changes of allele frequency on single locus due to genetic drift.
"""

import simuOpt, os, sys, time
from simuPOP import *

try:
    from simuPOP.plotter import VarPlotter
except:
    print("simuRPy import failed. Please check your rpy installation.")
    print("Allele Frequencies will not be plotted")
    useRPy = False
else:
    useRPy = True

def simuGeneticDrift(popSize=100, p=0.2, generations=100, replications=5):
    '''Simulate the Genetic Drift as a result of random mating.'''
    # diploid population, one chromosome with 1 locus
    # random mating with sex
    pop = Population(size=popSize, loci=[1])
    simu=Simulator(pop, rep=replications)

    if useRPy:
        plotter = VarPlotter('alleleFreq[0][0]', ylim=[0, 1], ylab='allele frequency',
            update=generations-1, saveAs='geneticDrift.png')
    else:
        plotter = NoneOp()

    # if number of generation is smaller than 200, step is 10 generations,
    # if it's between 200 and 500, set step to be 20 generations,
    # otherwise, step = 50 generations.
    if generations <= 200:
        s = 10
    elif 200 < generations <= 500:
        s = 20
    else:
        s = 50
        
    simu.evolve(
        # everyone initially will have the same allele frequency
        initOps = [
            InitSex(),
            InitGenotype(freq=[p, 1-p])
        ],
        matingScheme = RandomMating(),
        postOps = [
            Stat(alleleFreq=[0]),
            PyEval(r'"Generation %d:\t" % gen', reps = 0, step = s),
	        PyEval(r"'%.3f\t' % alleleFreq[0][0]", step = s),
	        PyOutput('\n', reps=-1, step = s),
	        plotter,
        ],
        gen = generations
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
        default=100,
        type=check_positive_int,
        help="Population Size")
    args.add_argument("--p",
        default=0.2,
        type=check_probability,
        help="Initial Allele Frequency")
    args.add_argument("--generations",
        default=100,
        type=check_positive_int,
        help="Number of Generations")
    args.add_argument("--replications",
        default=5,
        type=check_positive_int,
        help="Number of Replicates")
    args = args.parse_args()

    simuGeneticDrift(**vars(args))

    # wait ten seconds before exit
    if useRPy:
        print("Figure will be closed after 5 seconds.")
        time.sleep(5)
