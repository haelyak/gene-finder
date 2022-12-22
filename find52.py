import math
from findHelper import *

def find52():
    """Given coding and noncoding training data from E. coli, build a
    first order Markov model on codons and make predictions for
    sequences in Vibrio.
    """
    # load sequences
    ecoliCodeOrfL=[s.rstrip() for s in open("ecoliCodeTrainOrfs-mid.txt","r").readlines()]
    ecoliNoncodeOrfL=[s.rstrip() for s in open("ecoliNcTrainOrfs-mid.txt","r").readlines()]
    vibAllOrfL=[s.rstrip() for s in open("vibAllOrfs-mid.txt","r").readlines()]
    # run model
    ecoliCodeProbD=condProb(ecoliCodeOrfL)
    ecoliNoncodeProbD=condProb(ecoliNoncodeOrfL)
    llrD=makeLogLikelihoodRatioD(ecoliCodeProbD,ecoliNoncodeProbD)
    vibCodePredL,vibNoncodePredL=predict(vibAllOrfL,llrD)

    return vibCodePredL,vibNoncodePredL

def makeLogLikelihoodRatioD(codeProbD,noncodeProbD):
    """Make a dictionary of log likelihood ratios from coding and
    noncoding Markov models."""
    llrD={}
    for key in codeProbD:
        llrD[key] = math.log( codeProbD[key]/noncodeProbD[key] )
    return llrD

def count(orfL):
    """Takes in a list of open reading frames
    counts codons and codon pairs
    returns these counts as two dictionaries"""
    codonCountD, twoCodonCountD = initializeCountDicts()
    for x in range(len(orfL)):
        for y in range(len(orfL[x])-3):
            if y%3 == 0:
                twoCodonCountD[orfL[x][y:y+6]] += 1
                codonCountD[orfL[x][y:y+3]] += 1
    return codonCountD, twoCodonCountD

def condProb(orfL): 
    """Takes in a list of open reading frames
    returns a dictionary of conditional probabilities
    for each codon pair"""
    codonCountD, twoCodonCountD = count(orfL) 
    for key in twoCodonCountD:
        twoCodonCountD[key] = twoCodonCountD[key]/codonCountD[key[:3]]
    return twoCodonCountD

def logLikelihoodRatioSum(orf,llrD):
    """Takes in an open reading frame
    returns a sum of the log likelihood 
    ratio for every codon pair in the open
    reading frame"""
    runningSum = 0 
    for y in range(len(orf)-3):
            if y%3 == 0:
                runningSum += llrD[orf[y:y+6]]
    return runningSum


def predict(orfL,llrD):
    """Takes in a list of open reading frames
    and a log likelihood dictionary
    Returns a list of coding orfs 
    and a list of noncodoning orfs"""
    codingL = []
    noncodingL = []
    for x in range(len(orfL)):
        if logLikelihoodRatioSum(orfL[x], llrD) < 0:
            noncodingL += [orfL[x]]
        else:
            codingL += [orfL[x]]
    return codingL, noncodingL