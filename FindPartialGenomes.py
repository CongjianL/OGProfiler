import pickle as pic
import sys
import os
from collections import Counter

def MatricesLoad(MatriceName):
    with open(MatriceName, 'rb') as picFile:
        M = pic.load(picFile)
    return M


def OutputFile(fileName, content):
    with open(fileName, 'w') as f:
        f.writelines(content)


def ReadFile(FileName):
    with open(FileName) as file:
        data = file.readlines()
    return data


def StatGenomes(PredictionFile, SeqInfFile):
    SeqInf = MatricesLoad(SeqInfFile)
    GenomeStat = []
    for line in ReadFile(PredictionFile):
        gene1, gene2, Label = line.strip().split('\t')[:3]
        genome1, genome2 = gene1.split('|')[0], gene2.split('|')[0]
        RecodeGenome1, RecodeGenome2 = SeqInf['GenomeRecodeInf'][genome1], SeqInf['GenomeRecodeInf'][genome2]
        GenomeStat.extend([RecodeGenome1, RecodeGenome2])
    GenomeNumStat = sorted([(genome, num) for genome, num in Counter(GenomeStat).items()], key=lambda X:X[1], reverse=True)[:10]
    [print(Genome) for Genome, Num in GenomeNumStat]

if __name__ == '__main__':
    PFile = sys.argv[1:][0]
    SeqInf = sys.argv[1:][1]
    StatGenomes(PFile, SeqInf)