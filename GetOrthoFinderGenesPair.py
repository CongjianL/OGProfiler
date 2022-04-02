import sys
import os
from itertools import combinations as comb


def inputFile(fileName):
	with open(fileName) as f:
		data = f.readlines()
	return data


def outputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def GetGenesPair(OGList):
	pairwise = comb(OGList, 2)
	print(pairwise)
	pairwiseTab = '\n'.join(['%s\t%s' % (A, B) for A, B in pairwise])
	return pairwiseTab


def GetSingleCopyOG(OGFile, SingleOGs):
	SOGs =[g.strip('\n') for g in inputFile(SingleOGs)]
	OGs = inputFile(OGFile)
	SOGPairwise = ''
	for ogLine in OGs:
		OGName = ogLine.split(':')[0]
		if OGName in SOGs:
			genes = ogLine.split(':')[1].strip(' ').split(' ')
			genesID = [gene.split('|')[1] for gene in genes]
			ogPairs = GetGenesPair(genesID)
			SOGPairwise += ogPairs + '\n'
	return SOGPairwise


if __name__ == '__main__':
	SOGFile = sys.argv[1:][0]
	OGFileName = sys.argv[1:][1]
	OutputName = sys.argv[1:][2]
	pairs = GetSingleCopyOG(OGFileName, SOGFile)
	outputFile(OutputName, pairs)