import sys
import os


def InputFile(fileName):
	with open(fileName) as f:
		data = f.readlines()
	return data


def OutputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def GetOGPairs(OGStatFile):
	genesOGPair = {}
	for line in InputFile(OGStatFile):
		OGName = line.strip().split('\t')[0]
		genesList = [gene.split('|')[1] for gene in line.strip().split('\t')[-1].split(' ')]
		for gene in genesList:
			genesOGPair[gene] = OGName
	return genesOGPair


def StatisicOG(FNFiles, OGPair, labelType, output):
	FN_stat = ''
	for fn in InputFile(FNFiles):
		gene1 = fn.strip().split('\t')[0]
		gene2 = fn.strip().split('\t')[1]
		label = fn.strip().split('\t')[2]
		if label == labelType:
			if gene1 in OGPair.keys() and gene2 in OGPair.keys():
				OG1 = OGPair[gene1]
				OG2 = OGPair[gene2]
				if OG1 == OG2:
					FN_stat += '%s\t%s\t%s\tSameOG\n' % (gene1, gene2, labelType)
				else:
					FN_stat += '%s\t%s\t%s\tDiffOG\n' % (gene1, gene2, labelType)
					# print(gene1, OG1, gene2, OG2)
			elif gene1 not in OGPair.keys() and gene2 in OGPair.keys():
				print('%s Genes Not Found' % gene1)
			elif gene1 in OGPair.keys() and gene2 not in OGPair.keys():
				print('%s Genes Not Found' % gene2)
			else:
				print('%s %s Genes Not Found' % (gene1, gene2))
	OutputFile(output, FN_stat)


if __name__ == '__main__':
	OGStat = sys.argv[1:][0]
	FNFile = sys.argv[1:][1]
	labelType = sys.argv[1:][2]
	OutPut = sys.argv[1:][3]
	OGP = GetOGPairs(OGStat)
	StatisicOG(FNFile, OGP, labelType, OutPut)