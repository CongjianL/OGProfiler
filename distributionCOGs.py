import re
import sys

import igraph as ig


def inputFile(fileName):
	with open(fileName) as f:
		data = f.read()
	return data


def inputFilelines(fileName):
	with open(fileName) as f:
		data = f.readlines()
	return data


def COGNumCMapping(Cogsfile):
	Mapping = {}
	cogMapping = inputFilelines(Cogsfile)
	for line in cogMapping:
		COGNum = line.strip('\n').split('\t')[0]
		try:
			C = line.strip('\n').split('\t')[1]
			Mapping[COGNum] = C
		except IndexError:
			print(line)
	return Mapping


def MappingCOGForIDs(idMappingFile):
	idMapping = inputFile(idMappingFile)
	p = r'\w+\seggNOG\tCOG\w+'
	COGIterms = re.findall(p, idMapping)
	OGPair = {}
	for idsLine in COGIterms:
		ids = idsLine.strip('\n').split('\t')[0]
		COG = idsLine.strip('\n').split('\t')[2]
		OGPair[ids] = COG
	return OGPair


def GetGenesCOGs(OGCOGsPair, cogCPair, hhnName):
	hhn = ig.read(hhnName)
	OGNodes = hhn.vs.select(_degree=1)
	for OGNode in OGNodes:
		OGGenes = OGNode['geneIDs'].split(' ')
		COGs = []
		for gene in OGGenes:
			UniprotID = gene.split('|')[1]
			try:
				COGC = cogCPair[OGCOGsPair[UniprotID]]
				COGs.append(COGC)
			except KeyError:
				# pass
				print('Warning %s is not corresponding cog IDs' % UniprotID)
		COGset = set(COGs)
		if COGs:
			hhn.vs[OGNode.index]['COG'] = ' '.join(COGset)
	# print(len(COGs))
	# hhn.vs[OGNode.index]['COGNum'] = len(COGs)
	return hhn


if __name__ == '__main__':
	MappingFile = sys.argv[1:][0]
	COGsC = sys.argv[1:][1]
	hhnFile = sys.argv[1:][2]
	Output = sys.argv[1:][3]
	OG_COG_Pair = MappingCOGForIDs(MappingFile)
	COG_C = COGNumCMapping(COGsC)
	hhnNew = GetGenesCOGs(OG_COG_Pair, COG_C, hhnFile)
	hhnNew.write_gml(Output)
