import igraph
import sys


def OutputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def ReadFile(FileName):
	with open(FileName) as file:
		data = file.readlines()
	return data


def GetOGPairs(OGStatFile):
	genesOGPair = {}
	for line in ReadFile(OGStatFile):
		OGName = line.strip().split('\t')[0]
		genesList = [gene.split('|')[1] for gene in line.strip().split('\t')[-1].split(' ')]
		for gene in genesList:
			genesOGPair[gene] = OGName
	return genesOGPair


def GetSequencesIDPair(SequenceFile):
	IDPair = {}
	for line in ReadFile(SequenceFile):
		seqID = line.strip().split('\t')[0].split('|')[1]
		recodeID = line.strip().split('\t')[1]
		IDPair[seqID] = recodeID
	return IDPair


def GetRecodeID(idPair, OGFileSet):
	FNPair = []
	for line in ReadFile(OGFileSet):
		id1 = line.strip().split('\t')[0]
		id2 = line.strip().split('\t')[1]
		label = line.strip().split('\t')[2]
		try:
			recode1 = idPair[id1]
			recode2 = idPair[id2]
			FNPair.append((recode1, recode2, label))
		except KeyError:
			pass
	return FNPair


def IdentificationOGOfHHN(hmFile, FNs, fileName):
	nodesEvents = ''
	hhn = igraph.read(hmFile)
	for fn in FNs:
		node1 = hhn.vs.select(geneIDs=fn[0])
		node2 = hhn.vs.select(geneIDs=fn[1])
		if len(node1) > 0 and len(node2) > 0:
			shortest_path = node1[0].get_shortest_paths(node2[0])
			if shortest_path[0]:
				LCAIndex = sorted(shortest_path[0][1:-1], key=lambda X: int(hhn.vs[X]['genesNum']), reverse=True)[0]
				LCAE = hhn.vs[LCAIndex]['Event']
				neighb = hhn.vs[LCAIndex].neighbors()
				childs = []
				for child in neighb:
					if child['genesNum'] < hhn.vs[LCAIndex]['genesNum']:
						childs.append(str(child['genomesNum']))
				nodesEvents += '\t'.join(fn) + '\t' + LCAE + '\t' + str(hhn.vs[LCAIndex]['genomesNum']) + '\t' +'\t'.join(
					childs) + '\n'
			else:
				nodesEvents += '\t'.join(fn) + '\t' + 'Not Connected' + '\n'
		else:
			nodesEvents += '\t'.join(fn) + '\t' + 'Not Found' + '\n'
	OutputFile(fileName, nodesEvents)


if __name__ == '__main__':
	fnList = sys.argv[1:][0]
	seqIDs = sys.argv[1:][1]
	hm = sys.argv[1:][2]
	output = sys.argv[1:][3]
	IDPairs = GetSequencesIDPair(seqIDs)
	FNPairs = GetRecodeID(IDPairs, fnList)
	IdentificationOGOfHHN(hm, FNPairs, output)
