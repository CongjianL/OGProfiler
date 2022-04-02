import time
import igraph
import sys
import pickle as pic


def MatricesLoad(MatriceName):
	with open(MatriceName, 'rb') as picFile:
		M = pic.load(picFile)
	return M


def OutputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def TimeUsedComputation(StartTime, EndTime):
	timeUsed = EndTime - StartTime
	time_str = "time used: {:.0f}h {:.0f}m {:.0f}s"
	TimeUsed = time_str.format(timeUsed // 3600, (timeUsed % 3600) // 60, ((timeUsed % 3600) % 60) % 60)
	return TimeUsed


def ReadFile(FileName):
	with open(FileName) as file:
		data = file.readlines()
	return data


def BenchOGResult(OrthoGroups, OGsNotFound, RecodeSeqInf):
	# OGsNotFoundList = [line.strip().split('\t')[0] for line in ReadFile(OGsNotFound) if
	# 				   line.strip().split('\t')[3] == 'Not Found']
	BenchOGsPair = {}
	for line in ReadFile(OrthoGroups):
		ogName = line.strip().split(':')[0]
		genes = line.strip().split(':')[1].strip().split(' ')
		if len(genes) > 1:
			genesString = ' '.join(sorted([RecodeSeqInf[gene] for gene in genes]))
			BenchOGsPair[ogName] = genesString
	return BenchOGsPair


def GetSequencesIDPair(SeqInfFile):
	SeqInf = MatricesLoad(SeqInfFile)
	IDPair = {}
	for recode, original in SeqInf['GenesRecode'].items():
		IDPair[original] = recode
	# print(SeqInf['GenesRecode']['G42|g20008'], SeqInf['GenesRecode']['G48|g9369'])
	return IDPair


def GetHHNGeneIDs(graphFile):
	since = time.time()
	NodeIDGenesPair = {}
	for node in graphFile.vs:
		Index = node.index
		GeneIDs = node['geneIDs']
		NodeIDGenesPair[Index] = set(GeneIDs.strip().split(' '))
	done = time.time()
	Times = TimeUsedComputation(since, done)
	print('Get Nodes Pairwise %s' % Times)
	return NodeIDGenesPair


def FindOGsLocation(OrthgroupsFile, OGsNotFound, SeqInfFile, hhnFile, outPut):
	SeqRecodeInf = GetSequencesIDPair(SeqInfFile)
	OGsInf = BenchOGResult(OrthgroupsFile, OGsNotFound, SeqRecodeInf)
	hhn = igraph.read(hhnFile)
	IndexGenes = GetHHNGeneIDs(hhn)
	Stat = ''
	Singleton = ''
	for name, ogGenes in OGsInf.items():
		# OGSelect = hhn.vs.select(geneIDs=ogGenes)
		genesNum = len(ogGenes.split(' '))
		ogGenesSet = set(ogGenes.split(' '))
		genomeSet = set([gene.split('|')[0] for gene in ogGenes.split(' ')])
		OverlapOGs = []
		for NodeId, genesSet in IndexGenes.items():
			consist = genesSet & ogGenesSet
			if len(consist) > 0:
				GenesNum = hhn.vs[NodeId]['genesNum']
				pro1 = len(consist) / genesNum
				pro2 = len(consist)/ GenesNum
				OverlapOGs.append((NodeId, len(consist), pro1, pro2))
		if OverlapOGs:
			MaxOverlap = sorted(OverlapOGs, key=lambda X: (X[2], X[3]), reverse=True)[0]
			MaxNodes = hhn.vs[MaxOverlap[0]]
			childs = []
			for child in MaxNodes.neighbors():
				if child['genesNum'] < MaxNodes['genesNum']:
					childs.append(child['genesNum'])
			proportion = MaxOverlap[2]
			Event = MaxNodes['Event']
			MaxGenesNum = MaxNodes['genesNum']
			ChildsGenes = ' '.join(['%d' % num for num in sorted(childs)])
			proportion2 = MaxOverlap[1] / MaxGenesNum
			# MaxGenes = MaxNodes['geneIDs']
			Stat += '%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\n' % (
				name, Event, len(genomeSet), genesNum, MaxGenesNum, MaxOverlap[1], proportion, proportion2, ChildsGenes)
		else:
			genomeSet = set([gene.split('|')[0] for gene in ogGenes.split(' ')])
			if len(genomeSet) > 1:
				print(genomeSet)
			Singleton += '%s\t%d\t1\t%d\n' % (name, genesNum, len(genomeSet))
	# if len(OGSelect) > 0:
	# 	NodesEvents = OGSelect[0]['Event']
	# 	Stat += '%s\t%d\t%s\tFound\n' % (name, genesNum, NodesEvents)
	# 	if len(OGSelect) > 1:
	# 		print(len(OGSelect))
	# else:
	# 	Stat += '%s\t%d\tNone\tNot Found\n' % (name, genesNum)
	OutputFile('%s.single.txt' % outPut.split('.')[0], Singleton)
	OutputFile(outPut, Stat)


if __name__ == '__main__':
	OrthoFinder = sys.argv[1:][0]
	SeqInf = sys.argv[1:][1]
	hm = sys.argv[1:][2]
	NotFound = sys.argv[1:][3]
	out = sys.argv[1:][4]
	FindOGsLocation(OrthoFinder, NotFound, SeqInf, hm, out)
