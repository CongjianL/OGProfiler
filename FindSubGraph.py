import os
import time
from itertools import combinations
import igraph
import sys
import multiprocessing as mp
import progressbar
from collections import Counter
import pickle as pic


def MatricesLoad(MatriceName):
	with open(MatriceName, 'rb') as picFile:
		M = pic.load(picFile)
	return M


def OutputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def GetGenesIDs(selectedNode):
	neighbors = [selectedNode]
	childNode = []
	while neighbors:
		neighbors = GetGenesID(neighbors, childNode)
	return childNode


def GetGenesID(nodes, deletedNode):
	nodeSelected = []
	for node in nodes:
		if node['geneIDs'] == 'None' or node['geneIDs'] is None:
			for neighbor in node.neighbors():
				if neighbor['genesNum'] < node['genesNum']:
					nodeSelected.append(neighbor)
					deletedNode.append(neighbor)

	return nodeSelected


def ReadFile(FileName):
	with open(FileName) as file:
		data = file.readlines()
	return data


def TimeUsedComputation(StartTime, EndTime):
	timeUsed = EndTime - StartTime
	time_str = "time used: {:.0f}h {:.0f}m {:.0f}s"
	TimeUsed = time_str.format(timeUsed // 3600, (timeUsed % 3600) // 60, ((timeUsed % 3600) % 60) % 60)
	return TimeUsed


def GetOGPairs(OGStatFileName):
	OGGenesPair = {}
	genesOGPair = {}
	for line in ReadFile(OGStatFileName):
		OGName = line.strip().split('\t')[0]
		genesNum = line.strip().split('\t')[2]
		genesList = [gene.split('|')[1] for gene in line.strip().split('\t')[-1].split(' ')]
		OGGenesPair[OGName] = genesList
		for gene in genesList:
			genesOGPair[gene] = (OGName, genesNum)
	return genesOGPair, OGGenesPair


def GetSequencesIDPair(SeqInfFile):
	SeqInf = MatricesLoad(SeqInfFile)
	since = time.time()
	CDStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(since))
	print('[%s]Get ID Pairs' % CDStartTime)
	IDPair = {}
	for recode, original in SeqInf['GenesRecode'].items():
		# print(recode, original)
		IDPair[original.split('|')[1]] = recode
	done = time.time()
	print('Get ID Pairs %s' % TimeUsedComputation(since, done))
	return IDPair


def GetRecodeID(idPair, OGFileSet, OGsGroupPair, OGGenesPair, SeqInfFile, findLabel):
	since = time.time()
	SeqInf = MatricesLoad(SeqInfFile)
	CDStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(since))
	print('[%s]Recode Pairwise' % CDStartTime)
	OGsStat = []
	OGFN = ''
	for line in ReadFile(OGFileSet):
		genesID1 = line.strip().split('\t')[0]
		genesID2 = line.strip().split('\t')[1]
		label = line.strip().split('\t')[2]
		statue = line.strip().split('\t')[3]
		try:
			Event = line.strip().split('\t')[5]
		except IndexError:
			continue
		if statue == 'SameOG' and Event == 'III-1' and label == findLabel:
			genome1 = genesID1.split('|')[0]
			genome2 = genesID2.split('|')[0]
			if genome1 != genome2:
				recodeGenes1 = SeqInf['GenesRecode'][genesID1]
				og = OGsGroupPair[recodeGenes1.split('|')[1]]
				OGFN += '%s\t%s\t%s\n' % (og[0], genesID1, genesID2)
			# if og[0] == 'OG4544':
			# 	print(genesID1, genesID2)
				OGsStat.append(og)
	OGStatNum = \
		sorted([(og, num) for og, num in Counter(OGsStat).items() if num <= 10], key=lambda X: X[0][1], reverse=True)
	OGGenesPairRecode = []
	for OG, num in OGStatNum:
		name, genesNum = OG
		idsRecode = [idPair[ID] for ID in OGGenesPair[name]]
		OGGenesPairRecode.append((OG, idsRecode))
	OutputFile('OG_%s.txt' % findLabel, OGFN)
	return OGGenesPairRecode


def ListCombination(List):
	CombinationList = []
	for I in range(len(List) - 1):
		for J in range(I + 1, len(List)):
			CombinationList.append((List[I], List[J]))
	return CombinationList


def FindNodes(queue, hhn, seqPairsList):
	LCAIndexSet = []
	since = time.time()
	print(os.getpid())
	for pair in seqPairsList:
		q, r = pair
		node1 = hhn.vs.select(geneIDs=q)
		node2 = hhn.vs.select(geneIDs=r)
		if len(node1) > 0 and len(node2) > 0:
			shortest_path = node1[0].get_shortest_paths(node2[0])
			if shortest_path[0]:
				LCAIndex = sorted(shortest_path[0][1:-1], key=lambda X: int(hhn.vs[X]['genesNum']), reverse=True)[0]
				LCAIndexSet.append(LCAIndex)
	done = time.time()
	print(TimeUsedComputation(since, done))
	queue.put(LCAIndexSet)


def GetSet(queue, TaskNum):
	progressbar_widgets_set = [
		'Get Nodes Events:', progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(),
		' | ',
		progressbar.ETA()]
	bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=TaskNum)
	bar.start()
	done = 1
	IndexList = []
	while True:
		indexs = queue.get()
		if indexs:
			IndexList.extend(indexs)
		bar.update(done)
		print(done)
		if done >= TaskNum:
			break
		done += 1
	bar.finish()
	return IndexList


def FindSubGraphFN(hhnFile, MaxRecodeOG, threads, outputFile):
	hhn = igraph.read(hhnFile)
	print(len(MaxRecodeOG))
	seqPairsList = ListCombination(MaxRecodeOG)
	print(len(seqPairsList))
	seqPairsLists = [seqPairsList[pairIndex: pairIndex + 100] for pairIndex in range(0, len(seqPairsList), 100)]
	pro = mp.Pool(int(threads))
	MyQueue = mp.Manager().Queue()
	for seqPairs in seqPairsLists:
		# INDEX = FindNodes(MyQueue, hhn, seqPairs)
		# if INDEX:
		# 	LCAIndex.extend(INDEX)
		pro.apply_async(func=FindNodes, args=(MyQueue, hhn, seqPairs))
	LCAIndex = GetSet(MyQueue, len(seqPairsLists))
	LCAIndexSet = set(LCAIndex)
	pro.close()
	pro.join()
	LCAMaxIndex = sorted(list(LCAIndexSet), key=lambda X: int(hhn.vs[X]['genesNum']), reverse=True)[0]
	LCANode = hhn.vs[LCAMaxIndex]
	childNodes = GetGenesIDs(LCANode)
	childNodes.append(LCANode)
	Childs = [child.index for child in childNodes]
	MaxSub = hhn.subgraph(Childs)
	MaxSub.write_gml(outputFile)


def FindConnectedComponents(hhn, NumThreads):
	hhnComponents = hhn.components()
	CCountDic = {}
	NodesIDsDic = {}
	for num, cc in enumerate(hhnComponents):
		if 5 < hhn.subgraph(cc).vcount() < 1000:
			MaxGenesNum = sorted(hhn.subgraph(cc).vs['genesNum'], key=lambda X: int(X), reverse=True)[0]
			if int(MaxGenesNum) <= int(NumThreads):
				CCountDic[num] = cc
				TerminateNodes = hhn.subgraph(cc).vs.select(_degree=1)
				genesIDsList = []
				for Nodes in TerminateNodes:
					for geneID in Nodes['geneIDs'].strip().split(' '):
						genesIDsList.append(geneID)
				NodesIDsDic[num] = set(genesIDsList)
	return CCountDic, NodesIDsDic


def FindSmallCCNodes(RecodeOGs, hhnFile, genesNumThreads, subDir):
	try:
		os.mkdir(subDir)
	except FileExistsError:
		pass
	HHN = igraph.read(hhnFile)
	CCDic, IDDic = FindConnectedComponents(HHN, genesNumThreads)
	for og, genesList in RecodeOGs:
		ogName, genesNum = og
		genesSet = set(genesList)
		for CCNum, ids in IDDic.items():
			uninoSet = genesSet & ids
			if len(uninoSet) > 0:
				cc = CCDic[CCNum]
				SubGraph = HHN.subgraph(cc)
				subName = os.path.join(subDir, '%s_%d_%d.gml' % (ogName, CCNum, len(uninoSet)))
				SubGraph.write_gml(subName)


if __name__ == '__main__':
	OGPairList = sys.argv[1:][0]
	seqInfFile = sys.argv[1:][1]
	OGStatFile = sys.argv[1:][2]
	hm = sys.argv[1:][3]
	threads = sys.argv[1:][4]
	FindLabel = sys.argv[1:][5]
	output = sys.argv[1:][6]
	genesOGDic, OGGenesDic = GetOGPairs(OGStatFile)
	idPair = GetSequencesIDPair(seqInfFile)
	MaxRecode = GetRecodeID(idPair, OGPairList, genesOGDic, OGGenesDic, seqInfFile, FindLabel)
	# FindSubGraphFN(hm, MaxRecode, threads, output)
	FindSmallCCNodes(MaxRecode, hm, threads, output)
