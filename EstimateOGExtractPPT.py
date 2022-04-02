import os
import time
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


def ReadFile(FileName):
	with open(FileName) as file:
		data = file.readlines()
	return data


def TimeUsedComputation(StartTime, EndTime):
	timeUsed = EndTime - StartTime
	time_str = "time used: {:.0f}h {:.0f}m {:.0f}s"
	TimeUsed = time_str.format(timeUsed // 3600, (timeUsed % 3600) // 60, ((timeUsed % 3600) % 60) % 60)
	return TimeUsed


def GetOGPairs(OGStatFile):
	genesOGPair = {}
	for line in ReadFile(OGStatFile):
		OGName = line.strip().split('\t')[0]
		genesNum = line.strip().split('\t')[2]
		genesList = [gene.split('|')[1] for gene in line.strip().split('\t')[-1].split(' ')]
		for gene in genesList:
			genesOGPair[gene] = (OGName, genesNum)
	return genesOGPair


def GetSequencesIDPair(SeqInfFile):
	SeqInf = MatricesLoad(SeqInfFile)
	since = time.time()
	CDStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(since))
	print('[%s]Get ID Pairs' % CDStartTime)
	IDPair = {}
	for recode, original in SeqInf['GenesRecode'].items():
		if recode == 'G42|g20008' or recode == 'G48|g9369':
			print(recode, original)
		IDPair[original.split('|')[1]] = recode
	done = time.time()
	# print(SeqInf['GenesRecode']['G42|g20008'], SeqInf['GenesRecode']['G48|g9369'])
	print('Get ID Pairs %s' % TimeUsedComputation(since, done))
	return IDPair


def GetRecodeID(GenesID_RecodePair, OGFileSet, OG_GenesID_Pair):
	since = time.time()
	CDStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(since))
	print('[%s]Recode Pairwise' % CDStartTime)
	FNPair = []
	for line in ReadFile(OGFileSet):
		id1 = line.strip().split('\t')[0]
		id2 = line.strip().split('\t')[1]
		label = line.strip().split('\t')[2]
		if id1 not in OG_GenesID_Pair.keys() or id2 not in OG_GenesID_Pair.keys() or id1 not in GenesID_RecodePair.keys() \
				or id2 not in GenesID_RecodePair.keys():
			continue
		else:
			og1 = OG_GenesID_Pair[id1]
			og2 = OG_GenesID_Pair[id2]
			if og1[0] == og2[0]:
				statue = 'SameOG'
			else:
				statue = 'DiffOG'
			recode1 = GenesID_RecodePair[id1]
			recode2 = GenesID_RecodePair[id2]
			# if id1 == 'Q9DCL9' and id2 == 'Q10457':
			FNPair.append((id1, id2, recode1, recode2, label, statue))
	FNPairSplit = [FNPair[i:i + 1000] for i in range(0, len(FNPair), 1000)]
	done = time.time()
	print(FNPairSplit)
	return FNPairSplit


def IdentificationOGOfHHN(queue, hhn, FNs, CCNumDic):
	PID = str(os.getpid())
	since = time.time()
	nodesEvents = []
	nodesPathEvents = []
	for fn in FNs:
		node1 = hhn.vs.select(geneIDs=fn[0])
		node2 = hhn.vs.select(geneIDs=fn[1])
		if len(node1) > 0 and len(node2) > 0:
			shortest_path = node1[0].get_shortest_paths(node2[0])
			if shortest_path[0]:
				try:
					gene1CCNum = CCNumDic[fn[0]]
					gene2CCNum = CCNumDic[fn[1]]
				except KeyError:
					gene1CCNum = 0
					gene2CCNum = 0
				LCAIndex = sorted(shortest_path[0][1:-1], key=lambda X: int(hhn.vs[X]['genesNum']), reverse=True)[0]
				LCAE = hhn.vs[LCAIndex]['Event']
				# shortest_path1 = node1[0].get_shortest_paths(hhn.vs[LCAIndex])
				# shortest_path2 = node2[0].get_shortest_paths(hhn.vs[LCAIndex])
				# if shortest_path1[0][1:-1]:
				# 	Events = [hhn.vs[index]['Event'] for index in shortest_path1[0][1:-1]]
				# 	EventsCount = Counter(Events)
				# 	EventsStat = ['%s:%s' % (event, num) for event, num in EventsCount.items()]
				# 	EventsStatSort = sorted(EventsStat, key=lambda X: X.split(':')[-1], reverse=True)
				# 	LCAEs1 = ' '.join(EventsStatSort)
				# else:
				# 	LCAEs1 = 'No Path'
				# if shortest_path2[0][1:-1]:
				# 	Events = [hhn.vs[index]['Event'] for index in shortest_path2[0][1:-1]]
				# 	EventsCount = Counter(Events)
				# 	EventsStat = ['%s:%s' % (event, num) for event, num in EventsCount.items()]
				# 	EventsStatSort = sorted(EventsStat, key=lambda X: int(X.split(':')[-1]), reverse=True)
				# 	LCAEs2 = ' '.join(EventsStatSort)
				# else:
				# 	LCAEs2 = 'No Path'
				# childs = []
				# neighb = hhn.vs[LCAIndex].neighbors()
				# # ChildModularity = []
				# for child in neighb:
				# 	if child['genesNum'] < hhn.vs[LCAIndex]['genesNum']:
				# 		childs.append(child['genomesNum'])
						# ChildModularity.append(str(float(child['modularity'])))
				# Diff = int(abs((childs[0] - childs[1])))
				# Total = int(hhn.vs[LCAIndex]['genomesNum'])
				# DiffPro = Diff / Total
				# ChildModularityStr = ' '.join(ChildModularity)
				Modularity = hhn.vs[LCAIndex]['modularity']
				nodesEvents.append('%s\t%s\t%s\t%d\t%d\n' % ('\t'.join(fn), str(Modularity), LCAE, gene1CCNum, gene2CCNum))
				nodesPathEvents.append('%s\t%s\t%s\t%d\t%d\n' % ('\t'.join(fn), str(Modularity), LCAE, gene1CCNum, gene2CCNum))
			else:
				nodesEvents.append('%s\tNot Connect\n' % '\t'.join(fn))
				nodesPathEvents.append('%s\tNot Connect\n' % '\t'.join(fn))
		elif len(node1) > 0 and len(node1) == 0:
			nodesEvents.append('\t'.join(fn) + '\t' + 'Nodes were inseparable' + '\n')
			nodesPathEvents.append('\t'.join(fn) + '\t' + 'Nodes were inseparable' + '\n')
		elif len(node1) > 0 and len(node1) == 0:
			nodesEvents.append('\t'.join(fn) + '\t' + 'Nodes were inseparable' + '\n')
			nodesPathEvents.append('\t'.join(fn) + '\t' + 'Nodes were inseparable' + '\n')
		else:
			nodesEvents.append('\t'.join(fn) + '\t' + 'Nodes were inseparable' + '\n')
			nodesPathEvents.append('\t'.join(fn) + '\t' + 'Nodes were inseparable' + '\n')
	done = time.time()
	# print(len(LCAIndexList), len(nodesEvents), len(nodesPathEvents))
	# print(LCAIndexList, nodesEvents, nodesPathEvents)
	print('%s Find Node Event %s' % (PID, TimeUsedComputation(since, done)))
	queue.put((nodesEvents, nodesPathEvents))


def GetNodeEvent(queue, taskNum):
	progressbar_widgets_set = [
		'Get Nodes Events:', progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(),
		' | ',
		progressbar.ETA()]
	bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=taskNum)
	Events = ''
	NodesPathEvent = ''
	done = 0
	bar.start()
	while True:
		NodeEvent, NodesPathEvents = queue.get()
		for event, path in zip(NodeEvent, NodesPathEvents):
			Events += event
			NodesPathEvent += path
		done += 1
		if done >= taskNum:
			break
		bar.update(done)
	bar.finish()
	return Events, NodesPathEvent


def IdentificationOGOfHHNPP(hmFile, FNsList, Threads, fileName, NumThreads):
	pro = mp.Pool(int(Threads))
	MyQueue = mp.Manager().Queue()
	HHN = igraph.read(hmFile)
	CCNumDic = FindConnectedComponents(HHN, NumThreads)
	for FN in FNsList:
		pro.apply_async(func=IdentificationOGOfHHN, args=(MyQueue, HHN, FN, CCNumDic))
		# IdentificationOGOfHHN(MyQueue, HHN, FN, CCNumDic)
	NodesEvent, NodesPathEvent = GetNodeEvent(MyQueue, len(FNsList))
	pro.close()
	pro.join()
	OutputFile(fileName, NodesEvent)
	OutputFile('Path%s' % fileName, NodesPathEvent)


def FindConnectedComponents(hhn, NumThreads):
	since = time.time()
	hhnComponents = hhn.components()
	NodesIDsDic = {}
	for num, cc in enumerate(hhnComponents):
		# if 5 < hhn.subgraph(cc).vcount() < 1000:
		MaxGenesNum = sorted(hhn.subgraph(cc).vs['genesNum'], key=lambda X: int(X), reverse=True)[0]
		if int(MaxGenesNum) <= int(NumThreads):
			TerminateNodes = hhn.subgraph(cc).vs.select(_degree=1)
			for Nodes in TerminateNodes:
				for geneID in Nodes['geneIDs'].strip().split(' '):
					NodesIDsDic[geneID] = MaxGenesNum
	done = time.time()
	print('Get CC Numbers %s' % TimeUsedComputation(since, done))
	return NodesIDsDic



if __name__ == '__main__':
	fnList = sys.argv[1:][0]
	seqIDs = sys.argv[1:][1]
	OGStat = sys.argv[1:][2]
	hm = sys.argv[1:][3]
	threads = sys.argv[1:][4]
	output = sys.argv[1:][5]
	NumThreads = sys.argv[1:][6]
	IDPairs = GetSequencesIDPair(seqIDs)
	OGGroupPair = GetOGPairs(OGStat)
	FNPairs = GetRecodeID(IDPairs, fnList, OGGroupPair)
	# IdentificationOGOfHHNPP(hm, FNPairs, threads, output, NumThreads)
