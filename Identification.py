import sys
import pickle as pic


def ReadFile(FileName):
	with open(FileName) as file:
		data = file.readlines()
	return data


def MatricesLoad(MatriceName):
	with open(MatriceName, 'rb') as picFile:
		M = pic.load(picFile)
	return M


def GetSequencesIDPair(SeqInfFile):
	SeqInf = MatricesLoad(SeqInfFile)
	IDPair = {}
	for recode, original in SeqInf['GenesRecode'].items():
		# print(recode, original)
		IDPair[original.split('|')[1]] = recode
	return IDPair


def Comparison(Pairwise, Nodes, SeqInfFile):
	SeqInf = MatricesLoad(SeqInfFile)
	IDPair = GetSequencesIDPair(SeqInfFile)
	Pairset = set([line.strip() for line in ReadFile(Pairwise)])
	NodesRecode = set()
	for line in ReadFile(Nodes):
		gene1 = line.strip().split('\t')[0]
		gene2 = line.strip().split('\t')[1]
		label = line.strip().split('\t')[2]
		try:
			recode1, recode2 = SeqInf['GenesRecode'][gene1].split('|')[1], SeqInf['GenesRecode'][gene2].split('|')[1]
		except KeyError:
			continue
		NewLine = '\t'.join([recode1, recode2, label])
		NodesRecode.add(NewLine)
	Diff = Pairset - NodesRecode
	print(Diff)
	if Diff == set():
		print('Same')
	else:
		for i in Diff:
			try:
				print(i.split('\t')[0], i.split('\t')[1], IDPair[i.split('\t')[0]], IDPair[i.split('\t')[1]], i.split('\t')[2])
			except KeyError:
				continue



if __name__ == '__main__':
	Seq = sys.argv[1:][0]
	P1 = sys.argv[1:][1]
	P2 = sys.argv[1:][2]
	Comparison(P1, P2, Seq)
