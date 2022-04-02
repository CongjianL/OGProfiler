import sys


def InputFile(fileName):
	with open(fileName) as f:
		data = f.readlines()
	return data


def OutputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def GetOGPairDataSet(OGFile):
	G = set()
	RO = set()
	RP = set()
	for line in InputFile(OGFile):
		if 'TreeFamA' in line:
			q = line.strip().split('\t')[1]
			r = line.strip().split('\t')[2]
			label = line.strip().split('\t')[3]
			pairwise = '%s\t%s' % (q, r)
			pairwiseR = '%s\t%s' % (r, q)
			G.add(q)
			G.add(r)
			if label == 'TP' or label == 'FN':
				RO.add(pairwise)
				RO.add(pairwiseR)
			elif label == 'FP' or label == 'TN':
				RP.add(pairwise)
				RP.add(pairwiseR)
			else:
				print(line)
	return G, RO, RP


def GetPredictionOGSet(OGPairwise, OrthoSet):
	PredictionSet = set()
	for line in InputFile(OGPairwise):
		gene1 = line.strip().split('\t')[0]
		gene2 = line.strip().split('\t')[1]
		if gene1 in OrthoSet and gene2 in OrthoSet:
			pairwise = '%s\t%s' % (gene1, gene2)
			pairwiseR = '%s\t%s' % (gene2, gene1)
			PredictionSet.add(pairwise)
			PredictionSet.add(pairwiseR)
	return PredictionSet


def CalculatePR(predictionSet, OrthoPair, NonOrthoPair, Output):
	TP = predictionSet & OrthoPair
	FP = predictionSet & NonOrthoPair
	FN = OrthoPair - predictionSet
	TN = NonOrthoPair - predictionSet
	PPV = len(TP) / (len(TP) + len(FP))
	TPR = len(TP) / (len(TP) + len(FN))
	TPS = '\n'.join(['%s\tTP' % i for i in TP])
	FPS = '\n'.join(['%s\tFP' % i for i in FP])
	TNS = '\n'.join(['%s\tTN' % i for i in TN])
	FNS = '\n'.join(['%s\tFN' % i for i in FN])
	Estimate = '\n'.join([TPS, FPS, TNS, FNS])
	OutputFile(Output, Estimate)
	print('PPV: %f\t TRP: %f' % (PPV, TPR))


if __name__ == '__main__':
	OGPairDataSet = sys.argv[1:][0]
	OGPrediction = sys.argv[1:][1]
	OutPut = sys.argv[1:][2]
	GenesDataSet, Orthologs, NonOrthologs = GetOGPairDataSet(OGPairDataSet)
	predictionOGs = GetPredictionOGSet(OGPrediction, GenesDataSet)
	CalculatePR(predictionOGs, Orthologs, NonOrthologs, OutPut)
