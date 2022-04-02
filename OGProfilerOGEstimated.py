import sys


def OutputFile(fileName, content):
	with open(fileName, 'w') as f:
		f.writelines(content)


def ReadFile(FileName):
	with open(FileName) as file:
		data = file.readlines()
	return data


def BenchOGResult(OrthoGroups, FileType):
	SingleCopy = {}
	BenchOGsPair = {}
	OrthoSingle = ["OG0004241", "OG0004245", "OG0004247", "OG0004248", "OG0004258", "OG0004259", "OG0004266",
				   "OG0004274", "OG0004276", "OG0004281", "OG0004283", "OG0004289", "OG0004290", "OG0004292",
				   "OG0004295", "OG0004300", "OG0004301", "OG0004302", "OG0004306", "OG0004310", "OG0004316",
				   "OG0004319", "OG0004320", "OG0004326", "OG0004329", "OG0004331", "OG0004332", "OG0004333",
				   "OG0004334", "OG0004338", "OG0004341", "OG0004343", "OG0004344", "OG0004347", "OG0004351",
				   "OG0004354", "OG0004365", "OG0004372", "OG0004403", "OG0004407", "OG0004411", "OG0004416",
				   "OG0004422", "OG0004428", "OG0004431", "OG0004434", "OG0004437", "OG0004445", "OG0004451",
				   "OG0004452", "OG0004461", "OG0004462", "OG0004464", "OG0004468", "OG0004473", "OG0004475",
				   "OG0004476", "OG0004477", "OG0004478", "OG0004483", "OG0004485", "OG0004486", "OG0004489",
				   "OG0004491", "OG0004495", "OG0004500", "OG0004502", "OG0004503", "OG0004506", "OG0004507",
				   "OG0004513", "OG0004514", "OG0004519", "OG0004521", "OG0004522", "OG0004525", "OG0004526",
				   "OG0004530", "OG0004534", "OG0004535", "OG0004541", "OG0004543", "OG0004548", "OG0004555",
				   "OG0004556", "OG0004559", "OG0004561", "OG0004562", "OG0004563", "OG0004565", "OG0004566",
				   "OG0004569", "OG0004571", "OG0004574", "OG0004576", "OG0004578", "OG0004580", "OG0004583",
				   "OG0004584", "OG0004593", "OG0004595", "OG0004600", "OG0004601", "OG0004602", "OG0004603",
				   "OG0004610", "OG0004612", "OG0004615", "OG0004616", "OG0004617", "OG0004619", "OG0004622",
				   "OG0004624", "OG0004625", "OG0004628", "OG0004635", "OG0004637", "OG0004638", "OG0004639",
				   "OG0004641", "OG0004644", "OG0004645", "OG0004647", "OG0004649", "OG0004650", "OG0004653",
				   "OG0004657", "OG0004658", "OG0004659", "OG0004662", "OG0004663", "OG0004664", "OG0004665",
				   "OG0004667", "OG0004669", "OG0004670", "OG0004671", "OG0004674", "OG0004680", "OG0004682",
				   "OG0004684", "OG0004685", "OG0004688", "OG0004689", "OG0004691", "OG0004692", "OG0004694",
				   "OG0004695", "OG0004697", "OG0004699", "OG0004701", "OG0004702", "OG0004705", "OG0004707",
				   "OG0004709", "OG0004711", "OG0004712", "OG0004715", "OG0004717", "OG0004719", "OG0004721",
				   "OG0004722", "OG0004723", "OG0004724", "OG0004725", "OG0004726", "OG0004732", "OG0004738",
				   "OG0004739", "OG0004740", "OG0004741", "OG0004742", "OG0004743", "OG0004747", "OG0004748",
				   "OG0004753", "OG0004754", "OG0004755", "OG0004758", "OG0004759", "OG0004760", "OG0004761",
				   "OG0004764", "OG0004766", "OG0004767", "OG0004768", "OG0004769", "OG0004772", "OG0004773",
				   "OG0004774", "OG0004775", "OG0004781", "OG0004782", "OG0004785", "OG0004787", "OG0004790",
				   "OG0004796", "OG0004799", "OG0004806", "OG0004810", "OG0004811", "OG0004815", "OG0004816",
				   "OG0004818", "OG0004820", "OG0004822", "OG0004826", "OG0004828", "OG0004830", "OG0004831",
				   "OG0004836", "OG0004838", "OG0004840", "OG0004842", "OG0004843", "OG0004845", "OG0004847",
				   "OG0004849", "OG0004852", "OG0004853", "OG0004855", "OG0004856", "OG0004857", "OG0004858",
				   "OG0004863", "OG0004864", "OG0004865", "OG0004866", "OG0004873", "OG0004875", "OG0004877",
				   "OG0004878", "OG0004879", "OG0004882", "OG0004884", "OG0004889", "OG0004890", "OG0004896",
				   "OG0004898", "OG0004899", "OG0004900", "OG0004901", "OG0004903", "OG0004905", "OG0004906",
				   "OG0004907", "OG0004908", "OG0004909", "OG0004910", "OG0004911", "OG0004913", "OG0004914",
				   "OG0004915", "OG0004918", "OG0004920", "OG0004921", "OG0004923", "OG0004926", "OG0004927",
				   "OG0004928", "OG0004930", "OG0004931", "OG0004932", "OG0004933", "OG0004934", "OG0004935",
				   "OG0004938", "OG0004940", "OG0004943", "OG0004945", "OG0004947", "OG0004948", "OG0004951",
				   "OG0004953", "OG0004955", "OG0004956", "OG0004957", "OG0004958", "OG0004960", "OG0004961",
				   "OG0004963", "OG0004964", "OG0004966", "OG0004967", "OG0004969", "OG0004971", "OG0004974",
				   "OG0004979", "OG0004983", "OG0004984", "OG0004986", "OG0004989", "OG0004990", "OG0004991",
				   "OG0004992", "OG0004994", "OG0004996", "OG0004997", "OG0005001", "OG0005003", "OG0005004",
				   "OG0005005", "OG0005006", "OG0005007", "OG0005009", "OG0005010", "OG0005013", "OG0005014",
				   "OG0005016", "OG0005017", "OG0005019", "OG0005021", "OG0005031"]
	for line in ReadFile(OrthoGroups):
		if FileType == 'OrthoFinder':
			ogName = line.strip().split(':')[0]
			genes = line.strip().split(':')[1].strip().split(' ')
			if ogName in OrthoSingle:
				SingleCopy[ogName] = set(genes)
		else:
			ogName = line.strip().split('\t')[0]
			genomeNum = int(line.strip().split('\t')[1])
			genesNum = int(line.strip().split('\t')[2])
			genes = line.strip().split('\t')[-1].strip().split(' ')
			if genomeNum == genesNum and genomeNum == 10:
				SingleCopy[ogName] = set(genes)
		genesSet = set(genes)
		BenchOGsPair[ogName] = genesSet
	return BenchOGsPair, SingleCopy


def FindLargestOverlap(ogSet_Ortho, OGPairOGProfiler):
	OverlapList = []
	for ogName_OG, ogSet_OG in OGPairOGProfiler.items():
		consistent = ogSet_OG & ogSet_Ortho
		if len(consistent) > 0:
			proportions1 = len(consistent) / len(ogSet_Ortho)
			proportions2 = len(consistent) / len(ogSet_OG)
			OverlapList.append((ogName_OG, ogSet_OG, len(consistent), proportions1, proportions2))
	if OverlapList:
		LargestOverlap = sorted(OverlapList, key=lambda X: (X[3], X[4]), reverse=True)[0]
		Name, Set, Num, Consistent, ConsistentPro = LargestOverlap
		return Name, Set, Num, Consistent, ConsistentPro
	else:
		return None


def Comparison(OGPairOrthoFinder, OGPairOGProfiler):
	OGsNum = 0
	GenesNum = 0
	MissingGenesNum = 0
	ErrorGenesNum = 0
	CorrectGenesNum = 0
	OGWithoutMiss = []
	Estimated = ''
	for ogName_Ortho, ogSet_Ortho in OGPairOrthoFinder.items():
		if len(ogSet_Ortho) > 1:
			genomeNum = len(set([gene.split('|')[-1].split('_')[-1] for gene in ogSet_Ortho]))
			OGsNum += 1
			GenesNum += len(ogSet_Ortho)
			if FindLargestOverlap(ogSet_Ortho, OGPairOGProfiler) is None:
				Estimated += '%s\tOG Not Found\n' % ogName_Ortho
			else:
				SelectedOGName, SelectedOGSet, OverlapNum, ConsistentProportion1, ConsistentProportion2 = FindLargestOverlap(
					ogSet_Ortho, OGPairOGProfiler)
				Missing = ogSet_Ortho - SelectedOGSet
				Error = SelectedOGSet - ogSet_Ortho
				if len(Missing) == 0 and len(Error) == 0:
					OGWithoutMiss.append(ogName_Ortho)
				elif len(Missing) > 0 and len(Error) == 0:
					MissingGenesNum += len(Missing)
				elif len(Missing) == 0 and len(Error) > 0:
					ErrorGenesNum += len(Error)
				else:
					MissingGenesNum += len(Missing)
					ErrorGenesNum += len(Error)
				CorrectGenesNum += OverlapNum
				Correct_Num = OverlapNum
				Miss_Num = len(Missing)
				Error = len(Error)
				Estimated += '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n' % (ogName_Ortho, SelectedOGName, genomeNum, len(ogSet_Ortho), len(SelectedOGSet), Correct_Num, Miss_Num, Error)
	print('Total Numbers of OGs: %d\nTotal Numbers of Genes: %d' % (OGsNum, GenesNum))
	print('The Numbers of Correctly Assigned Genes: %d' % CorrectGenesNum)
	print('The Numbers of Missing Assigned Genes: %d' % MissingGenesNum)
	print('The Numbers of Erroneously Assigned Genes: %d' % ErrorGenesNum)
	print('The Numbers of orthogroups without any errors: %d' % len(OGWithoutMiss))
	return Estimated


def EstimatedOGProfiler(orthofinder, ogprofiler, outName, compareType):
	OGPairOrthoFinder, SingleCopyOrtho = BenchOGResult(orthofinder, 'OrthoFinder')
	OGPairOGProfiler, SingleCopyOGProfiler = BenchOGResult(ogprofiler, 'OGProfiler')
	if compareType == 'Single':
		Estimated = Comparison(SingleCopyOrtho, OGPairOGProfiler)
	else:
		Estimated = Comparison(OGPairOrthoFinder, OGPairOGProfiler)
	OutputFile(outName, Estimated)


if __name__ == '__main__':
	OGProfilerPrediction = sys.argv[1:][0]
	OrthoFinderPrediction = sys.argv[1:][1]
	OUT = sys.argv[1:][2]
	Type = sys.argv[1:][3]
	EstimatedOGProfiler(OrthoFinderPrediction, OGProfilerPrediction, OUT, Type)
