import multiprocessing as mp
import os
import pickle as pic
import subprocess
import sys
import time
from optparse import OptionParser, OptionGroup
import igraph
import leidenalg as la
import numpy as np
import progressbar
import pyfasta as pyf
import scipy.sparse as sparse
from scipy.optimize import curve_fit
from scipy.special import comb
import random
from itertools import combinations


def get_parameters():
    usage = "OGProfiler.py -i <input dir> -o graph name <Options>"
    opt = OptionParser(usage=usage)
    group0 = OptionGroup(opt, "General options")
    group0.add_option(
        '-i', '--in', type=str, dest='input_dir', default=False,
        help='specify the directory including all genome files')
    group0.add_option(
        '-o', '--out', type=str, dest='output_graph', default=False, help="specify the graph file name")
    group0.add_option(
        '-s', '--search_method', type=str, dest='search_method', default='diamond', help='Sequence searching methods'
                                                                                         ': blastp, mmseqs, diamond. Both mmseqs and diamond are sensitive mode (default: diamond)')
    group0.add_option(
        '-e', '--evalue', type=float, dest='e_values', default=1e-5, help="Sequences Searching E values thread")
    group0.add_option(
        '-t', '--threads', type=int, dest='search_threads', default=8, help='Sequence searching threads. default: 8')
    group0.add_option(
        '-c', '--continue', action='store_true', dest='continued', default=False,
        help="Continue after an unfinished work. default: False")
    group1 = OptionGroup(opt, "Basic graph construction")
    group1.add_option(
        '-d', '--distance', type=str, dest='distance_defined', default='lrb',
        help='definition of the most distant hit for sequences hits; lrb: the most lower values for RBH (), '
             'arb: only select the RBH, ar:only select the RH')
    group1.add_option(
        '-w', '--weight', type=str, dest='weight_type', default='NBS',
        help='Edge weight type for community detection. Default None. Options: NBS, normalized bit-scores; '
             'BS, original bit-Scores; Sim, similarity')
    group1.add_option(
        '-m', '--community_method', type=str, dest='community_detect_method', default='rber',
        help='leidenalg community detection contains 4 algorithm based on th modularity. mvp: Modularity Vertex Partition\n'
             'rbcv: RB Configuration Vertex Partition\nrber: RBERVertexPartition\ncmp: CPM Vertex Partition\n')
    group1.add_option(
        '-a', '--network_threads', type=int, dest='network_threads', default=8,
        help='parallel running threads of network analysis')
    group1.add_option(
        '-g', '--gamma_coefficient', type=float, dest='gamma_coefficient', default=1.0,
        help='coefficient of primary resolution parameters of community detection')
    group1.add_option(
        '--so', '--species_overlap', type=int, dest='species_overlap', default=0,
        help='Species overlap threads for inference of speciation events in homologs hierarchical networks')
    opt.add_option_group(group0)
    opt.add_option_group(group1)
    options, args = opt.parse_args()
    inputDir = options.input_dir
    outputFile = options.output_graph
    continued = options.continued
    searchMethod = options.search_method
    e_values = options.e_values
    searchThreads = options.search_threads
    distance = options.distance_defined
    weightTypes = options.weight_type
    communityDetectMethod = options.community_detect_method
    GraphThreads = options.network_threads
    gammaCoefficient = options.gamma_coefficient
    speciesOverlap = options.species_overlap
    parametersDict = dict(
        inputDir=inputDir, outputFile=outputFile, searchMethod=searchMethod, e_values=e_values,
        searchThreads=searchThreads, distance=distance, weightTypes=weightTypes,
        communityDetectMethod=communityDetectMethod, GraphThreads=GraphThreads, continued=continued,
        gammaCoefficient=gammaCoefficient, speciesOverlap=speciesOverlap)
    return parametersDict


# prepare genomes sequences file methods
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


def MatricesDumpy(matrix, matrixType, name, path):
    with open(os.path.join(path, '%s%s.pic' % (matrixType, name)), 'wb') as picFile:
        pic.dump(matrix, picFile, protocol=3)


def MatricesLoad(matrixType, name, path):
    with open(os.path.join(path, '%s%s.pic' % (matrixType, name)), 'rb') as picFile:
        M = pic.load(picFile)
    return M


def GetInputGenomeInf(InputPath):
    InputGenomes = [
        GenomeFile for GenomeFile in os.listdir(InputPath) if
        GenomeFile.split('.')[-1] in ['fa', 'faa', 'fas', 'fasta']]
    return InputGenomes


def RecodeGenomeSeq(GenomesFiles, GenomePath, WorkDir):
    GenomeRecodePath = os.path.join(WorkDir, 'GenomeSeq')
    try:
        os.mkdir(GenomeRecodePath)
    except FileExistsError:
        pass
    SeqInfPair = {
        'GenomeRecodeInf': {}, 'GenesRecode': {}, 'SeqLengthInf': {}, 'SpeciesGeneNum': {}, 'GenomeToUsed': [],
        'SequencesRecode': {}}
    seqRecode = ''
    for GenomeNum, Genome in enumerate(GenomesFiles):
        SeqInfPair['GenomeRecodeInf']['G%d' % GenomeNum] = Genome
        FastaGenome = os.path.join(GenomePath, Genome)
        FastaFile = pyf.Fasta(FastaGenome)
        SeqNumbers = len(list(FastaFile.keys()))
        SeqInfPair['GenomeToUsed'].append('G%d' % GenomeNum)
        SeqInfPair['SpeciesGeneNum']['G%d' % GenomeNum] = SeqNumbers
        NewSeq = ''
        for GeneNum, SeqID in enumerate(FastaFile):
            NewSeq += '%s\n%s\n' % ('>G%d|g%d' % (GenomeNum, GeneNum), FastaFile[SeqID])
            SeqInfPair['GenesRecode']['G%d|g%d' % (GenomeNum, GeneNum)] = SeqID.split(' ')[0]
            seqRecode += '%s\tG%d|g%d\n' % (SeqID.split(' ')[0], GenomeNum, GeneNum)
            SeqInfPair['SeqLengthInf']['G%d|g%d' % (GenomeNum, GeneNum)] = len(FastaFile[SeqID][:])
            SeqInfPair['SequencesRecode'][SeqID.split(' ')[0]] = FastaFile[SeqID][:]
        OutputFile(os.path.join(GenomeRecodePath, 'G%d.fa' % GenomeNum), NewSeq)
    OutputFile('SequenceIDs.txt', seqRecode)
    [os.remove(os.path.join(GenomePath, file)) for file in os.listdir(GenomePath) if
     file.split('.')[-1] in ['flat', 'gdx']]
    return SeqInfPair, GenomeRecodePath


# homologous searching methods
def makeBlastDB(referGenomes, referGenomesPath, methods, WorkingDirectory):
    DBPath = os.path.join(WorkingDirectory, 'localDB')
    try:
        os.mkdir(DBPath)
    except FileExistsError:
        pass
    DBList = []
    for Refer in referGenomes:
        refer = os.path.join(referGenomesPath, Refer)
        db = os.path.join(DBPath, 'db%s' % Refer)
        if methods == 'blastp':
            MakeDBCmd = ' '.join(['makeblastdb', '-dbtype', 'prot', '-in', refer, '-out', db, '>', '/dev/null'])
        elif methods == 'diamond':
            MakeDBCmd = ' '.join(['diamond', 'makedb', '--in', refer, '--db', db, '--threads', '1', '>', '/dev/null'])
        else:
            print('Error! Not support searching methods')
            sys.exit(0)
        DBList.append(db)
        process = subprocess.Popen(MakeDBCmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stdin=subprocess.PIPE)
        process.wait()
    return DBList


def diamond(queue, queryGenome, queryGenomePath, BlastResultDir, DB, e_value):
    BlastQuery = os.path.join(queryGenomePath, queryGenome)
    BlastFilePath = os.path.join(
        BlastResultDir, '%s-%s.out' % (queryGenome.split('.')[0], DB.split('db')[-1].split('.')[0]))
    command = ' '.join([
        'diamond', 'blastp', '--more-sensitive', '-p', '1', '-q', BlastQuery, '-d', '%s.dmnd' % DB,
        '--evalue', str(e_value), '-f', '6', '--out', BlastFilePath, '--quiet'])
    pro = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pro.wait()
    queue.put((command, pro.returncode))


def mmseqs(queue, queryGenome, queryGenomePath, BlastResultDir, DB, e_value):
    BlastQuery = os.path.join(queryGenomePath, queryGenome)
    BlastFilePath = os.path.join(
        BlastResultDir, '%s-%s.out' % (queryGenome.split('.')[0], DB.split('/')[-1].split('.')[0]))
    tempPath = os.path.join(BlastResultDir, '%s-%s' % (queryGenome.split('.')[0], DB.split('/')[-1].split('.')[0]))
    command = ' '.join([
        'mmseqs', 'easy-search', BlastQuery, DB, BlastFilePath, tempPath, '--threads', '1', '-v', '1', '--format-mode',
        '0', '--remove-tmp-files', '-s', '7.5', '-e', str(e_value)])
    pro = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    queue.put((command, pro.returncode))


def blastp(queue, queryGenome, queryGenomePath, BlastResultDir, DB, e_value):
    BlastQuery = os.path.join(queryGenomePath, queryGenome)
    BlastFilePath = os.path.join(
        BlastResultDir, '%s-%s.out' % (queryGenome.split('.')[0], DB.split('db')[-1].split('.')[0]))
    command = ' '.join(
        ['blastp', '-outfmt', '6', '-query', BlastQuery, '-db', DB, '-evalue', str(e_value), '-out', BlastFilePath])
    pro = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pro.wait()
    queue.put((command, pro.returncode))


def RunBlastSearchParallel(queryGenomes, queryPath, searchMethod, DBList, BlastResultsPath, E_Values, thread):
    Queue = mp.Manager().Queue()
    pros = mp.Pool(processes=int(thread))
    for query in queryGenomes:
        for DB in DBList:
            if searchMethod == 'blastp':
                pros.apply_async(func=blastp, args=(Queue, query, queryPath, BlastResultsPath, DB, E_Values))
            elif searchMethod == 'diamond':
                pros.apply_async(func=diamond, args=(Queue, query, queryPath, BlastResultsPath, DB, E_Values))
            elif searchMethod == 'mmseqs':
                pros.apply_async(func=mmseqs, args=(Queue, query, queryPath, BlastResultsPath, DB, E_Values))
            else:
                pass
    PrintParallelBar(Queue, 'Homologs Searching:', len(queryGenomes))
    pros.close()
    pros.join()


def SequencesSearchBlast(InputGenomePath, searchingMethod, WorkingDirectory, EValues, searchThreads):
    InputGenomes = os.listdir(InputGenomePath)
    if searchingMethod == 'mmseqs':
        DBs = [os.path.join(InputGenomePath, genomes) for genomes in os.listdir(InputGenomePath)]
    else:
        DBs = makeBlastDB(InputGenomes, InputGenomePath, searchingMethod, WorkingDirectory)
    BlastResultsPath = os.path.join(WorkingDirectory, 'BlastResults')
    try:
        os.mkdir(BlastResultsPath)
    except FileExistsError:
        pass
    RunBlastSearchParallel(
        InputGenomes, InputGenomePath, searchingMethod, DBs, BlastResultsPath, EValues, searchThreads)
    return BlastResultsPath


# Get original matrices between genomes pair
def ReadBlastResults(BlastFileName, SeqInformation):
    """
    queryID: query sequences id
    referID: refer sequences id
    queryLength: length of query sequences
    referLength: length of refer sequences
    LengthProduct: product of  length of query and refer sequences
    :param SeqInformation: sequences information: {seqID: length of sequences}
    :param BlastFileName: file name list of all-against-all blast results
    :return: BlastHitList: handled blast out list
    """
    SimilarityPair = {}
    bitScoresMatrices = {}
    QueryGenome = BlastFileName.split('/')[-1].split('.')[0].split('-')[0]
    referGenome = BlastFileName.split('/')[-1].split('.')[0].split('-')[1]
    LengthSparseMatrix = sparse.lil_matrix(
        (SeqInformation['SpeciesGeneNum'][QueryGenome], SeqInformation['SpeciesGeneNum'][referGenome]))
    BitScoreSparseMatrix = LengthSparseMatrix.copy()
    BlastHitLines = ReadFile(BlastFileName)
    for hitLine in BlastHitLines:
        try:
            queryID = hitLine.strip('\n').split('\t')[0]
            referID = hitLine.strip('\n').split('\t')[1]
            mathLength = int(hitLine.strip('\n').split('\t')[3])
            queryGeneIndex = int(queryID.split('|')[-1].split('g')[-1])
            referGeneIndex = int(referID.split('|')[-1].split('g')[-1])
            queryLength = SeqInformation['SeqLengthInf'][queryID]
            referLength = SeqInformation['SeqLengthInf'][referID]
            LengthProduct = queryLength * referLength
            queryCover = (mathLength / queryLength) * 100
            referCover = (mathLength / referLength) * 100
            Similarity = float(hitLine.strip('\n').split('\t')[2])
            bitscores = float(hitLine.strip('\n').split('\t')[-1])
            coverage = queryCover if queryCover < referCover else referCover
            if queryID != referID:
                LengthSparseMatrix[queryGeneIndex, referGeneIndex] = LengthProduct
                BitScoreSparseMatrix[queryGeneIndex, referGeneIndex] = bitscores
                SimilarityPair['%s-%s' % (queryID, referID)] = Similarity
                bitScoresMatrices['%s-%s' % (queryID, referID)] = bitscores
        except ValueError:
            print(hitLine)
    return SimilarityPair, bitScoresMatrices, LengthSparseMatrix, BitScoreSparseMatrix


class BSN:
    @staticmethod
    def RetainTopData(LengthMatrix, BitScoresMatrix):
        """
        :param LengthMatrix: Length product of query and hit
        :param BitScoresMatrix: BitScores of query and hit
        :return TopLength: sequences length in top 5% bit scores bins from different length bins:
        :return TopBitScores BitScores in top 5% bit scores bins from different length bins:
        """
        LengthArray = [
            LengthMatrix[rowIndex, colIndex] for rowIndex, colIndex in
            zip(LengthMatrix.nonzero()[0], LengthMatrix.nonzero()[1])]
        BitScoresArray = [
            BitScoresMatrix[rowIndex, colIndex]
            for rowIndex, colIndex in zip(BitScoresMatrix.nonzero()[0], BitScoresMatrix.nonzero()[1])]
        LengthSorted = sorted([(length, bit) for length, bit in zip(LengthArray, BitScoresArray)], key=lambda X: X[0])
        HitNumbers = len(LengthSorted)
        if HitNumbers < 100:
            LengthBinRaw = [element[0] for element in LengthSorted]
            BitBinRaw = [element[1] for element in LengthSorted]
            return LengthBinRaw, BitBinRaw
        LengthTop = []
        BitTop = []
        scale = 1000 if HitNumbers > 5000 else (200 if HitNumbers > 1000 else 20)
        for BinsIndex in range(0, len(LengthSorted), scale):
            bitBin = [binLen[-1] for binLen in LengthSorted[BinsIndex:BinsIndex + scale + 1]]
            cutoff = np.percentile(np.array(bitBin), 95)
            for BinElement in LengthSorted[BinsIndex:BinsIndex + scale + 1]:
                if BinElement[-1] >= cutoff:
                    LengthTop.append(BinElement[0])
                    BitTop.append(BinElement[1])
        Top5Length = np.array(LengthTop)
        Top5BitScores = np.array(BitTop)
        return Top5Length, Top5BitScores

    @staticmethod
    def FitnessFunction(x, a, b):
        y = a * np.log10(x) + b
        return y

    @staticmethod
    def GetFitnessPara(xData, yData):
        A, B = curve_fit(BSN.FitnessFunction, xData, np.log10(yData))[0]
        return A, B

    @staticmethod
    def NormalizedFunction(LengthProduct, RawBitScores, ParaA, ParaB):
        """
        NormalizedBitscores: normalized bit-scores
        :param RawBitScores: raw bit-scores
        :param LengthProduct: length product of hit and query sequences
        :param RawBitScores: raw bit-scores
        :param ParaA: fitness function coefficient a
        :param ParaB: fitness function constant b
        :return nbs:normalized bit-scores
        """
        NormalizedBitScores = RawBitScores / ((10 ** ParaB) * (LengthProduct.power(ParaA)))
        return NormalizedBitScores


# Get connected matrices
class MatrixHandle:
    def __init__(self, MatrixDirectory, SpeciesInf, iGenome):
        self.MatrixDirectory = MatrixDirectory
        self.SpeciesInf = SpeciesInf
        self.genomeI = iGenome

    def GetSpecieNBMatrix(self, BlastHitDir):
        """
        normalized bit-scores
        :param BlastHitDir: Blast results files directory
        """
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            BlastFile = '%s-%s.out' % (self.genomeI, genomeJ)
            matrixName = '%s-%s' % (self.genomeI, genomeJ)
            BlastFileName = os.path.join(BlastHitDir, BlastFile)
            Similarities, bitScoresDic, LengthMatrix, BitScoreMatrix = ReadBlastResults(BlastFileName, self.SpeciesInf)
            MatricesDumpy(Similarities, 'similarity', matrixName, self.MatrixDirectory)
            MatricesDumpy(bitScoresDic, 'bitScores', matrixName, self.MatrixDirectory)
            if BitScoreMatrix.count_nonzero() > 1:
                LengthTopArray, BitTopArray = BSN.RetainTopData(LengthMatrix, BitScoreMatrix)
                paraA, ParaB = BSN.GetFitnessPara(LengthTopArray, BitTopArray)
                NormalBitMatrix = BSN.NormalizedFunction(LengthMatrix, BitScoreMatrix, paraA, ParaB)
                NormalBitMatrix[np.isnan(NormalBitMatrix)] = 0
                MatricesDumpy(sparse.lil_matrix(NormalBitMatrix), 'NB', matrixName, self.MatrixDirectory)
                MatricesDumpy(LengthMatrix.tolil(), 'length', matrixName, self.MatrixDirectory)
            else:
                print('Warning! There are no hits in all-to-all blast except the same genes:%s-%s' % (
                    self.SpeciesInf['GenomeRecodeInf'][self.genomeI], self.SpeciesInf['GenomeRecodeInf'][genomeJ]))
                NormalBitMatrix = sparse.lil_matrix(
                    (self.SpeciesInf['SpeciesGeneNum'][self.genomeI], self.SpeciesInf['SpeciesGeneNum'][genomeJ]))
                LengthMatrix = NormalBitMatrix.copy()
                MatricesDumpy(sparse.lil_matrix(NormalBitMatrix), 'NB', matrixName, self.MatrixDirectory)
                MatricesDumpy(LengthMatrix.tolil(), 'length', matrixName, self.MatrixDirectory)

    def GetBestHit(self, accepted=1e-3):
        """
        BestHitI: Best Hit index for each genes from Genome I in original normalized matrix
        Matrices: normalized matrix for each species pairs
        BitHitBestIndexList: best hit for each genome pair
        """
        BestHitI = -1 * np.ones(self.SpeciesInf['SpeciesGeneNum'][self.genomeI])
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            if self.genomeI == genomeJ:
                continue
            rowMatrix = MatricesLoad('NB', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            BestHitIndexList = []
            BestHitRowList = []
            for rowNum in range(self.SpeciesInf['SpeciesGeneNum'][self.genomeI]):
                bitSoresObject = rowMatrix.getrowview(rowNum)
                if bitSoresObject.nnz > 0:
                    MaxValues = max(bitSoresObject.data[0])
                    BestHitI[rowNum] = MaxValues if MaxValues > BestHitI[rowNum] else BestHitI[rowNum]
                    BestHitsColIndex = [
                        colIndex for colIndex, bitValue in zip(bitSoresObject.rows[0], bitSoresObject.data[0])
                        if bitValue > MaxValues - accepted]
                    BestHitIndexList.extend(BestHitsColIndex)
                    BestHitRowList.extend(np.full(len(BestHitsColIndex), rowNum, dtype=np.dtype(int)))
            BestHitMatrixForDump = sparse.csr_matrix(
                (np.ones(len(BestHitRowList)), (BestHitRowList, BestHitIndexList)),
                shape=rowMatrix.get_shape())
            MatricesDumpy(BestHitMatrixForDump, 'BH', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
        # paralogs
        rowMatrix = MatricesLoad('NB', '%s-%s' % (self.genomeI, self.genomeI), self.MatrixDirectory)
        BestHitIndexList = []
        BestHitRowList = []
        for rowNum in range(self.SpeciesInf['SpeciesGeneNum'][self.genomeI]):
            bitSoresObject = rowMatrix.getrowview(rowNum)
            if bitSoresObject.nnz > 0:
                BestHitsColIndex = [
                    colIndex for colIndex, bitValue in zip(bitSoresObject.rows[0], bitSoresObject.data[0])
                    if bitValue > BestHitI[rowNum] - accepted]
                BestHitIndexList.extend(BestHitsColIndex)
                BestHitRowList.extend(np.full(len(BestHitsColIndex), rowNum, dtype=np.dtype(int)))
        BestHitMatrixForDump = sparse.csr_matrix(
            (np.ones(len(BestHitRowList)), (BestHitRowList, BestHitIndexList)),
            shape=rowMatrix.get_shape())
        MatricesDumpy(BestHitMatrixForDump, 'BH', '%s-%s' % (self.genomeI, self.genomeI), self.MatrixDirectory)

    def RBHMatrix(self):
        """
        calculate reciprocal best hit for each genome pairs
        BHs: best hit matrix
        RBHs: reciprocal best hit matrices for  each genome pairs
        """
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            if genomeJ == self.genomeI:
                RBHs = sparse.csr_matrix(
                    (self.SpeciesInf['SpeciesGeneNum'][self.genomeI], self.SpeciesInf['SpeciesGeneNum'][genomeJ]))
            else:
                BHs = MatricesLoad('BH', '%s-%s' % (self.genomeI, genomeJ), MatrixDir)
                RBH = MatricesLoad('BH', '%s-%s' % (genomeJ, self.genomeI), MatrixDir)
                RBHs = BHs.multiply(RBH.transpose())
            MatricesDumpy(RBHs, 'RBHs', '%s-%s' % (self.genomeI, genomeJ), MatrixDir)

    def SortedConnected(self):
        """
        RBHM: reciprocal best hits matrix
        NBMatrices: Normalized Bit-scores matrices
        SeqInformation: sequences information
        ConnectedMatrixPair : connected genes pair for each genome pairs
        """
        GenomeIRBMin = 1e9 * np.ones(self.SpeciesInf['SpeciesGeneNum'][self.genomeI])
        bestHit = np.zeros(self.SpeciesInf['SpeciesGeneNum'][self.genomeI])
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            if genomeJ == self.genomeI:
                continue
            NBMetricIJ = MatricesLoad('NB', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            NBMetricIJCsr = NBMetricIJ.tocsr()
            RBHMIJ = MatricesLoad('RBHs', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            for rowNum in range(self.SpeciesInf['SpeciesGeneNum'][self.genomeI]):
                if NBMetricIJ.getrowview(rowNum).nnz > 0:
                    bestHit[rowNum] = max(bestHit[rowNum], max(NBMetricIJ.getrowview(rowNum).data[0]))
                if RBHMIJ[rowNum].nnz > 0:
                    GenomeIRBMin[rowNum] = min(
                        min(NBMetricIJCsr[rowNum, RBHMIJ[rowNum].indices].data), GenomeIRBMin[rowNum])
        Indices = GenomeIRBMin > 1e8
        GenomeIRBMin[Indices] = bestHit[Indices] + 1e-6
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            NBMetricIJ = MatricesLoad('NB', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            NBMetricIJCsr = NBMetricIJ.tocsr()
            RetainedHitIJ = []
            RetainedHitRow = []
            RetainedHitCol = []
            for rowNum in range(self.SpeciesInf['SpeciesGeneNum'][self.genomeI]):
                if NBMetricIJCsr[rowNum].nnz > 0:
                    for values, index in zip(NBMetricIJCsr[rowNum].data, NBMetricIJCsr[rowNum].indices):
                        if values >= GenomeIRBMin[rowNum]:
                            RetainedHitCol.append(index)
                            RetainedHitRow.append(rowNum)
                            RetainedHitIJ.append(values)
            ConnectedMatrix = sparse.csr_matrix(
                (RetainedHitIJ, (RetainedHitRow, RetainedHitCol)), shape=NBMetricIJCsr.get_shape())
            MatricesDumpy(ConnectedMatrix, 'CM', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)

    def SortedConnectedRBHs(self):
        """
        Only Selected reciprocal best hits for SSN construction
        RBHMIJ: Reciprocal best hits matrix from Species I and Species J pairs
        RetainedHitIJ: normalized bit-scores values of retained hit
        RetainedHitRow: row index of retained hit
        RetainedHitCol: col index of retained hit
        """
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            RBHMIJ = MatricesLoad('RBHs', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            NBMetricIJ = MatricesLoad('NB', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            RBHMIJCsr = RBHMIJ.tocsr()
            RetainedHitIJ = []
            RetainedHitRow = []
            RetainedHitCol = []
            for rowNum in range(self.SpeciesInf['SpeciesGeneNum'][self.genomeI]):
                if self.genomeI == genomeJ:
                    for index in RBHMIJCsr[rowNum].indices:
                        RetainedHitCol.append(index)
                        RetainedHitRow.append(rowNum)
                        values = NBMetricIJ[rowNum, index]
                        RetainedHitIJ.append(values)
                else:
                    if RBHMIJCsr[rowNum].nnz > 0:
                        for index in RBHMIJCsr[rowNum].indices:
                            RetainedHitCol.append(index)
                            RetainedHitRow.append(rowNum)
                            values = NBMetricIJ[rowNum, index]
                            RetainedHitIJ.append(values)
                ConnectedMatrix = sparse.csr_matrix(
                    (RetainedHitIJ, (RetainedHitRow, RetainedHitCol)), shape=NBMetricIJ.get_shape())
                MatricesDumpy(ConnectedMatrix, 'CM', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)

    def SortedConnectedRHs(self):
        """Only Selected reciprocal hits for SSN construction
        NBMetricIJ: Normalized bit-scores matrix of pair between species I and species J
        RetainedHitIJ: normalized bit-scores values of retained hit
        RetainedHitRow: row index of retained hit
        RetainedHitCol: col index of retained hit
        """
        for genomeJ in self.SpeciesInf['GenomeToUsed']:
            NBMetricIJ = MatricesLoad('NB', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            NBMetricIJCsr = NBMetricIJ.tocsr()
            NBMetricJI = MatricesLoad('NB', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)
            NBMetricJICsr = NBMetricJI.toscr()
            RetainedHitIJ = []
            RetainedHitRow = []
            RetainedHitCol = []
            for rowNum in range(self.SpeciesInf['SpeciesGeneNum'][self.genomeI]):
                if self.genomeI == genomeJ:
                    for index in NBMetricIJCsr[rowNum].indices:
                        RetainedHitCol.append(index)
                        RetainedHitRow.append(rowNum)
                        values = NBMetricIJ[rowNum, index]
                        RetainedHitIJ.append(values)
                else:
                    for values, index in zip(NBMetricIJCsr[rowNum].data, NBMetricIJCsr[rowNum].indices):
                        if NBMetricJICsr[rowNum, index] > 0:
                            RetainedHitCol.append(index)
                            RetainedHitRow.append(rowNum)
                            RetainedHitIJ.append(values)
                ConnectedMatrix = sparse.csr_matrix(
                    (RetainedHitIJ, (RetainedHitRow, RetainedHitCol)), shape=NBMetricIJ.get_shape())
                MatricesDumpy(ConnectedMatrix, 'CM', '%s-%s' % (self.genomeI, genomeJ), self.MatrixDirectory)


def GetBHMatrixAll(queue, SpeciesI, MostD, searchResultPath, recodeInf, matrixDirectory):
    MatricesObject = MatrixHandle(matrixDirectory, recodeInf, SpeciesI)
    MatricesObject.GetSpecieNBMatrix(searchResultPath)
    if MostD == 'lrb' or MostD == 'arb':
        MatricesObject.GetBestHit()
    else:
        pass
    queue.put((SpeciesI, 0))


def GetBHMatrixParallel(ResultPath, mostDistances, reInf, matDirectory, matrixThreads):
    MyQueue = mp.Manager().Queue()
    P = mp.Pool(processes=int(matrixThreads))
    for iSpecies in reInf['GenomeToUsed']:
        P.apply_async(func=GetBHMatrixAll, args=(MyQueue, iSpecies, mostDistances, ResultPath, reInf, matDirectory))
    PrintParallelBar(MyQueue, 'Get Matrices', len(reInf['GenomeToUsed']))
    P.close()
    P.join()


def PrintParallelBar(queue, parallelType, taskNum):
    progressbar_widgets_set = [
        parallelType, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(), ' | ',
        progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=taskNum)
    bar.start()
    doneNum = 0
    while True:
        cmd, completeNum = queue.get()
        if completeNum == 0:
            pass
        else:
            print('%s Error!' % cmd)
        doneNum += 1
        bar.update(doneNum)
        if doneNum >= taskNum:
            break
    bar.finish()


def GetConnectedMatrixSpeciesI(queue, SpeciesI, MostD, recodeInf, matrixDirectory):
    MatricesAbout = MatrixHandle(matrixDirectory, recodeInf, SpeciesI)
    if MostD == 'lrb':
        MatricesAbout.RBHMatrix()
        MatricesAbout.SortedConnected()
    elif MostD == 'arb':
        MatricesAbout.SortedConnectedRBHs()
    elif MostD == 'ar':
        MatricesAbout.SortedConnectedRHs()
    else:
        print('Error: %s is not defined method for the most distances' % MostD)
        sys.exit(0)
    queue.put((SpeciesI, 0))


def GetConnectedMatrixParallel(inf, mostDistance, matDir, nt):
    Queue = mp.Manager().Queue()
    pro = mp.Pool(processes=int(nt))
    for iGenome in inf['GenomeToUsed']:
        pro.apply_async(func=GetConnectedMatrixSpeciesI, args=(Queue, iGenome, mostDistance, inf, matDir))
    PrintParallelBar(Queue, 'Construction of Connected Matrices:', len(inf['GenomeToUsed']))
    pro.close()
    pro.join()


def BuildSSN_parallel(MatrixDirectory, speciesInf, buildThreads):
    MyQueue = mp.Manager().Queue()
    pro = mp.Pool(processes=int(buildThreads))
    for i in range(len(speciesInf['GenomeToUsed'])):
        for j in range(i, len(speciesInf['GenomeToUsed'])):
            iSpecies, jSpecies = speciesInf['GenomeToUsed'][i], speciesInf['GenomeToUsed'][j]
            pro.apply_async(
                func=GetConnections, args=(MyQueue, iSpecies, jSpecies, MatrixDirectory, speciesInf))
    pairNum = comb(len(speciesInf['GenomeToUsed']), 2) + len(speciesInf['GenomeToUsed'])
    ssn = UnionGraphs(MyQueue, pairNum)
    pro.close()
    pro.join()
    return ssn


def GetConnections(queues, SpeciesI, SpeciesJ, MatrixDirectory, SpeciesInf):
    """
    Getting edges and its attribution based on genes connected matrices
    :param SpeciesInf: sequences information record
    :return: graphs of species pairwise
    """
    rawGraph = igraph.Graph()
    connectMatrix = MatricesLoad('CM', '%s-%s' % (SpeciesI, SpeciesJ), MatrixDirectory)
    SP = MatricesLoad('similarity', '%s-%s' % (SpeciesI, SpeciesJ), MatrixDirectory)
    BP = MatricesLoad('bitScores', '%s-%s' % (SpeciesI, SpeciesJ), MatrixDirectory)
    NodeName = []
    EdgeList = []
    SimList = []
    BSList = []
    NBSList = []
    NBSRList = []
    SimRList = []
    BSRList = []
    if SpeciesI != SpeciesJ:
        connectMatrixR = MatricesLoad('CM', '%s-%s' % (SpeciesJ, SpeciesI), MatrixDirectory).transpose()
        SPR = MatricesLoad('similarity', '%s-%s' % (SpeciesJ, SpeciesI), MatrixDirectory)
        BPR = MatricesLoad('bitScores', '%s-%s' % (SpeciesJ, SpeciesI), MatrixDirectory)
        connectMatrixAdd = connectMatrix + connectMatrixR
        for rowNum in range(SpeciesInf['SpeciesGeneNum'][SpeciesI]):
            sourceNode = '%s|g%d' % (SpeciesI, rowNum)
            if sourceNode not in NodeName:
                NodeName.append(sourceNode)
            if connectMatrixAdd[rowNum].nnz > 0:
                for Index in connectMatrixAdd[rowNum].indices:
                    NBS = connectMatrix[rowNum, Index]
                    NBS1 = connectMatrixR[rowNum, Index]
                    if NBS > 0 and NBS1 > 0:
                        NBSValues = np.float64(NBS).item()
                        NBSR = np.float64(NBS1).item()
                        Similarity = SP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        SimilarityR = SPR['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                        BitScores = BP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        BitScoresR = BPR['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                    elif NBS > 0 and NBS1 == 0:
                        NBSValues = np.float64(NBS).item()
                        NBSR = 0
                        Similarity = SP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        SimilarityR = 0
                        BitScores = BP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        BitScoresR = 0
                    elif NBS == 0 and NBS1 > 0:
                        NBSValues = np.float64(NBS1).item()
                        NBSR = 0
                        Similarity = SPR['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                        SimilarityR = 0
                        BitScores = BPR['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                        BitScoresR = 0
                    else:
                        NBSValues = 0
                        NBSR = 0
                        Similarity = 0
                        SimilarityR = 0
                        BitScores = 0
                        BitScoresR = 0
                    targetNode = '%s|g%d' % (SpeciesJ, Index)
                    if targetNode not in NodeName:
                        NodeName.append(targetNode)
                    EdgeList.append([sourceNode, targetNode])
                    SimList.append(Similarity)
                    SimRList.append(SimilarityR)
                    BSList.append(BitScores)
                    BSRList.append(BitScoresR)
                    NBSList.append(NBSValues)
                    NBSRList.append(NBSR)
    else:
        for rowNum in range(SpeciesInf['SpeciesGeneNum'][SpeciesI]):
            sourceNode = '%s|g%d' % (SpeciesI, rowNum)
            if sourceNode not in NodeName:
                NodeName.append(sourceNode)
            if connectMatrix[rowNum].nnz > 0:
                for Index in connectMatrix[rowNum].indices:
                    targetNode = '%s|g%d' % (SpeciesJ, Index)
                    if [targetNode, sourceNode] not in EdgeList:
                        if connectMatrix[Index, rowNum] > 0:
                            NBS = connectMatrix[rowNum, Index]
                            NBS1 = connectMatrix[Index, rowNum]
                            Similarity = SP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            SimilarityR = SP['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                            BitScores = BP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            BitScoresR = BP['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                        else:
                            NBS = connectMatrix[rowNum, Index]
                            NBS1 = 0
                            Similarity = SP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            SimilarityR = 0
                            BitScores = BP['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            BitScoresR = 0
                        if targetNode not in NodeName:
                            NodeName.append(targetNode)
                        EdgeList.append([sourceNode, targetNode])
                        SimList.append(Similarity)
                        SimRList.append(SimilarityR)
                        BSList.append(BitScores)
                        BSRList.append(BitScoresR)
                        NBSList.append(NBS)
                        NBSRList.append(NBS1)
    rawGraph.add_vertices(NodeName)
    attributesDic = dict(Sim=SimList, BS=BSList, NBS=NBSList, NBSRev=NBSRList, BSRev=BSRList, SimRev=SimRList)
    rawGraph.add_edges(EdgeList, attributes=attributesDic)
    queues.put(rawGraph)


def UnionGraphs(graphLists, taskNum):
    """
    Union graphs of genome pairwise
    :param graphLists: graph queue
    :param taskNum: task numbers
    :return: sequences similarity network
    """
    DoneTask = 0
    progressbar_widgets_set = [
        'Construction of SSN:', progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(), ' | ',
        progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=taskNum)
    bar.start()
    graphList = []
    while True:
        graph = graphLists.get()
        if graph.is_named():
            graphList.append(graph)
        DoneTask += 1
        bar.update(DoneTask)
        if DoneTask >= taskNum:
            break
    bar.finish()
    ssn = igraph.union(graphList, byname=True)
    return ssn


def MappingGammaForCC(coefficientG, lengthCC):
    coefficientG = (10 ** (int(np.log10(lengthCC) - 4))) * coefficientG
    index1 = float(np.log10(lengthCC))
    gamma = float(coefficientG) * (10 ** (-index1 + 1))
    return gamma


# HHN build Methods
class HHN(igraph.Graph):
    def __init__(self, ssn=None, hhn=None, **attr):
        super(HHN, self).__init__(**attr)
        self.ssn = ssn
        self.hhn = hhn

    def splitConnectedComponents(self):
        connectedComponents = [('%d' % num, cc) for num, cc in enumerate(self.ssn.components()) if len(cc) >= 2]
        return connectedComponents

    def SpecificBoard(self, Nodes, detectionMethods, edgeAttribute, gc, startNum, stopNum):
        communityGraph = self.ssn.subgraph(Nodes)
        left = 0
        right = MappingGammaForCC(gc, len(Nodes))
        methodOptional = dict(
            mvp=la.ModularityVertexPartition, rbcv=la.RBConfigurationVertexPartition,
            rber=la.RBERVertexPartition, cmp=la.CPMVertexPartition)
        n = 1
        p = None
        while True:
            partition = la.find_partition(
                communityGraph, partition_type=methodOptional[detectionMethods],
                n_iterations=10, weights=edgeAttribute, resolution_parameter=right)
            if len(partition.subgraphs()) > stopNum:
                right = right
                break
            elif startNum <= len(partition.subgraphs()) <= stopNum:
                p = partition.subgraphs()
                break
            else:
                left = right
                right = right + (right / 2)
            n += 1
        return left, right, n, p, partition.q

    def BipartiteGraphs(self, communityNode, detectionMethods, edgeAttribute, gc, startNum, stopNum):
        communityGraph = self.ssn.subgraph(communityNode)
        methodOptional = dict(
            mvp=la.ModularityVertexPartition, rbcv=la.RBConfigurationVertexPartition,
            rber=la.RBERVertexPartition, cmp=la.CPMVertexPartition)
        if len(communityNode) >= 1000:
            left, right, dn, p, quality = self.SpecificBoard(communityNode, detectionMethods, edgeAttribute, gc,
                                                             startNum,
                                                             stopNum)
        else:
            left, right, dn, p = 0, 1, 1, None
        if p is None:
            n = 1
            ResolutionParameter = (left + right) / 2
            partition = la.find_partition(
                communityGraph, partition_type=methodOptional[detectionMethods],
                n_iterations=10, weights=edgeAttribute, resolution_parameter=ResolutionParameter)
            while n <= 1000:
                partition = la.find_partition(
                    communityGraph, partition_type=methodOptional[detectionMethods],
                    n_iterations=10, weights=edgeAttribute, resolution_parameter=ResolutionParameter)
                if len(partition.subgraphs()) < startNum:
                    left = ResolutionParameter
                    ResolutionParameter = (left + right) / 2
                elif len(partition.subgraphs()) > stopNum:
                    right = ResolutionParameter
                    ResolutionParameter = (left + right) / 2
                else:
                    break
                n += 1
            partitionG = partition.subgraphs()
            quality = partition.q
        else:
            partitionG = p
            quality = 'None'
        return partitionG, quality

    def RunCommunityDetection(self, connected, modularise, detectedMethod, weightType, gc):
        """
        implement community detection
        :param connected: connected list
        :param detectedMethod: communities detection methods
        :param weightType: edge weight type for community detection
        :param modularise: connected components contains community on corresponding iteration
        :return: communities from current iteration
        """
        communitiesList = []
        for ID, modular in modularise:
            sn = ID
            snAttribution, modularGenomeNum = GetAttribution(self.ssn.subgraph(modular), True)
            if modularGenomeNum > 1 and len(modular) >= 10000:
                stopNum = 50
                startNum = 20
                Communities, Quality = self.BipartiteGraphs(modular, detectedMethod, weightType, gc, startNum,
                                                            stopNum)
                if len(Communities) >= 2:
                    tns = ''
                    for num, community in enumerate(Communities):
                        CCsCommunities = [int(node) for node in community.vs['id']]
                        tn = '%s-%d' % (ID, num)
                        tns += tn + ' '
                        communitiesList.append((tn, CCsCommunities))
                    connectedInf = ','.join([sn, tns, snAttribution, str(Quality)])
                else:
                    connectedInf = ','.join([sn, '+', snAttribution, '0'])
            elif modularGenomeNum > 1 and 2 < len(modular) < 10000:
                stopNum = 2
                startNum = 2
                Communities, Quality = self.BipartiteGraphs(modular, detectedMethod, weightType, gc, startNum,
                                                            stopNum)
                if len(Communities) == 2:
                    tns = ''
                    for num, community in enumerate(Communities):
                        CCsCommunities = [int(node) for node in community.vs['id']]
                        tn = '%s-%d' % (ID, num)
                        tns += tn + ' '
                        communitiesList.append((tn, CCsCommunities))
                    connectedInf = ','.join([sn, tns, snAttribution, str(Quality)])
                else:
                    connectedInf = ','.join([sn, '+', snAttribution, '0'])
            elif modularGenomeNum > 1 and len(modular) == 2:
                tns = '%s-0 %s-1' % (ID, ID)
                connectedInf = ','.join([sn, tns, snAttribution, '0'])
                communitiesList.append(('%s-0' % ID, [modular[0]]))
                communitiesList.append(('%s-1' % ID, [modular[1]]))
            else:
                connectedInf = ','.join([sn, '+', snAttribution, '0'])
            connected.append(connectedInf)
        return communitiesList

    def GetEvolutionEvents(self, overlapThreads):
        """
        Mark evolution events for nodes on network
        :param overlapThreads: numbers of species overlap
        :return: network
        """
        SelectedVertex = self.hhn.vs.select(_degree_gt=1)
        for vertex in SelectedVertex:
            sgs = set(vertex['genomeIDs'].split(' '))
            vertexGenesNum = vertex['genesNum']
            gcs = []
            for node in vertex.neighbors():
                if int(node['genesNum']) < int(vertexGenesNum):
                    gcs.append(set(node['genomeIDs'].strip().split(' ')))
            if len(gcs) > 2:
                vertex['Event'] = 'III-3'
            else:
                cs1, cs2 = gcs[0], gcs[1]
                if cs1 & cs2 == sgs:
                    vertex['Event'] = 'II'
                elif cs1 & cs2 == set():
                    vertex['Event'] = 'I'
                else:
                    if len(cs1 & cs2) <= overlapThreads and len(cs1 & cs2) < len(sgs):
                        vertex['Event'] = 'I'
                    else:
                        if cs1 == sgs or cs2 == sgs:
                            vertex['Event'] = 'III-1'
                        else:
                            vertex['Event'] = 'III-2'
        SelectedVertex0 = self.hhn.vs.select(_degree=0)
        for vertex in SelectedVertex0:
            sourceNodeGenome = vertex['genomesNum']
            sourceNodeGenes = vertex['genesNum']
            if sourceNodeGenome == sourceNodeGenes:
                vertex['Event'] = 'I'
            else:
                pass
        return self.hhn

    def ExtractOG(self, GenomeNum):
        adLists = []
        if GenomeNum == 1:
            ogNodes = self.hhn.vs.select(Event=None)
        elif GenomeNum == 0:
            ogNodes = self.ssn.vs.select(_degree=0)
        else:
            ogNodes = self.hhn.vs.select(Event='I', genomesNum=GenomeNum)
        nodesDeleted = []
        OGInf = {}
        for num, ogNode in enumerate(ogNodes):
            if GenomeNum == 1:
                genesIDList = ogNode['geneIDs'].split(' ')
                nodeCC = None
            elif GenomeNum == 0:
                genesIDList = ogNode['name'].split(' ')
                nodeCC = None
            else:
                genesIDList, deletedVertex = GetGenesIDs(ogNode)
                nodeCC = [n['name'] for n in deletedVertex]
                nodeCC.append(ogNode['name'])
                if OneNodeCoalescence(ogNode):
                    parent, adjacent = OneNodeCoalescence(ogNode)[0]
                    genesIDList.extend(adjacent['geneIDs'].strip().split(' '))
                    deletedVertex.extend([ogNode, adjacent])
                    nodeCC.extend([parent['name'], adjacent['name']])
                    ogNode = parent
                else:
                    ogNode = ogNode
                nodesDeleted.extend(deletedVertex)
            genesIDs = ' '.join([i for i in set(genesIDList)])
            OGInf[ogNode['name']] = (GenomeNum, len(genesIDList), genesIDs, nodeCC)
        if GenomeNum > 1:
            self.hhn.delete_vertices(nodesDeleted)
        return OGInf, self.hhn, adLists


def OneNodeCoalescence(NodeObject):
    NoneNode = []
    for neighbor in NodeObject.neighbors():
        if neighbor['genesNum'] > NodeObject['genesNum']:
            parent = neighbor
            for child in parent.neighbors():
                if child['genesNum'] < parent['genesNum'] and child['name'] != NodeObject['name']:
                    neighborhood = child
                    if neighborhood['Event'] is None or neighborhood['Event'] == 'None':
                        NoneNode.append((parent, neighborhood))
    if len(NoneNode) > 1:
        return []
    else:
        return NoneNode


def GetConnectedList(Queues, taskNum, turns):
    taskDoneNum = 0
    progressbar_widgets_set = [
        'Split connected components %d:' % turns, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(),
        ' | ',
        progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=taskNum)
    secondQueueList = []
    NodeList = []
    EdgeList = []
    genesIDsList = []
    genomeIDsList = []
    genesNumList = []
    genomeNumList = []
    modularityList = []
    bar.start()
    while True:
        connects, secondQueue = Queues.get()
        if connects:
            for connected in connects:
                sn, tnsGenes, GenesIDs, GenesNum, genomes, genomeNum, modularity = connected.strip().split(',')
                if sn not in NodeList:
                    NodeList.append(sn)
                    genesIDsList.append(GenesIDs)
                    genomeIDsList.append(genomes)
                    genesNumList.append(int(GenesNum))
                    genomeNumList.append(int(genomeNum))
                    modularityList.append(modularity)
                if '-' in tnsGenes:
                    for tn in tnsGenes.strip().split(' '):
                        if [sn, tn] not in EdgeList:
                            EdgeList.append([sn, tn])
        if secondQueue:
            secondQueueList.extend(secondQueue)
        taskDoneNum += 1
        bar.update(taskDoneNum)
        if taskDoneNum >= taskNum:
            break
    connectedList = [
        NodeList, EdgeList, genesIDsList, genomeIDsList, genesNumList, genomeNumList, modularityList]
    bar.finish()
    return connectedList, secondQueueList


def SplitTasks1(CCLists):
    TimeUsedU = []
    TimeUsedS = []
    TimeUsedM = []
    TimeUsedF = []
    CCGroups = []
    for ids, cc in CCLists:
        if len(cc) >= 10 ** 5:
            TimeUsedU.append((ids, cc))
        elif 10 ** 5 > len(cc) >= 10 ** 4:
            TimeUsedS.append((ids, cc))
        elif 10 ** 4 > len(cc) >= 10 ** 3:
            TimeUsedM.append((ids, cc))
        elif 10 ** 3 > len(cc) >= 0:
            TimeUsedF.append((ids, cc))
        else:
            pass
    if TimeUsedU:
        random.shuffle(TimeUsedU)
        [CCGroups.append(TimeUsedU[i:i + 1]) for i in range(0, len(TimeUsedU), 1)]
    if TimeUsedS:
        random.shuffle(TimeUsedS)
        [CCGroups.append(TimeUsedS[i:i + 15]) for i in range(0, len(TimeUsedS), 15)]
    if TimeUsedM:
        random.shuffle(TimeUsedM)
        [CCGroups.append(TimeUsedM[i:i + 150]) for i in range(0, len(TimeUsedM), 150)]
    if TimeUsedF:
        random.shuffle(TimeUsedF)
        [CCGroups.append(TimeUsedF[i:i + 1500]) for i in range(0, len(TimeUsedF), 1500)]
    return CCGroups, None


def SplitTasks2(CCLists):
    TimeUsedM = []
    TimeUsedF = []
    CCGroups = []
    for ids, cc in CCLists:
        if 300 <= len(cc):
            CCGroups.append([(ids, cc)])
        elif 300 > len(cc) >= 100:
            TimeUsedM.append((ids, cc))
        else:
            TimeUsedF.append((ids, cc))
    random.shuffle(TimeUsedM)
    [CCGroups.append(TimeUsedM[i:i + 20]) for i in range(0, len(TimeUsedM), 20)]
    random.shuffle(TimeUsedF)
    [CCGroups.append(TimeUsedF[i:i + 200]) for i in range(0, len(TimeUsedF), 200)]
    return CCGroups, 'c'


def GetAttribution(network, getType):
    geneIDRaw = []
    genesNum = 0
    genomeIDs = set()
    for recodeName in network.vs['name']:
        genome = recodeName.split('|')[0]
        genomeIDs.add(genome)
        genesNum += 1
        if getType:
            geneIDRaw.append(recodeName)
    genomesNum = len(genomeIDs)
    if geneIDRaw:
        geneIDs = sorted(geneIDRaw)
        attribute = ','.join([' '.join(geneIDs), str(genesNum), ' '.join(genomeIDs), str(genomesNum)])
    else:
        attribute = ','.join([str(genesNum), ' '.join(genomeIDs), str(genomesNum)])
    return attribute, genomesNum


def GetConnectedOfCCsV2(queues, graph, queueList, detectedM, modularWeight, GammaCoefficients):
    HHGC = HHN(graph)
    ConnectedList = []
    ccc = queueList
    while ccc:
        ccc = HHGC.RunCommunityDetection(ConnectedList, ccc, detectedM, modularWeight, GammaCoefficients)
    ConnectedTuple = tuple(ConnectedList)
    queues.put((ConnectedTuple, []))


def GetConnectedOfCCs(queues, graph, connectedComponent, detectedM, modularWeight, GammaCoefficient):
    HHGC = HHN(graph)
    ConnectedList = []
    secondList = HHGC.RunCommunityDetection(ConnectedList, connectedComponent, detectedM, modularWeight,
                                            GammaCoefficient)
    ConnectedTuple = tuple(ConnectedList)
    queues.put((ConnectedTuple, secondList))


def ConstructedHnCC(ssnGraph, detectedM, modularWeight, parallelNum, gc):
    """
    :param: ccc: connected components communities
    :return:hnl: hierarchical network communities list
    """
    hnn = HHN(ssnGraph)
    ConnectedComponents = hnn.splitConnectedComponents()
    turn = 1
    NodeListALl = []
    EdgeListAll = []
    genesIDsListAll = []
    genomeIDsListAll = []
    genesNumListAll = []
    genomeNumListAll = []
    modularityListAll = []
    CommunityAll, continued = SplitTasks1(ConnectedComponents)
    while CommunityAll:
        MyQueues = mp.Manager().Queue()
        if len(CommunityAll) <= int(parallelNum):
            parallelNum = len(CommunityAll)
        else:
            parallelNum = parallelNum
        pro = mp.Pool(processes=int(parallelNum))
        for cc in CommunityAll:
            if continued is None:
                pro.apply_async(
                    func=GetConnectedOfCCs, args=(MyQueues, ssnGraph, cc, detectedM, modularWeight, gc))
            else:
                pro.apply_async(
                    func=GetConnectedOfCCsV2, args=(MyQueues, ssnGraph, cc, detectedM, modularWeight, gc))
        connects, SecondList = GetConnectedList(MyQueues, len(CommunityAll), turn)
        if SecondList:
            MaxCC = sorted(SecondList, key=lambda X: len(X[1]), reverse=True)[0][1]
            print('Max Nodes Numbers of CC %d ' % len(MaxCC))
            if len(MaxCC) < 5000:
                CommunityAll, continued = SplitTasks2(SecondList)
            else:
                CommunityAll, continued = SplitTasks1(SecondList)
        else:
            CommunityAll = SecondList
        pro.close()
        pro.join()
        NodeList, EdgeList, genesIDsList, genomeIDsList, genesNumList, genomeNumList, modularityList = connects
        NodeListALl.extend(NodeList)
        EdgeListAll.extend(EdgeList)
        genesIDsListAll.extend(genesIDsList)
        genomeIDsListAll.extend(genomeIDsList)
        genesNumListAll.extend(genesNumList)
        genomeNumListAll.extend(genomeNumList)
        modularityListAll.extend(modularityList)
        turn += 1
    hhm = igraph.Graph()
    attributesDic = dict(
        geneIDs=genesIDsListAll, genomeIDs=genomeIDsListAll, genesNum=genesNumListAll, genomesNum=genomeNumListAll,
        modularity=modularityListAll)
    hhm.add_vertices(NodeListALl, attributes=attributesDic)
    hhm.add_edges(EdgeListAll)
    return hhm


def ExtractOGSorted(hhn, ssn, seqInf):
    inputGenomeNum = len(seqInf['GenomeToUsed'])
    genomeNum = inputGenomeNum
    ogsAll = {}
    adListsAll = []
    hhnO = HHN(ssn, hhn)
    Since = time.time()
    Ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(Since))
    progressbar_widgets_set = [
        'Extract OGs: %s ' % Ts, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(),
        ' | ',
        progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=inputGenomeNum + 1)
    bar.start()
    num = 0
    while genomeNum >= 0:
        ogsInf, hhn, adLists = hhnO.ExtractOG(genomeNum)
        adListsAll.extend(adLists)
        hhnO = HHN(ssn, hhn)
        genomeNum -= 1
        num += 1
        bar.update(num)
        for k, v in ogsInf.items():
            ogsAll[k] = v
    bar.finish()
    done = time.time()
    print('Extract OGs %s' % TimeUsedComputation(Since, done))
    for AD in adListsAll:
        if AD['name'] in ogsAll.keys():
            del ogsAll[AD['name']]
    return hhn, ogsAll


def GetGenesID(nodes, genesIds, deletedNode):
    nodeSelected = []
    for node in nodes:
        if node['Event'] is not None and node.degree() > 1:
            for neighbor in node.neighbors():
                if neighbor['genesNum'] < node['genesNum']:
                    nodeSelected.append(neighbor)
                    deletedNode.append(neighbor)
        else:
            for recode in node['geneIDs'].split(' '):
                genesIds.append(recode)
    return nodeSelected


def GetGenesIDs(selectedNode):
    neighbors = [selectedNode]
    genesIDs = []
    childNode = []
    while neighbors:
        neighbors = GetGenesID(neighbors, genesIDs, childNode)
    return genesIDs, childNode


def GetSeqs(SeqIDs, SeqInformation):
    fastaFormat = ''
    for seqID in SeqIDs:
        seq = SeqInformation['SequencesRecode'][seqID]
        fastaFormat += '>%s\n%s\n' % (seqID, seq)
    return fastaFormat


def WriteOGFiles(seqInf, hhn, ogInf, wd):
    num = 0
    ogStatic = ''
    genesNum = 0
    genesSet = set()
    handle_og = []
    since = time.time()
    recode_node_dic = {}
    Orthogroups_Sequences_dir = os.path.join(wd, 'Orthogroups_Sequences')
    try:
        os.mkdir(Orthogroups_Sequences_dir)
    except FileExistsError:
        pass
    for ogName, inf in ogInf.items():
        recodeGenes = inf[2].split(' ')
        genes = ' '.join([seqInf['GenesRecode'][recode] for recode in recodeGenes])
        output_og_name = '{}{:0>6d}'.format('OG', num + 1)
        recode_node_dic[ogName] = output_og_name
        ogStatic += '%s\t%d\t%d\t%s\n' % (output_og_name, inf[0], inf[1], genes)
        og_seq = ''.join(['>%s\n%s\n' % (gene, seqInf['SequencesRecode'][gene]) for gene in genes.split(' ')])
        OutputFile(os.path.join(Orthogroups_Sequences_dir, '%s.fasta' % output_og_name), og_seq)
        genesNum += inf[1]
        num += 1
        [genesSet.add(gene) for gene in recodeGenes]
        if inf[-1] is not None and inf[1] > 1:
            graphs = hhn.subgraph(inf[-1])
            handle_og.append(graphs)
    done = time.time()
    OutputFile(os.path.join(wd, 'OGFile_coalescence_SameGenome.txt'), ogStatic)
    print('Total Genes Set Numbers: %d' % len(genesSet))
    print('Total Genes Numbers: %d' % genesNum)
    print('Write OGs %s' % TimeUsedComputation(since, done))
    return recode_node_dic, handle_og


class EstimateOrthologs:
    @staticmethod
    def WriteOGPairwise(seqInf, ogPairs, wd):
        since = time.time()
        ogPairsS = ''
        for ogPair in ogPairs:
            ogPairsS += '%s\t%s\t%s\t%s\n' % (
                seqInf['GenesRecode'][ogPair[0]].split('|')[1], seqInf['GenesRecode'][ogPair[1]].split('|')[1],
                ogPair[0], ogPair[1])
        done = time.time()
        print('Write OG Pairwise %s ' % TimeUsedComputation(since, done))
        OGPairName = os.path.join(wd, 'OGPairwise_coalescence_SameGenome.txt')
        OutputFile(OGPairName, ogPairsS)

    @staticmethod
    def splitTaskForOG(OGSubList):
        OGGroups = []
        TemList = []
        for Index, OGSub in enumerate(OGSubList):
            selectNodes = OGSub.vs.select(_degree=1)
            nodeCombs = combinations(selectNodes, 2)
            nodeCombs_to_list = list(nodeCombs)
            TemList.extend([(OGSub, nodes) for nodes in nodeCombs_to_list])
            if len(TemList) >= 50000 or Index == (len(OGSubList) - 1):
                OGGroups.append(TemList.copy())
                TemList.clear()
        return OGGroups


def SelectedOGFromOneNode(OGList):
    OGPairList = []
    for OG in OGList:
        for OG2 in OGList:
            if OG.split('|')[0] != OG2.split('|')[0]:
                OGPairList.append((OG, OG2))
    return OGPairList


def GetOGPairwise(queue, selectSubList):
    ogPairs = []
    for selectSub in selectSubList:
        selectNodes = selectSub.vs.select(_degree=1)
        nodeComb = combinations(selectNodes, 2)
        for q, r in list(nodeComb):
            shortest_path = q.get_shortest_paths(r)
            LCAIndex = sorted(shortest_path[0][1:-1], key=lambda X: int(selectSub.vs[X]['genesNum']), reverse=True)[0]
            LCAE = selectSub.vs[LCAIndex]['Event']
            queryGenes = q['geneIDs'].strip().split(' ')
            referGenes = r['geneIDs'].strip().split(' ')
            ogPairs.extend(SelectedOGFromOneNode(queryGenes))
            ogPairs.extend(SelectedOGFromOneNode(referGenes))
            if LCAE == 'I':
                for query in queryGenes:
                    for refer in referGenes:
                        if query.split('|')[0] != refer.split('|')[0]:
                            ogPairs.append((query, refer))
            else:
                pass
    queue.put(ogPairs)


def GetPairwise(Queue, TaskNum):
    progressbar_widgets_set = [
        'Get Pairwise: ', progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(),
        ' | ',
        progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=TaskNum)
    MyPairsList = []
    doneNum = 0
    bar.start()
    while True:
        MyPairs = Queue.get()
        for pairwise in MyPairs:
            MyPairsList.append(pairwise)
        doneNum += 1
        bar.update(doneNum)
        if doneNum >= TaskNum:
            break
    bar.finish()
    return MyPairsList


def GetOrthologsFromOGsPP(SeqInformation, OGGraphs, WorkingD, Nthreads):
    since = time.time()
    OGPairsListAll = splitTaskForOG(OGGraphs)
    RunThreads = int(Nthreads) if int(Nthreads) <= len(OGPairsListAll) else len(OGPairsListAll)
    ogPairs_queue = mp.Manager().Queue()
    pool2 = mp.Pool(processes=RunThreads)
    for OGPairsList in OGPairsListAll:
        pool2.apply_async(func=GetOGPairwise, args=(ogPairs_queue, OGPairsList))
    ogPairs = GetPairwise(ogPairs_queue, len(OGPairsListAll))
    pool2.close()
    pool2.join()
    done = time.time()
    EstimateOrthologs.WriteOGPairwise(SeqInformation, ogPairs, WorkingD)
    print('OGPairs %s' % TimeUsedComputation(since, done))


def splitTaskForOG(OGSubList):
    TimeUsedS = []
    TimeUsedM = []
    TimeUsedF = []
    TimeUsedU = []
    OGGroups = []
    for OGSub in OGSubList:
        if 500 <= OGSub.vcount():
            TimeUsedS.append(OGSub)
        elif 200 <= OGSub.vcount() < 500:
            TimeUsedM.append(OGSub)
        elif 100 <= OGSub.vcount() < 200:
            TimeUsedF.append(OGSub)
        else:
            TimeUsedU.append(OGSub)
    if TimeUsedS:
        TimeUsedS = sorted(TimeUsedS, key=lambda X: X.vcount(), reverse=True)
        [OGGroups.append(TimeUsedS[i:i + 1]) for i in range(0, len(TimeUsedS), 1)]
    if TimeUsedM:
        random.shuffle(TimeUsedM)
        [OGGroups.append(TimeUsedM[i:i + 10]) for i in range(0, len(TimeUsedM), 10)]
    if TimeUsedF:
        random.shuffle(TimeUsedF)
        [OGGroups.append(TimeUsedF[i:i + 100]) for i in range(0, len(TimeUsedF), 100)]
    if TimeUsedU:
        random.shuffle(TimeUsedU)
        [OGGroups.append(TimeUsedU[i:i + 1000]) for i in range(0, len(TimeUsedU), 1000)]
    return OGGroups


def output_final_graph(seqInfDic, final_graph, output_graph_name, og_node_recode_dic):
    for index in final_graph.vs.indices:
        final_graph_node = final_graph.vs[index]
        node_name = final_graph_node['name']
        if node_name in og_node_recode_dic.keys():
            final_graph_node['name'] = node_recode_dic[node_name]
        else:
            final_graph_node['name'] = 'None'
        recode_ids = ' '.join([seqInfDic['GenesRecode'][ids] for ids in final_graph_node['geneIDs'].split(' ')])
        recode_genome_ids = ' '.join(
            [seqInfDic['GenomeRecodeInf'][genome_ids] for genome_ids in final_graph_node['genomeIDs'].split(' ')])
        final_graph_node['geneIDs'] = recode_ids
        final_graph_node['genomeIDs'] = recode_genome_ids
    final_graph.write_gml(output_graph_name)


if __name__ == '__main__':
    parameters = get_parameters()
    InputDir = parameters['inputDir']
    continuedMethod = parameters['continued']
    E_values = parameters['e_values']
    threads = parameters['searchThreads']
    seqSearchingMethod = parameters['searchMethod']
    OutputGraph = parameters['outputFile']
    weight = parameters['weightTypes']
    if weight not in ['NBS', 'Sim', 'BS']:
        print('Error! Attribute does not exist :%s\nCurrent Version only provide Sim ,BS , NBS' % weight)
        sys.exit(0)
    NetworkThreads = parameters['GraphThreads']
    detectionMethod = parameters['communityDetectMethod']
    MostDistance = parameters['distance']
    Overlaps = parameters['speciesOverlap']
    WorkingDir = os.path.join(InputDir, 'WorkingDirectory')
    gc = parameters['gammaCoefficient']
    parameter = ' '.join(sys.argv)
    print('Current parameter %s' % parameter)
    print('Whole process PID%d' % os.getpid())
    try:
        os.mkdir(WorkingDir)
    except FileExistsError:
        pass
    # ----------------Step1--------------------------------------
    # prepare genome data
    processSinceTime = time.time()
    if continuedMethod:
        if os.path.isfile(os.path.join(WorkingDir, 'ssn.gml')):
            ssnR = igraph.read(os.path.join(WorkingDir, 'ssn.gml'))
            SeqInf = MatricesLoad('SeqInf', 'Information', WorkingDir)
        else:
            if os.path.exists(os.path.join(WorkingDir, 'BlastResults')) and os.path.isfile(
                    os.path.join(WorkingDir, 'SeqInfInformation.pic')):
                BlastResultPath = os.path.join(WorkingDir, 'BlastResults')
                SeqInf = MatricesLoad('SeqInf', 'Information', WorkingDir)
                SSNSince = time.time()
                SSNStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(SSNSince))
                print('[%s]SSN Construction' % SSNStartTime)
                MatrixDir = os.path.join(WorkingDir, 'matrices')
                try:
                    os.mkdir(MatrixDir)
                except FileExistsError:
                    pass
                GetBHMatrixParallel(BlastResultPath, MostDistance, SeqInf, MatrixDir, NetworkThreads)
                GetConnectedMatrixParallel(SeqInf, MostDistance, MatrixDir, NetworkThreads)
                ssnR = BuildSSN_parallel(MatrixDir, SeqInf, NetworkThreads)
                SSNDone = time.time()
                SSNEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(SSNDone))
                SSNTimeStr = TimeUsedComputation(SSNSince, SSNDone)
                print('[%s]SSN constructed %s' % (SSNEndTime, SSNTimeStr))
                ssnNetwork = os.path.join(WorkingDir, 'ssn.gml')
                ssnR.write_gml(ssnNetwork)
                ssnR = igraph.read(ssnNetwork)
            else:
                # ----------------Step2--------------------------------------
                # Running sequences homology searching
                print('Sequences Similarity Network Not Found!\n It Will Reconstruct SSN')
                prepareSince = time.time()
                prepareStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(prepareSince))
                print('[%s]Prepare genome data' % prepareStartTime)
                GenomeFiles = GetInputGenomeInf(InputDir)
                SeqInf, BlastInput = RecodeGenomeSeq(GenomeFiles, InputDir, WorkingDir)
                MatricesDumpy(SeqInf, 'SeqInf', 'Information', WorkingDir)
                prepareDone = time.time()
                prepareEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(prepareDone))
                prepareTimeStr = TimeUsedComputation(prepareSince, prepareDone)
                print('[%s]Done\n%s' % (prepareEndTime, prepareTimeStr))
                searchSince = time.time()
                searchStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(searchSince))
                print('[%s]Running sequences homology searching' % searchStartTime)
                BlastResultPath = SequencesSearchBlast(BlastInput, seqSearchingMethod, WorkingDir, E_values, threads)
                searchDone = time.time()
                searchEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(searchDone))
                searchTimeStr = TimeUsedComputation(searchSince, searchDone)
                print('[%s]Searching Done\n%s' % (searchEndTime, searchTimeStr))
                # ----------------Step3--------------------------------------
                # Construction of original sequences similarity network based on the sequences searching results
                SSNSince = time.time()
                SSNStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(SSNSince))
                print('[%s]SSN Construction' % SSNStartTime)
                MatrixDir = os.path.join(WorkingDir, 'matrices')
                try:
                    os.mkdir(MatrixDir)
                except FileExistsError:
                    pass
                GetBHMatrixParallel(BlastResultPath, MostDistance, SeqInf, MatrixDir, NetworkThreads)
                GetConnectedMatrixParallel(SeqInf, MostDistance, MatrixDir, NetworkThreads)
                ssnR = BuildSSN_parallel(MatrixDir, SeqInf, NetworkThreads)
                SSNDone = time.time()
                SSNEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(SSNDone))
                SSNTimeStr = TimeUsedComputation(SSNSince, SSNDone)
                print('[%s]SSN constructed %s' % (SSNEndTime, SSNTimeStr))
                ssnNetwork = os.path.join(WorkingDir, 'ssn.gml')
                ssnR.write_gml(ssnNetwork)
                ssnR = igraph.read(ssnNetwork)
    else:
        processSinceTime = time.time()
        prepareSince = time.time()
        prepareStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(prepareSince))
        print('[%s]Prepare genome data' % prepareStartTime)
        GenomeFiles = GetInputGenomeInf(InputDir)
        SeqInf, BlastInput = RecodeGenomeSeq(GenomeFiles, InputDir, WorkingDir)
        MatricesDumpy(SeqInf, 'SeqInf', 'Information', WorkingDir)
        prepareDone = time.time()
        prepareEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(prepareDone))
        prepareTimeStr = TimeUsedComputation(prepareSince, prepareDone)
        print('[%s]Done\n%s' % (prepareEndTime, prepareTimeStr))
        # ----------------Step2--------------------------------------
        # Running sequences homology searching
        searchSince = time.time()
        searchStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(searchSince))
        print('[%s]Running sequences homology searching' % searchStartTime)
        BlastResultPath = SequencesSearchBlast(BlastInput, seqSearchingMethod, WorkingDir, E_values, threads)
        searchDone = time.time()
        searchEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(searchDone))
        searchTimeStr = TimeUsedComputation(searchSince, searchDone)
        print('[%s]Searching Done\n%s' % (searchEndTime, searchTimeStr))
        # ----------------Step3--------------------------------------
        # Construction of original sequences similarity network based on the sequences searching results
        SSNSince = time.time()
        SSNStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(SSNSince))
        print('[%s]SSN Construction' % SSNStartTime)
        MatrixDir = os.path.join(WorkingDir, 'matrices')
        try:
            os.mkdir(MatrixDir)
        except FileExistsError:
            pass
        MatrixSince = time.time()
        MatrixStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(MatrixSince))
        print('[%s]Get RBH Matrices' % MatrixStartTime)
        GetBHMatrixParallel(BlastResultPath, MostDistance, SeqInf, MatrixDir, NetworkThreads)
        GetConnectedMatrixParallel(SeqInf, MostDistance, MatrixDir, NetworkThreads)
        MatrixDone = time.time()
        MatrixEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(MatrixDone))
        MatrixTimeStr = TimeUsedComputation(MatrixSince, MatrixDone)
        print('[%s]Matrices Construction Done!' % MatrixEndTime)
        ssnR = BuildSSN_parallel(MatrixDir, SeqInf, NetworkThreads)
        SSNDone = time.time()
        SSNEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(SSNDone))
        SSNTimeStr = TimeUsedComputation(SSNSince, SSNDone)
        print('[%s]SSN constructed %s' % (SSNEndTime, SSNTimeStr))
        ssnNetwork = os.path.join(WorkingDir, 'ssn.gml')
        ssnR.write_gml(ssnNetwork)
        ssnR = igraph.read(ssnNetwork)
    # ----------------Step4--------------------------------------
    # community detection of ssn based on the leidenalg algorithm for construction of hierarchical network
    if continuedMethod:
        if os.path.isfile(os.path.join(WorkingDir, 'hm.gml')):
            hg = igraph.read(os.path.join(WorkingDir, 'hm.gml'))
        else:
            CDSince = time.time()
            CDStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(CDSince))
            print('[%s]Community Detection' % CDStartTime)
            hg = ConstructedHnCC(ssnR, detectionMethod, weight, NetworkThreads, gc)
            hg.write_gml(os.path.join(WorkingDir, 'hm.gml'))
            CDDone = time.time()
            CDEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(CDDone))
            CDTimeStr = TimeUsedComputation(CDSince, CDDone)
            print('[%s]Community Detection Done %s' % (CDEndTime, CDTimeStr))
    else:
        CDSince = time.time()
        CDStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(CDSince))
        print('[%s]Community Detection' % CDStartTime)
        hg = ConstructedHnCC(ssnR, detectionMethod, weight, NetworkThreads, gc)
        hg.write_gml(os.path.join(WorkingDir, 'hm.gml'))
        CDDone = time.time()
        CDEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(CDDone))
        CDTimeStr = TimeUsedComputation(CDSince, CDDone)
        print('[%s]Community Detection Done\n%s' % (CDEndTime, CDTimeStr))
    # ----------------Step5--------------------------------------
    # mark the evolution events on th homologs hierarchical network
    MarkSince = time.time()
    MarkStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(MarkSince))
    print('[%s]Mark the evolution events on HHNC and output final file' % MarkStartTime)
    hm = HHN(ssnR, hg)
    hhng = hm.GetEvolutionEvents(Overlaps)
    hhng.write_gml(os.path.join(WorkingDir, 'hmm.gml'))
    hhngO = hhng.copy()
    hhngF, OGDic = ExtractOGSorted(hhng, ssnR, SeqInf)
    MarkDone = time.time()
    MarkEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(MarkDone))
    MarkTimeStr = TimeUsedComputation(MarkSince, MarkDone)
    node_recode_dic, OGForPairwise = WriteOGFiles(SeqInf, hhngO, OGDic, WorkingDir)
    output_final_graph(SeqInf, hhngF, OutputGraph, node_recode_dic)
    print('[%s]Done %s' % (MarkEndTime, MarkTimeStr))
    PairsTime = time.time()
    PairsStartTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(PairsTime))
    print('[%s]Get Pairwise for all ogs' % PairsStartTime)
    GetOrthologsFromOGsPP(SeqInf, OGForPairwise, WorkingDir, NetworkThreads)
    PairsDone = time.time()
    PairsEndTime = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(PairsDone))
    PairsTimeStr = TimeUsedComputation(MarkSince, PairsDone)
    print('[%s]Get Pairwise Done %s' % (PairsEndTime, PairsTimeStr))
    processEndTime = time.time()
    processUseTime = TimeUsedComputation(processSinceTime, processEndTime)
    print('Whole process %s' % processUseTime)
