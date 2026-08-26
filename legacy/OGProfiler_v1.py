import functools
import multiprocessing as mp
import os
import pickle as pic
import random
import subprocess
import sys
import time
from itertools import combinations, groupby
from optparse import OptionParser, OptionGroup
import ete3
import igraph
import leidenalg as la
import numpy as np
import progressbar
import scipy.sparse as sparse
from scipy.optimize import curve_fit
from scipy.special import comb
from operator import methodcaller


def get_parameters():
    usage = "OGProfiler.py -i <input dir> <Options>"
    opt = OptionParser(usage=usage)
    group0 = OptionGroup(opt, "General options")
    group0.add_option(
        '-i', '--in', type=str, dest='input_dir', default=False,
        help='specify the directory including all genome files')
    group0.add_option(
        '-o', '--out', type=str, dest='output_dir', default=os.getcwd(),
        help=f"specify a output directory default: {os.getcwd()}")
    group0.add_option(
        '-x', '--extension', type=str, dest='extension', default='fna',
        help=f"Extension of genome files; default: fna")
    group0.add_option(
        '-s', '--search_method', dest='search_method', choices=['diamond', 'mmseqs', 'blastp'], default='diamond',
        help='Homologs searching methods: blastp, mmseqs, diamond. Both mmseqs and diamond are sensitive mode ('
             'default: diamond)')
    group0.add_option(
        '-e', '--evalue', type=float, dest='e_values', default=0.001,
        help="cut-off of E value")
    group0.add_option(
        '-t', '--threads', type=int, dest='search_threads', default=8,
        help='Homologs searching threads. default: 8')
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
        help='leidenalg community detection contains 4 algorithm based on th modularity. mvp: Modularity Vertex '
             'Partition \nrbcv: RB Configuration Vertex Partition\nrber: RBERVertexPartition\ncmp: CPM Vertex '
             'Partition\n')
    group1.add_option(
        '-a', '--network_threads', type=int, dest='network_threads', default=8,
        help='parallel running threads of network analysis')
    group1.add_option(
        '-g', '--gamma_coefficient', type=float, dest='gamma_coefficient', default=1.0,
        help='coefficient of primary resolution parameters of community detection')
    group1.add_option(
        '--so', '--species_overlap', type=int, dest='species_overlap', default=0,
        help='Species overlap threads for inference of speciation events in homologs hierarchical networks')
    group1.add_option(
        '-r', '--refined', action='store_true', dest='refined', default=False,
        help='refined network structure based on the trees')
    opt.add_option_group(group0)
    opt.add_option_group(group1)
    options, args = opt.parse_args()
    input_dir = options.input_dir
    output_dir = options.output_dir
    extension = options.extension
    continued = options.continued
    search_method = options.search_method
    e_values = options.e_values
    search_threads = options.search_threads
    distance = options.distance_defined
    weight_types = options.weight_type
    community_detect_method = options.community_detect_method
    graph_threads = options.network_threads
    gamma_coefficient = options.gamma_coefficient
    species_overlap = options.species_overlap
    refine = options.refined
    parameters_dict = dict(
        input_dir=input_dir, output_dir=output_dir, extension=extension, search_method=search_method, e_values=e_values,
        search_threads=search_threads, distance=distance, weight_types=weight_types,
        community_detect_method=community_detect_method, graph_threads=graph_threads, continued=continued,
        gamma_coefficient=gamma_coefficient, species_overlap=species_overlap, refined=refine)
    return parameters_dict


# prepare genomes sequences file methods
# ------------------------------------------------------------------------------------------------------------------
def output_file(fileName, content):
    with open(fileName, 'w') as f:
        f.writelines(content)


def parse_fasta(fasta_name):
    seqDict = {}
    with open(fasta_name) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = header.__next__()[1:].strip().split(' ')[0]
            seq = "".join(s.strip() for s in faiter.__next__())
            if header in seqDict:
                sys.exit('FASTA contains multiple entries with the same name')
            else:
                seqDict[header] = seq
    return seqDict


def time_used(info=''):
    def timer(function):
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            start = time.perf_counter() if sys.version[0] == '3' else time.clock()
            results = function(*args, **kwargs)
            end = time.perf_counter() if sys.version[0] == '3' else time.clock()
            time_use = end - start
            print(
                f'[{info}]: {time_use // 3600:.0f}h {(time_use % 3600) // 60:.0f}m {((time_use % 3600) % 60) % 60:.4f}s')
            return results

        return wrapper

    return timer


def read_file(file_name):
    with open(file_name) as file:
        data = file.readlines()
    return data


def matrices_dumpy(matrix, matrixType, name, path):
    with open(os.path.join(path, '%s%s.pic' % (matrixType, name)), 'wb') as picFile:
        pic.dump(matrix, picFile, protocol=3)


def matrices_load(matrixType, name, path):
    with open(os.path.join(path, '%s%s.pic' % (matrixType, name)), 'rb') as picFile:
        m = pic.load(picFile)
    return m


def get_input_genome_inf(genome_path, suffix):
    input_genomes = [
        genome_file for genome_file in os.listdir(genome_path) if
        genome_file.split('.')[-1] == suffix]
    if len(input_genomes) < 4:
        print("Error! the number of genomes were less than 4!\n Please check your datasets or suffix of genome file")
        sys.exit(0)
    else:
        return input_genomes


# homologous searching methods
# ------------------------------------------------------------------------------------------------------------------
def make_blast_db(referGenomes, referGenomesPath, methods, WorkingDirectory):
    db_path = os.path.join(WorkingDirectory, 'localDB')
    try:
        os.mkdir(db_path)
    except FileExistsError:
        pass
    db_list = []
    for Refer in referGenomes:
        refer = os.path.join(referGenomesPath, Refer)
        db = os.path.join(db_path, 'db%s' % Refer)
        if methods == 'blastp':
            make_db_cmd = ' '.join(['makeblastdb', '-dbtype', 'prot', '-in', refer, '-out', db, '>', '/dev/null'])
        elif methods == 'diamond':
            make_db_cmd = ' '.join(['diamond', 'makedb', '--in', refer, '--db', db, '--threads', '1', '>', '/dev/null'])
        else:
            print('Error! Not support searching methods')
            sys.exit(0)
        db_list.append(db)
        process = subprocess.Popen(make_db_cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stdin=subprocess.PIPE)
        process.wait()
    return db_list


def diamond(queue, queryGenome, queryGenomePath, BlastResultDir, DB, e_value):
    blast_query = os.path.join(queryGenomePath, queryGenome)
    blast_file_path = os.path.join(
        BlastResultDir, '%s-%s.out' % (queryGenome.split('.')[0], DB.split('db')[-1].split('.')[0]))
    command = ' '.join([
        'diamond', 'blastp', '--more-sensitive', '-p', '1', '-q', blast_query, '-d', '%s.dmnd' % DB,
        '--evalue', str(e_value), '-f', '6', '--out', blast_file_path, '--quiet'])
    pro = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pro.wait()
    queue.put((blast_file_path, pro.returncode))


def mmseqs(queue, queryGenome, queryGenomePath, BlastResultDir, DB, e_value):
    blast_query = os.path.join(queryGenomePath, queryGenome)
    blast_file_path = os.path.join(
        BlastResultDir, '%s-%s.out' % (queryGenome.split('.')[0], DB.split('/')[-1].split('.')[0]))
    temp_path = os.path.join(BlastResultDir, '%s-%s' % (queryGenome.split('.')[0], DB.split('/')[-1].split('.')[0]))
    command = ' '.join([
        'mmseqs', 'easy-search', blast_query, DB, blast_file_path, temp_path, '--threads', '1', '-v', '1',
        '--format-mode', '0', '--remove-tmp-files', '-s', '7.5', '-e', str(e_value)])
    pro = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pro.wait()
    queue.put((blast_file_path, pro.returncode))


def blastp(queue, queryGenome, queryGenomePath, BlastResultDir, DB, e_value):
    blast_query = os.path.join(queryGenomePath, queryGenome)
    blast_file_path = os.path.join(
        BlastResultDir, '%s-%s.out' % (queryGenome.split('.')[0], DB.split('db')[-1].split('.')[0]))
    command = ' '.join(
        ['blastp', '-outfmt', '6', '-query', blast_query, '-db', DB, '-evalue', str(e_value), '-out', blast_file_path])
    pro = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pro.wait()
    queue.put((blast_file_path, pro.returncode))


def run_blast_search_parallel(queryGenomes, queryPath, searchMethod, DBList, BlastResultsPath, E_Values, thread):
    queue = mp.Manager().Queue()
    pros = mp.Pool(processes=int(thread))
    for query in queryGenomes:
        for DB in DBList:
            if searchMethod == 'blastp':
                pros.apply_async(func=blastp, args=(queue, query, queryPath, BlastResultsPath, DB, E_Values))
            elif searchMethod == 'diamond':
                pros.apply_async(func=diamond, args=(queue, query, queryPath, BlastResultsPath, DB, E_Values))
            elif searchMethod == 'mmseqs':
                pros.apply_async(func=mmseqs, args=(queue, query, queryPath, BlastResultsPath, DB, E_Values))
            else:
                pass
    pros.close()
    print_parallel_bar(queue, 'Homologs Searching:', len(queryGenomes) * len(DBList))
    pros.join()


# Get original matrices between genomes pair
# ------------------------------------------------------------------------------------------------------------------
def read_blast_results(BlastFileName, SeqInformation):
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
    similarity_pair = {}
    bit_scores_matrices = {}
    query_genome = BlastFileName.split('/')[-1].split('.')[0].split('-')[0]
    refer_genome = BlastFileName.split('/')[-1].split('.')[0].split('-')[1]
    length_sparse_matrix = sparse.lil_matrix(
        (SeqInformation['SpeciesGeneNum'][query_genome], SeqInformation['SpeciesGeneNum'][refer_genome]))
    bit_score_sparse_matrix = length_sparse_matrix.copy()
    blast_hit_lines = read_file(BlastFileName)
    for hitLine in blast_hit_lines:
        line_items = hitLine.strip('\n').split('\t')
        try:
            query_id = line_items[0]
            refer_id = line_items[1]
            query_gene_index = int(query_id.split('|')[-1].split('g')[-1])
            refer_gene_index = int(refer_id.split('|')[-1].split('g')[-1])
            q_match = int(line_items[7]) - int(line_items[6]) + 1
            r_match = int(line_items[9]) - int(line_items[8]) + 1
            query_length = SeqInformation['SeqLengthInf'][query_id]
            refer_length = SeqInformation['SeqLengthInf'][refer_id]
            length_product = query_length * refer_length
            query_cover = (q_match / query_length) * 100
            refer_cover = (r_match / refer_length) * 100
            similarity = float(line_items[2])
            bit_cores = float(line_items[-1])
            coverage = min([query_cover, refer_cover])
            if query_id != refer_id:
                length_sparse_matrix[query_gene_index, refer_gene_index] = length_product
                bit_score_sparse_matrix[query_gene_index, refer_gene_index] = bit_cores
                similarity_pair['%s-%s' % (query_id, refer_id)] = similarity
                bit_scores_matrices['%s-%s' % (query_id, refer_id)] = bit_cores
        except ValueError:
            print(hitLine)
    return similarity_pair, bit_scores_matrices, length_sparse_matrix, bit_score_sparse_matrix


class BSN:
    @staticmethod
    def retain_top_data(length_matrix, bit_score_matrix):
        """
        :param length_matrix: Length product of query and hit
        :param bit_score_matrix: BitScores of query and hit
        :return TopLength: sequences length in top 5% bit scores bins from different length bins:
        :return TopBitScores BitScores in top 5% bit scores bins from different length bins:
        """
        length_array = [
            length_matrix[row, col] for row, col in
            zip(length_matrix.nonzero()[0], length_matrix.nonzero()[1])]
        bit_scores_array = [
            bit_score_matrix[row, col]
            for row, col in zip(bit_score_matrix.nonzero()[0], bit_score_matrix.nonzero()[1])]
        length_sorted = sorted([(length, bit)
                                for length, bit in zip(length_array, bit_scores_array)], key=lambda X: X[0])
        hit_numbers = len(length_sorted)
        if hit_numbers < 100:
            length_bin_raw = [element[0] for element in length_sorted]
            bit_bin_raw = [element[1] for element in length_sorted]
            return length_bin_raw, bit_bin_raw
        length_top = []
        bit_top = []
        scale = 1000 if hit_numbers > 5000 else (200 if hit_numbers > 1000 else 20)
        for bin_index in range(0, len(length_sorted), scale):
            bit_bin = [binLen[-1] for binLen in length_sorted[bin_index:bin_index + scale + 1]]
            cutoff = np.percentile(np.array(bit_bin), 95)
            for bin_element in length_sorted[bin_index:bin_index + scale + 1]:
                if bin_element[-1] >= cutoff:
                    length_top.append(bin_element[0])
                    bit_top.append(bin_element[1])
        top5_length = np.array(length_top)
        top5_bit_scores = np.array(bit_top)
        return top5_length, top5_bit_scores

    @staticmethod
    def fitness_function(x, a, b):
        y = a * np.log10(x) + b
        return y

    @staticmethod
    def get_fitness_para(xData, yData):
        a, b = curve_fit(BSN.fitness_function, xData, np.log10(yData))[0]
        return a, b

    @staticmethod
    def normalized_function(length_product, bit_score, a, b):
        normalized_bit_scores = bit_score / ((10 ** b) * (length_product.power(a)))
        return normalized_bit_scores


# Get connected matrices
class MatrixHandle:
    def __init__(self, matrix_dir, sequences_inf, genomeI):
        self.matrix_dir = matrix_dir
        self.sequences_inf = sequences_inf
        self.genomeI = genomeI

    def get_specie_nb_matrix(self, blast_hit_dir):
        """
        normalized bit-scores
        :param blast_hit_dir: Blast results files directory
        """
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            genome_pair = f'{self.genomeI}-{genomeJ}'
            blast_file_name = os.path.join(blast_hit_dir, f'{genome_pair}.out')
            sim_dic, bit_dic, length_matrix, bit_matrix = read_blast_results(blast_file_name, self.sequences_inf)
            matrices_dumpy(sim_dic, 'similarity', genome_pair, self.matrix_dir)
            matrices_dumpy(bit_dic, 'bitScores', genome_pair, self.matrix_dir)
            if bit_matrix.count_nonzero() > 1:
                length_top_array, bit_top_array = BSN.retain_top_data(length_matrix, bit_matrix)
                a, b = BSN.get_fitness_para(length_top_array, bit_top_array)
                nb_matrix = BSN.normalized_function(length_matrix, bit_matrix, a, b)
                nb_matrix[np.isnan(nb_matrix)] = 0
                matrices_dumpy(sparse.lil_matrix(nb_matrix), 'NB', genome_pair, self.matrix_dir)
                matrices_dumpy(length_matrix.tolil(), 'length', genome_pair, self.matrix_dir)
            else:
                print('Warning! There are no hits in all-to-all blast except the same genes:%s-%s' % (
                    self.sequences_inf['GenomeRecodeInf'][self.genomeI],
                    self.sequences_inf['GenomeRecodeInf'][genomeJ]))
                nb_matrix = sparse.lil_matrix(
                    (self.sequences_inf['SpeciesGeneNum'][self.genomeI],
                     self.sequences_inf['SpeciesGeneNum'][genomeJ]))
                length_matrix = nb_matrix.copy()
                matrices_dumpy(sparse.lil_matrix(nb_matrix), 'NB', genome_pair, self.matrix_dir)
                matrices_dumpy(length_matrix.tolil(), 'length', genome_pair, self.matrix_dir)

    def get_best_hit(self, accepted=1e-3):
        """
        best_hit_i: Best Hit index for each genes from Genome I in original normalized matrix
        Matrices: normalized matrix for each species pairs
        BitHitBestIndexList: best hit for each genome pair
        """
        best_hit_i = -1 * np.ones(self.sequences_inf['SpeciesGeneNum'][self.genomeI])
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            if self.genomeI == genomeJ:
                continue
            row_matrix = matrices_load('NB', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            best_hit_index_list = []
            best_hit_row_list = []
            for row_num in range(self.sequences_inf['SpeciesGeneNum'][self.genomeI]):
                bit_sores_object = row_matrix.getrowview(row_num)
                if bit_sores_object.nnz > 0:
                    max_values = max(bit_sores_object.data[0])
                    best_hit_i[row_num] = max_values if max_values > best_hit_i[row_num] else best_hit_i[row_num]
                    best_hits_col_index = [
                        colIndex for colIndex, bitValue in zip(bit_sores_object.rows[0], bit_sores_object.data[0])
                        if bitValue > max_values - accepted]
                    best_hit_index_list.extend(best_hits_col_index)
                    best_hit_row_list.extend(np.full(len(best_hits_col_index), row_num, dtype=int))
            best_hit_matrix_for_dump = sparse.csr_matrix(
                (np.ones(len(best_hit_row_list)), (best_hit_row_list, best_hit_index_list)),
                shape=row_matrix.get_shape())
            matrices_dumpy(best_hit_matrix_for_dump, 'BH', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
        # paralogs
        row_matrix = matrices_load('NB', '%s-%s' % (self.genomeI, self.genomeI), self.matrix_dir)
        best_hit_index_list = []
        best_hit_row_list = []
        for rowNum in range(self.sequences_inf['SpeciesGeneNum'][self.genomeI]):
            bit_sores_object = row_matrix.getrowview(rowNum)
            if bit_sores_object.nnz > 0:
                best_hits_col_index = [
                    colIndex for colIndex, bitValue in zip(bit_sores_object.rows[0], bit_sores_object.data[0])
                    if bitValue > best_hit_i[rowNum] - accepted]
                best_hit_index_list.extend(best_hits_col_index)
                best_hit_row_list.extend(np.full(len(best_hits_col_index), rowNum, dtype=int))
        best_hit_matrix_for_dump = sparse.csr_matrix(
            (np.ones(len(best_hit_row_list)), (best_hit_row_list, best_hit_index_list)),
            shape=row_matrix.get_shape())
        matrices_dumpy(best_hit_matrix_for_dump, 'BH', '%s-%s' % (self.genomeI, self.genomeI), self.matrix_dir)

    def rbh_matrix(self):
        """
        calculate reciprocal best hit for each genome pairs
        bh: best hit matrix
        rbh: reciprocal best hit matrices for  each genome pairs
        """
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            if genomeJ == self.genomeI:
                rbh = sparse.csr_matrix(
                    (self.sequences_inf['SpeciesGeneNum'][self.genomeI], self.sequences_inf['SpeciesGeneNum'][genomeJ]))
            else:
                bh_1 = matrices_load('BH', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
                bh_2 = matrices_load('BH', '%s-%s' % (genomeJ, self.genomeI), self.matrix_dir)
                rbh = bh_1.multiply(bh_2.transpose())
            matrices_dumpy(rbh, 'RBHs', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)

    def sorted_connected(self):
        """
        RBHM: reciprocal best hits matrix
        NBMatrices: Normalized Bit-scores matrices
        SeqInformation: sequences information
        ConnectedMatrixPair : connected genes pair for each genome pairs
        """
        genome_i_rbh_min = 1e9 * np.ones(self.sequences_inf['SpeciesGeneNum'][self.genomeI])
        best_hit = np.zeros(self.sequences_inf['SpeciesGeneNum'][self.genomeI])
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            if genomeJ == self.genomeI:
                continue
            nb_metric_ij = matrices_load('NB', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            nb_metric_ij_csr = nb_metric_ij.tocsr()
            rbh_metric_ij = matrices_load('RBHs', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            for row_num in range(self.sequences_inf['SpeciesGeneNum'][self.genomeI]):
                if nb_metric_ij.getrowview(row_num).nnz > 0:
                    best_hit[row_num] = max(best_hit[row_num], max(nb_metric_ij.getrowview(row_num).data[0]))
                if rbh_metric_ij[row_num].nnz > 0:
                    genome_i_rbh_min[row_num] = min(
                        min(nb_metric_ij_csr[row_num, rbh_metric_ij[row_num].indices].data), genome_i_rbh_min[row_num])
        indices = genome_i_rbh_min > 1e8
        genome_i_rbh_min[indices] = best_hit[indices] + 1e-6
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            nb_metric_ij = matrices_load('NB', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            nb_metric_ij_csr = nb_metric_ij.tocsr()
            retained_hit_ij = []
            retained_hit_row = []
            retained_hit_col = []
            for row_num in range(self.sequences_inf['SpeciesGeneNum'][self.genomeI]):
                if nb_metric_ij_csr[row_num].nnz > 0:
                    for values, index in zip(nb_metric_ij_csr[row_num].data, nb_metric_ij_csr[row_num].indices):
                        if values >= genome_i_rbh_min[row_num]:
                            retained_hit_col.append(index)
                            retained_hit_row.append(row_num)
                            retained_hit_ij.append(values)
            connected_matrix = sparse.csr_matrix(
                (retained_hit_ij, (retained_hit_row, retained_hit_col)), shape=nb_metric_ij_csr.get_shape())
            matrices_dumpy(connected_matrix, 'CM', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)

    def sorted_connected_rbh(self):
        """
        Only Selected reciprocal best hits for SSN construction
        retained_hit_ij: normalized bit-scores values of retained hit
        retained_hit_row: row index of retained hit
        retained_hit_col: col index of retained hit
        """
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            rbh_metric_ij = matrices_load('RBHs', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            nb_metric_ij = matrices_load('NB', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            rbh_ij_csr = rbh_metric_ij.tocsr()
            retained_hit_ij = []
            retained_hit_row = []
            retained_hit_col = []
            for row_num in range(self.sequences_inf['SpeciesGeneNum'][self.genomeI]):
                if self.genomeI == genomeJ:
                    for index in rbh_ij_csr[row_num].indices:
                        retained_hit_col.append(index)
                        retained_hit_row.append(row_num)
                        values = nb_metric_ij[row_num, index]
                        retained_hit_ij.append(values)
                else:
                    if rbh_ij_csr[row_num].nnz > 0:
                        for index in rbh_ij_csr[row_num].indices:
                            retained_hit_col.append(index)
                            retained_hit_row.append(row_num)
                            values = nb_metric_ij[row_num, index]
                            retained_hit_ij.append(values)
                connected_matrix = sparse.csr_matrix(
                    (retained_hit_ij, (retained_hit_row, retained_hit_col)), shape=nb_metric_ij.get_shape())
                matrices_dumpy(connected_matrix, 'CM', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)

    def sorted_connected_rh(self):
        """Only Selected reciprocal hits for SSN construction
        nb_metric_ij: Normalized bit-scores matrix of pair between species I and species J
        retained_hit_ij: normalized bit-scores values of retained hit
        retained_hit_row: row index of retained hit
        retained_hit_col: col index of retained hit
        """
        for genomeJ in self.sequences_inf['GenomeToUsed']:
            nb_metric_ij = matrices_load('NB', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            nb_metric_ij_csr = nb_metric_ij.tocsr()
            nb_metric_ji = matrices_load('NB', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)
            nb_metric_ji_csr = nb_metric_ji.tocsr()
            retained_hit_ij = []
            retained_hit_row = []
            retained_hit_col = []
            for row_num in range(self.sequences_inf['SpeciesGeneNum'][self.genomeI]):
                if self.genomeI == genomeJ:
                    for index in nb_metric_ij_csr[row_num].indices:
                        retained_hit_col.append(index)
                        retained_hit_row.append(row_num)
                        values = nb_metric_ij[row_num, index]
                        retained_hit_ij.append(values)
                else:
                    for values, index in zip(nb_metric_ij_csr[row_num].data, nb_metric_ij_csr[row_num].indices):
                        if nb_metric_ji_csr[row_num, index] > 0:
                            retained_hit_col.append(index)
                            retained_hit_row.append(row_num)
                            retained_hit_ij.append(values)
                connected_matrix = sparse.csr_matrix(
                    (retained_hit_ij, (retained_hit_row, retained_hit_col)), shape=nb_metric_ij.get_shape())
                matrices_dumpy(connected_matrix, 'CM', '%s-%s' % (self.genomeI, genomeJ), self.matrix_dir)


def get_bh_matrix_all(queue, SpeciesI, MostD, searchResultPath, recodeInf, matrixDirectory):
    matrices_object = MatrixHandle(matrixDirectory, recodeInf, SpeciesI)
    matrices_object.get_specie_nb_matrix(searchResultPath)
    if MostD == 'lrb' or MostD == 'arb':
        matrices_object.get_best_hit()
    else:
        pass
    queue.put((SpeciesI, 0))


def get_bh_matrix_parallel(ResultPath, mostDistances, reInf, matDirectory, matrixThreads):
    my_queue = mp.Manager().Queue()
    p = mp.Pool(processes=int(matrixThreads))
    for iSpecies in reInf['GenomeToUsed']:
        p.apply_async(func=get_bh_matrix_all, args=(my_queue, iSpecies, mostDistances, ResultPath, reInf, matDirectory))
    print_parallel_bar(my_queue, 'Get Matrices', len(reInf['GenomeToUsed']))
    p.close()
    p.join()


def print_parallel_bar(queue, parallelType, taskNum):
    t = time.time()
    ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t))
    progressbar_widgets_set = [
        '[%s]%s: ' % (ts, parallelType), progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=taskNum)
    bar.start()
    done_num = 0
    while True:
        cmd, complete_num = queue.get()
        if complete_num == 0:
            pass
        else:
            print('%s Error!' % cmd)
        done_num += 1
        bar.update(done_num)
        if done_num >= taskNum:
            break
    bar.finish()


def get_connected_matrix_species_i(queue, SpeciesI, MostD, recodeInf, matrixDirectory):
    matrices_about = MatrixHandle(matrixDirectory, recodeInf, SpeciesI)
    if MostD == 'lrb':
        matrices_about.rbh_matrix()
        matrices_about.sorted_connected()
    elif MostD == 'arb':
        matrices_about.sorted_connected_rbh()
    elif MostD == 'ar':
        matrices_about.sorted_connected_rh()
    else:
        print('Error: %s is not defined method for the most distances' % MostD)
        pass
    queue.put((SpeciesI, 0))


def get_connected_matrix_parallel(inf, mostDistance, matDir, nt):
    queue = mp.Manager().Queue()
    pro = mp.Pool(processes=int(nt))
    for iGenome in inf['GenomeToUsed']:
        pro.apply_async(func=get_connected_matrix_species_i, args=(queue, iGenome, mostDistance, inf, matDir))
    pro.close()
    print_parallel_bar(queue, 'Construction of Connected Matrices:', len(inf['GenomeToUsed']))
    pro.join()


def get_connections(queues, SpeciesI, SpeciesJ, MatrixDirectory, SpeciesInf):
    raw_graph = igraph.Graph()
    connect_matrix = matrices_load('CM', '%s-%s' % (SpeciesI, SpeciesJ), MatrixDirectory)
    sp = matrices_load('similarity', '%s-%s' % (SpeciesI, SpeciesJ), MatrixDirectory)
    bp = matrices_load('bitScores', '%s-%s' % (SpeciesI, SpeciesJ), MatrixDirectory)
    node_name = []
    edge_list = []
    sim_list = []
    bs_list = []
    nbs_list = []
    NBSRList = []
    SimRList = []
    BSRList = []
    if SpeciesI != SpeciesJ:
        connectMatrixR = matrices_load('CM', '%s-%s' % (SpeciesJ, SpeciesI), MatrixDirectory).transpose()
        SPR = matrices_load('similarity', '%s-%s' % (SpeciesJ, SpeciesI), MatrixDirectory)
        BPR = matrices_load('bitScores', '%s-%s' % (SpeciesJ, SpeciesI), MatrixDirectory)
        connectMatrixAdd = connect_matrix + connectMatrixR
        for rowNum in range(SpeciesInf['SpeciesGeneNum'][SpeciesI]):
            sourceNode = '%s|g%d' % (SpeciesI, rowNum)
            if sourceNode not in node_name:
                node_name.append(sourceNode)
            if connectMatrixAdd[rowNum].nnz > 0:
                for Index in connectMatrixAdd[rowNum].indices:
                    NBS = connect_matrix[rowNum, Index]
                    NBS1 = connectMatrixR[rowNum, Index]
                    if NBS > 0 and NBS1 > 0:
                        NBSValues = np.float64(NBS).item()
                        NBSR = np.float64(NBS1).item()
                        Similarity = sp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        SimilarityR = SPR['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                        BitScores = bp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        BitScoresR = BPR['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                    elif NBS > 0 and NBS1 == 0:
                        NBSValues = np.float64(NBS).item()
                        NBSR = 0
                        Similarity = sp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                        SimilarityR = 0
                        BitScores = bp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
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
                    if targetNode not in node_name:
                        node_name.append(targetNode)
                    edge_list.append([sourceNode, targetNode])
                    sim_list.append(Similarity)
                    SimRList.append(SimilarityR)
                    bs_list.append(BitScores)
                    BSRList.append(BitScoresR)
                    nbs_list.append(NBSValues)
                    NBSRList.append(NBSR)
    else:
        for rowNum in range(SpeciesInf['SpeciesGeneNum'][SpeciesI]):
            sourceNode = '%s|g%d' % (SpeciesI, rowNum)
            if sourceNode not in node_name:
                node_name.append(sourceNode)
            if connect_matrix[rowNum].nnz > 0:
                for Index in connect_matrix[rowNum].indices:
                    targetNode = '%s|g%d' % (SpeciesJ, Index)
                    if [targetNode, sourceNode] not in edge_list:
                        if connect_matrix[Index, rowNum] > 0:
                            NBS = connect_matrix[rowNum, Index]
                            NBS1 = connect_matrix[Index, rowNum]
                            Similarity = sp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            SimilarityR = sp['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                            BitScores = bp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            BitScoresR = bp['%s|g%d-%s|g%d' % (SpeciesJ, Index, SpeciesI, rowNum)]
                        else:
                            NBS = connect_matrix[rowNum, Index]
                            NBS1 = 0
                            Similarity = sp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            SimilarityR = 0
                            BitScores = bp['%s|g%d-%s|g%d' % (SpeciesI, rowNum, SpeciesJ, Index)]
                            BitScoresR = 0
                        if targetNode not in node_name:
                            node_name.append(targetNode)
                        edge_list.append([sourceNode, targetNode])
                        sim_list.append(Similarity)
                        SimRList.append(SimilarityR)
                        bs_list.append(BitScores)
                        BSRList.append(BitScoresR)
                        nbs_list.append(NBS)
                        NBSRList.append(NBS1)
    raw_graph.add_vertices(node_name)
    attributesDic = dict(Sim=sim_list, BS=bs_list, NBS=nbs_list, NBSRev=NBSRList, BSRev=BSRList, SimRev=SimRList)
    raw_graph.add_edges(edge_list, attributes=attributesDic)
    queues.put(raw_graph)


def build_ssn_parallel(MatrixDirectory, speciesInf, buildThreads, workingDir):
    my_queue = mp.Manager().Queue()
    pro = mp.Pool(processes=int(buildThreads))
    for i in range(len(speciesInf['GenomeToUsed'])):
        for j in range(i, len(speciesInf['GenomeToUsed'])):
            i_species, j_species = speciesInf['GenomeToUsed'][i], speciesInf['GenomeToUsed'][j]
            pro.apply_async(
                func=get_connections, args=(my_queue, i_species, j_species, MatrixDirectory, speciesInf))
    pro.close()
    pair_num = comb(len(speciesInf['GenomeToUsed']), 2) + len(speciesInf['GenomeToUsed'])
    ssn = UnionGraphs(my_queue, pair_num)
    pro.join()
    ssn.write_gml(os.path.join(workingDir, 'ssn.gml'))
    return ssn


def UnionGraphs(graphLists, taskNum):
    """
    Union graphs of genome pairwise
    :param graphLists: graph queue
    :param taskNum: task numbers
    :return: sequences similarity network
    """
    DoneTask = 0
    T = time.time()
    Ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(T))
    progressbar_widgets_set = [
        '[%s]Construction of SSN: ' % Ts, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer()]
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


def MappingGammaForCC(coefficient_g, lengthCC):
    coefficient_g = (10 ** (int(np.log10(lengthCC) - 4))) * coefficient_g
    index1 = float(np.log10(lengthCC))
    gamma = float(coefficient_g) * (10 ** (-index1 + 1))
    return gamma


# HHN build Methods
# ------------------------------------------------------------------------------------------------------------------
class HHN(igraph.Graph):
    def __init__(self, ssn=None, hhn=None, **attr):
        super(HHN, self).__init__(**attr)
        self.ssn = ssn
        self.hhn = hhn

    def splitConnectedComponents(self):
        connected_components = [('%d' % num, cc) for num, cc in enumerate(self.ssn.components()) if len(cc) >= 2]
        return connected_components

    def SpecificBoard(self, Nodes, detectionMethods, edgeAttribute, gammaC, startNum, stopNum):
        community_graph = self.ssn.subgraph(Nodes)
        left = 0
        right = MappingGammaForCC(gammaC, len(Nodes))
        method_optional = dict(
            mvp=la.ModularityVertexPartition, rbcv=la.RBConfigurationVertexPartition,
            rber=la.RBERVertexPartition, cmp=la.CPMVertexPartition)
        n = 1
        p = None
        while True:
            partition = la.find_partition(
                community_graph, partition_type=method_optional[detectionMethods],
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

    def BipartiteGraphs(self, communityNode, detectionMethods, edgeAttribute, GammaCo, startNum, stopNum):
        community_graph = self.ssn.subgraph(communityNode)
        method_optional = dict(
            mvp=la.ModularityVertexPartition, rbcv=la.RBConfigurationVertexPartition,
            rber=la.RBERVertexPartition, cmp=la.CPMVertexPartition)
        if len(communityNode) >= 1000:
            left, right, dn, p, quality = self.SpecificBoard(communityNode, detectionMethods, edgeAttribute, GammaCo,
                                                             startNum,
                                                             stopNum)
        else:
            left, right, dn, p = 0, 1, 1, None
        if p is None:
            n = 1
            resolution_parameter = (left + right) / 2
            partition = la.find_partition(
                community_graph, partition_type=method_optional[detectionMethods],
                n_iterations=10, weights=edgeAttribute, resolution_parameter=resolution_parameter)
            while n <= 1000:
                partition = la.find_partition(
                    community_graph, partition_type=method_optional[detectionMethods],
                    n_iterations=10, weights=edgeAttribute, resolution_parameter=resolution_parameter)
                if len(partition.subgraphs()) < startNum:
                    left = resolution_parameter
                    resolution_parameter = (left + right) / 2
                elif len(partition.subgraphs()) > stopNum:
                    right = resolution_parameter
                    resolution_parameter = (left + right) / 2
                else:
                    break
                n += 1
            partition_g = partition.subgraphs()
            quality = partition.q
        else:
            partition_g = p
            quality = 'None'
        return partition_g, quality

    def RunCommunityDetection(self, connected, modularise, detectedMethod, weightType, GC):
        """
        implement community detection
        :param connected: connected list
        :param detectedMethod: communities detection methods
        :param weightType: edge weight type for community detection
        :param modularise: connected components contains community on corresponding iteration
        :param GC: gamma coefficient
        :return: communities from current iteration
        """
        communities_list = []
        for ID, modular in modularise:
            sn = ID
            sn_attribution, modular_genome_num = GetAttribution(self.ssn.subgraph(modular), True)
            if modular_genome_num > 1 and len(modular) >= 10000:
                stop_num = 50
                start_num = 20
                communities, quality = self.BipartiteGraphs(modular, detectedMethod, weightType, GC, start_num,
                                                            stop_num)
                if len(communities) >= 2:
                    tns = ''
                    for num, community in enumerate(communities):
                        c_cs_communities = [int(node) for node in community.vs['index']]
                        tn = '%s-%d' % (ID, num)
                        tns += tn + ' '
                        communities_list.append((tn, c_cs_communities))
                    connected_inf = ','.join([sn, tns, sn_attribution, str(quality)])
                else:
                    connected_inf = ','.join([sn, '+', sn_attribution, '0'])
            elif modular_genome_num > 1 and 2 < len(modular) < 10000:
                stop_num = 2
                start_num = 2
                communities, quality = self.BipartiteGraphs(modular, detectedMethod, weightType, GC, start_num,
                                                            stop_num)
                if len(communities) == 2:
                    tns = ''
                    for num, community in enumerate(communities):
                        c_cs_communities = [int(node) for node in community.vs['index']]
                        tn = '%s-%d' % (ID, num)
                        tns += tn + ' '
                        communities_list.append((tn, c_cs_communities))
                    connected_inf = ','.join([sn, tns, sn_attribution, str(quality)])
                else:
                    connected_inf = ','.join([sn, '+', sn_attribution, '0'])
            elif modular_genome_num > 1 and len(modular) == 2:
                tns = '%s-0 %s-1' % (ID, ID)
                connected_inf = ','.join([sn, tns, sn_attribution, '0'])
                communities_list.append(('%s-0' % ID, [modular[0]]))
                communities_list.append(('%s-1' % ID, [modular[1]]))
            else:
                connected_inf = ','.join([sn, '+', sn_attribution, '0'])
            connected.append(connected_inf)
        return communities_list

    def GetEvolutionEvents(self, overlapThreads, workingDir):
        """
        Mark evolution events for nodes on network
        :param workingDir: working directory
        :param overlapThreads: numbers of species overlap
        :return: network
        """
        selected_vertex = self.hhn.vs.select(_degree_gt=1)
        for vertex in selected_vertex:
            sgs = set(vertex['genomeIDs'].split(' '))
            vertex_genes_num = vertex['genesNum']
            gcs = []
            for node in vertex.neighbors():
                if int(node['genesNum']) < int(vertex_genes_num):
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
        selected_vertex0 = self.hhn.vs.select(_degree=0)
        for vertex in selected_vertex0:
            source_node_genome = vertex['genomesNum']
            source_node_genes = vertex['genesNum']
            if source_node_genome == source_node_genes:
                vertex['Event'] = 'I'
            else:
                pass
        self.hhn.write_gml(os.path.join(workingDir, 'hmm.gml'))
        return self.hhn

    def ExtractOG(self, GenomeNum, event, max_num):
        if GenomeNum == 1:
            ogNodes = self.hhn.vs.select(Event='None')
        elif GenomeNum == 0:
            ogNodes = self.ssn.vs.select(_degree=0)
        else:
            if max_num > 0:
                ogNodes = self.hhn.vs.select(Event=event, genomesNum=GenomeNum, genesNum_lt=max_num, genesNum_gt=4)
            else:
                ogNodes = self.hhn.vs.select(Event=event, genomesNum=GenomeNum)
        nodesDeleted = []
        OGInf = {}
        if GenomeNum > 0:
            ogNodes = sorted(ogNodes, key=lambda X: X['genesNum'], reverse=True)
        for ogNode in ogNodes:
            if ogNode['name'] not in [n['name'] for n in nodesDeleted]:
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
                    # if OneNodeCoalescence(ogNode):
                    #     parent, adjacent = OneNodeCoalescence(ogNode)[0]
                    #     genesIDList.extend(adjacent['geneIDs'].strip().split(' '))
                    #     deletedVertex.extend([ogNode, adjacent])
                    #     nodeCC.extend([parent['name'], adjacent['name']])
                    #     ogNode = parent
                    # else:
                    #     ogNode = ogNode
                    nodesDeleted.extend(deletedVertex)
                genesIDs = ' '.join([i for i in set(genesIDList)])
                OGInf[ogNode['name']] = (GenomeNum, len(genesIDList), genesIDs, nodeCC)
        if GenomeNum > 1:
            self.hhn.delete_vertices(nodesDeleted)
        return OGInf, self.hhn


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
    T = time.time()
    Ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(T))
    progressbar_widgets_set = [
        '[%s]Split connected components %d: ' % (Ts, turns), progressbar.Percentage(), progressbar.Bar('#'),
        progressbar.Timer()]
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


def get_connected_of_ccs(queues, graph, connectedComponent, detectedM, modularWeight, GammaCoefficient):
    HHGC = HHN(graph)
    ConnectedList = []
    secondList = HHGC.RunCommunityDetection(ConnectedList, connectedComponent, detectedM, modularWeight,
                                            GammaCoefficient)
    connected_tuple = tuple(ConnectedList)
    queues.put((connected_tuple, secondList))


def ExtractOGSorted(hhn, ssn, seqInf, evolution_event, stop_at, max_numbers):
    inputGenomeNum = len(seqInf['GenomeToUsed'])
    genomeNum = inputGenomeNum
    ogsAll = {}
    hhnO = HHN(ssn, hhn)
    Ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
    progressbar_widgets_set = [
        '[%s]Extract OGs: ' % Ts, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=inputGenomeNum + 1)
    bar.start()
    num = 0
    while genomeNum >= stop_at:
        ogsInf, hhn = hhnO.ExtractOG(genomeNum, evolution_event, max_numbers)
        hhnO = HHN(ssn, hhn)
        genomeNum = genomeNum - 1
        num += 1
        bar.update(num)
        for k, v in ogsInf.items():
            ogsAll[k] = v
    bar.finish()
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


# refined nodes event by building trees
# --------------------------------------------------------------------------------------------------------------------
def find_fixing_nodes(ssn, hmm, sequences_inf, Orthogroups_Sequences_dir):
    max_genome = len(sequences_inf['GenomeToUsed']) * 10
    hhn, nodes_dic = ExtractOGSorted(hmm, ssn, sequences_inf, 'III-1', 2, max_genome)
    og_name = []
    num = 0
    for og, og_inf in nodes_dic.items():
        genes = og_inf[2].split(' ')
        og_seq = ''.join(
            ['>%s\n%s\n' % (gene, sequences_inf['SequencesRecode'][sequences_inf['GenesRecode'][gene]]) for gene in
             genes])
        output_og_name = '{}{:0>6d}'.format('Nodes', num + 1)
        hhn.vs.select(name=og)[0]['name'] = output_og_name
        og_file_name = os.path.join(Orthogroups_Sequences_dir, '%s.fasta' % output_og_name)
        output_file(og_file_name, og_seq)
        og_name.append(output_og_name)
        num += 1
    return hhn, og_name


def alignments(queue, og_file, alignments_path):
    og_aln = os.path.join(alignments_path, '%s.aln.fasta' % os.path.basename(og_file).split('.')[0])
    # og_trimal = os.path.join(alignments_path, '%s.trimal.fasta' % os.path.basename(og_file).split('.')[0])
    aln_cmd = ' '.join(['mafft', '--anysymbol', og_file, '>', og_aln])
    aln_pro = subprocess.Popen(aln_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    aln_pro.wait()
    # trimal_cmd = ' '.join(['trimal', '-in', og_aln, '-out', og_trimal, '-automated1']) trimal_pro =
    # subprocess.Popen(trimal_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    # trimal_pro.wait()
    queue.put((aln_cmd, aln_pro.returncode))


def fast_tree(queue, alignments_files, genes_trees_out_dir):
    prefix = os.path.join(genes_trees_out_dir, os.path.basename(alignments_files).split('.')[0])
    treefile = '%s.nwk' % prefix
    build_tree_command = ' '.join(['/home/licongjian/fastTree/FastTree', '-lg', alignments_files, '>', treefile])
    pro = subprocess.Popen(build_tree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pro.wait()
    queue.put((build_tree_command, pro.returncode))


def alignmentsPP(nodes_fasta_list, seq_path, align_path, align_threads):
    aln_build_pool = mp.Pool(int(align_threads))
    myQueue = mp.Manager().Queue()
    for orthologs in nodes_fasta_list:
        orthologs_file = os.path.join(seq_path, '%s.fasta' % orthologs)
        aln_build_pool.apply_async(alignments, args=(myQueue, orthologs_file, align_path))
    print_parallel_bar(myQueue, 'Alignments :', len(nodes_fasta_list))
    aln_build_pool.close()
    aln_build_pool.join()


def ConstructionTrees(aln_list, aln_path, trees_path, tree_parallel_threads):
    trees_build_pool = mp.Pool(int(tree_parallel_threads))
    myQueue = mp.Manager().Queue()
    for aln in aln_list:
        aln_file = os.path.join(aln_path, '%s.aln.fasta' % aln)
        trees_build_pool.apply_async(fast_tree, args=(myQueue, aln_file, trees_path))
    print_parallel_bar(myQueue, 'Building Trees', len(aln_list))
    trees_build_pool.close()
    trees_build_pool.join()


def DecipherTrees(queue, newick_tree_file):
    tree = ete3.Tree(newick_tree_file)
    Root = tree.get_midpoint_outgroup()
    tree.set_outgroup(Root)
    attributions = []
    multi_branch = []
    for node in tree.traverse("preorder"):
        genesIDs = node.get_leaf_names() if len(node.get_leaf_names()) > 0 else [node.name]
        if node.is_root():
            nodeName = os.path.basename(newick_tree_file).split('.')[0]
        else:
            nodeName = ' '.join(sorted(genesIDs))
        genomes_parent_set = set([gene.split('|')[0] for gene in genesIDs])
        genomes = ' '.join(sorted(list(genomes_parent_set)))
        genomeNum = len(genomes_parent_set)
        genesNum = len(genesIDs)
        child = node.children
        if 3 > len(child) > 0:
            child1 = child[0]
            child2 = child[1]
            leaf_1 = child1.get_leaf_names() if len(child1.get_leaf_names()) > 0 else [child1.name]
            leaf_2 = child2.get_leaf_names() if len(child2.get_leaf_names()) > 0 else [child2.name]
            leaf_1_genome_set = set([gene.split('|')[0] for gene in leaf_1])
            leaf_2_genome_set = set([gene.split('|')[0] for gene in leaf_2])
            if leaf_1_genome_set & leaf_2_genome_set == genomes_parent_set:
                event = 'II'
            elif leaf_1_genome_set & leaf_2_genome_set == set():
                event = 'I'
            else:
                if leaf_1_genome_set == genomes_parent_set or leaf_2_genome_set == genomes_parent_set:
                    event = 'III-1'
                else:
                    event = 'III-2'
            edge1 = '-'.join([nodeName, ' '.join(sorted(leaf_1))])
            edge2 = '-'.join([nodeName, ' '.join(sorted(leaf_2))])
            edges = ';'.join([edge1, edge2])
        elif len(child) >= 3:
            [multi_branch.append(n.name) for n in child]
            edges = 'None'
            event = 'None'
        else:
            edges = 'None'
            event = 'None'
        if node.is_leaf() and node.name in multi_branch:
            pass
        else:
            genes = ' '.join(genesIDs)
            attr = ','.join([nodeName, edges, genes, genomes, str(genesNum), str(genomeNum), event])
            attributions.append(attr)
    queue.put(attributions)


def GetTreesStructurePP(NodesRefined, g, treePath, cpus):
    trees = mp.Pool(int(cpus))
    myQueue = mp.Manager().Queue()
    for node in NodesRefined:
        treeName = os.path.join(treePath, '%s.nwk' % node)
        trees.apply_async(func=DecipherTrees, args=(myQueue, treeName,))
        # DecipherTrees(myQueue, treeName, max_num_dic)
    BuildGraph(myQueue, g, len(NodesRefined))
    trees.close()
    trees.join()


def BuildGraph(queue, original_graph, taskNum):
    NodeListALl = []
    EdgeListAll = []
    genesIDsListAll = []
    genomeIDsListAll = []
    genesNumListAll = []
    genomeNumListAll = []
    Events = []
    T = time.time()
    Ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(T))
    progressbar_widgets_set = [
        '[%s]Deciphering trees: ' % Ts, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer()]
    bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=taskNum)
    bar.start()
    DoneTask = 0
    while True:
        attribute = queue.get()
        for attr in attribute:
            nodeName, edges, genesID, genomeID, genesNum, genomeNum, events = attr.split(',')
            if edges == 'None':
                pass
            else:
                edge1 = [edges.split(';')[0].split('-')[0], edges.split(';')[0].split('-')[1]]
                edge2 = [edges.split(';')[1].split('-')[0], edges.split(';')[1].split('-')[1]]
                EdgeListAll.append(edge1)
                EdgeListAll.append(edge2)
            if 'Node' in nodeName:
                original_graph.vs.select(name=nodeName)[0]['Event'] = events
            else:
                NodeListALl.append(nodeName)
                genesIDsListAll.append(genesID)
                genomeIDsListAll.append(genomeID)
                genesNumListAll.append(int(genesNum))
                genomeNumListAll.append(int(genomeNum))
                Events.append(events)
        DoneTask += 1
        bar.update(DoneTask)
        if DoneTask >= taskNum:
            break
    bar.finish()
    attributesDic = dict(
        geneIDs=genesIDsListAll, genomeIDs=genomeIDsListAll, genesNum=genesNumListAll, genomesNum=genomeNumListAll,
        Event=Events)
    original_graph.add_vertices(NodeListALl, attributes=attributesDic)
    original_graph.add_edges(EdgeListAll)


def check_point(nodes_list, path, file_type):
    nodes = []
    for node in nodes_list:
        tree = os.path.join(path, '%s.%s' % (node, file_type))
        if os.path.isfile(tree):
            if os.path.getsize(tree) > 0:
                pass
            else:
                nodes.append(node)
        else:
            nodes.append(node)
    return nodes


def RefinedNodesEvents(ssn, hmm, sequences_inf, refined_dir, refined_threads):
    seq_dir = os.path.join(refined_dir, 'seq')
    aln_dir = os.path.join(refined_dir, 'aln')
    tree_dir = os.path.join(refined_dir, 'tree')
    for dir_ in [refined_dir, seq_dir, aln_dir, tree_dir]:
        os.makedirs(dir_, exist_ok=True)
    hhn, og_name = find_fixing_nodes(ssn, hmm, sequences_inf, seq_dir)
    hhn.write_gml(os.path.join(refined_dir, 'hnn.gml'))
    re_aln_og = check_point(og_name, aln_dir, 'aln.fasta')
    re_tree_og = check_point(og_name, tree_dir, 'nwk')
    if re_aln_og:
        alignmentsPP(re_aln_og, seq_dir, aln_dir, refined_threads)
    if re_tree_og:
        ConstructionTrees(re_tree_og, aln_dir, tree_dir, refined_threads)
    GetTreesStructurePP(og_name, hhn, tree_dir, refined_threads)
    try:
        hhn.vs['id'] = list(map(str, hhn.vs['id']))
    except KeyError:
        pass
    hhn.write_gml(os.path.join(refined_dir, 'hhn.gml'))
    return hhn


# output orthologous groups and hhm graph file
# ----------------------------------------------------------------------------------------------------------------
def WriteOGFiles(seqInf, hhn, ogInf, wd):
    num = 0
    ogStatic = ''
    genesNum = 0
    genesSet = set()
    handle_og = []
    recode_node_dic = {}
    Orthogroups_Sequences_dir = os.path.join(wd, 'Orthogroups_Sequences')
    os.makedirs(Orthogroups_Sequences_dir, exist_ok=True)
    for ogName, inf in ogInf.items():
        recodeGenes = inf[2].split(' ')
        genes = ' '.join([seqInf['GenesRecode'][recode] for recode in recodeGenes])
        output_og_name = '{}{:0>6d}'.format('OG', num + 1)
        recode_node_dic[ogName] = output_og_name
        ogStatic += '%s\t%d\t%d\t%s\n' % (output_og_name, inf[0], inf[1], genes)
        og_seq = ''.join(['>%s\n%s\n' % (gene, seqInf['SequencesRecode'][gene]) for gene in genes.split(' ')])
        output_file(os.path.join(Orthogroups_Sequences_dir, '%s.fasta' % output_og_name), og_seq)
        genesNum += inf[1]
        num += 1
        [genesSet.add(gene) for gene in recodeGenes]
        if inf[-1] is not None and inf[1] > 1:
            graphs = hhn.subgraph(inf[-1])
            handle_og.append(graphs)
    output_file(os.path.join(wd, 'OGFile_coalescence_SameGenome.txt'), ogStatic)
    print('Total Genes Set Numbers: %d' % len(genesSet))
    print('Total Genes Numbers: %d' % genesNum)
    return recode_node_dic, handle_og


def output_final_graph(seqInfDic, final_graph, output_graph_name, og_node_recode_dic):
    for index in final_graph.vs.indices:
        final_graph_node = final_graph.vs[index]
        node_name = final_graph_node['name']
        if node_name in og_node_recode_dic.keys():
            final_graph_node['name'] = og_node_recode_dic[node_name]
        else:
            final_graph_node['name'] = 'None'
        recode_ids = ' '.join([seqInfDic['GenesRecode'][ids] for ids in final_graph_node['geneIDs'].split(' ')])
        recode_genome_ids = ' '.join(
            [seqInfDic['GenomeRecodeInf'][genome_ids] for genome_ids in final_graph_node['genomeIDs'].split(' ')])
        final_graph_node['geneIDs'] = recode_ids
        final_graph_node['genomeIDs'] = recode_genome_ids
    final_graph.write_gml(output_graph_name)


# estimated for pairwise genes
# ----------------------------------------------------------------------------------------------------------------
class EstimateOrthologs:
    @staticmethod
    def WriteOGPairwise(seqInf, ogPairs, wd):
        ogPairsS = ''
        for ogPair in ogPairs:
            ogPairsS += '%s\t%s\t%s\t%s\n' % (
                seqInf['GenesRecode'][ogPair[0]].split('|')[1], seqInf['GenesRecode'][ogPair[1]].split('|')[1],
                ogPair[0], ogPair[1])
        OGPairName = os.path.join(wd, 'OGPairwise_coalescence_SameGenome.txt')
        output_file(OGPairName, ogPairsS)

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
    T = time.time()
    Ts = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(T))
    progressbar_widgets_set = [
        '[%s]Get Pairwise: ' % Ts, progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer()]
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


def GetOrthologsFromOGsPP(SeqInformation, OGGraphs, WorkingD, NumThread):
    OGPairsListAll = splitTaskForOG(OGGraphs)
    RunThreads = int(NumThread) if int(NumThread) <= len(OGPairsListAll) else len(OGPairsListAll)
    ogPairs_queue = mp.Manager().Queue()
    pool2 = mp.Pool(processes=RunThreads)
    for OGPairsList in OGPairsListAll:
        pool2.apply_async(func=GetOGPairwise, args=(ogPairs_queue, OGPairsList))
    ogPairs = GetPairwise(ogPairs_queue, len(OGPairsListAll))
    pool2.close()
    pool2.join()
    EstimateOrthologs.WriteOGPairwise(SeqInformation, ogPairs, WorkingD)


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


class Call:
    @staticmethod
    def read_generate_file(file_name):
        if os.path.basename(file_name).split('.')[-1] == 'pic':
            with open(file_name, 'rb') as picFile:
                obj = pic.load(picFile)
        else:
            obj = igraph.read(file_name)
        return obj

    @staticmethod
    def make_dirs(working_dir):
        dirs_dic = dict(
            genomes=os.path.join(working_dir, "GenomeSeq"),
            homology=os.path.join(working_dir, "BlastResults"),
            sim_matrix=os.path.join(working_dir, "matrices"),
            refined=os.path.join(working_dir, "refined")
        )
        [os.makedirs(d, exist_ok=True) for d in dirs_dic.values()]

    @staticmethod
    @time_used('Genomes Loading')
    def recode_genome_seq(genome_path, work_dir, suffix):
        genome_files = get_input_genome_inf(genome_path, suffix)
        genome_recode_path = os.path.join(work_dir, 'GenomeSeq')
        seq_inf_pair = {
            'GenomeRecodeInf': {}, 'GenesRecode': {}, 'SeqLengthInf': {}, 'SpeciesGeneNum': {}, 'GenomeToUsed': [],
            'SequencesRecode': {}}
        seq_recode = ''
        for GenomeNum, Genome in enumerate(genome_files):
            seq_inf_pair['GenomeRecodeInf']['G%d' % GenomeNum] = Genome
            fasta_genome = os.path.join(genome_path, Genome)
            fasta_file = parse_fasta(fasta_genome)
            seq_numbers = len(list(fasta_file.keys()))
            seq_inf_pair['GenomeToUsed'].append('G%d' % GenomeNum)
            seq_inf_pair['SpeciesGeneNum']['G%d' % GenomeNum] = seq_numbers
            new_seq = ''
            for GeneNum, SeqID in enumerate(fasta_file):
                new_seq += '%s\n%s\n' % ('>G%d|g%d' % (GenomeNum, GeneNum), fasta_file[SeqID])
                seq_inf_pair['GenesRecode']['G%d|g%d' % (GenomeNum, GeneNum)] = SeqID.split(' ')[0]
                seq_recode += '%s\tG%d|g%d\n' % (SeqID.split(' ')[0], GenomeNum, GeneNum)
                seq_inf_pair['SeqLengthInf']['G%d|g%d' % (GenomeNum, GeneNum)] = len(fasta_file[SeqID][:])
                seq_inf_pair['SequencesRecode'][SeqID.split(' ')[0]] = fasta_file[SeqID][:]
            output_file(os.path.join(genome_recode_path, 'G%d.fa' % GenomeNum), new_seq)
        output_file(os.path.join(work_dir, 'SequenceIDs.txt'), seq_recode)
        [os.remove(os.path.join(genome_path, file)) for file in os.listdir(genome_path) if
         file.split('.')[-1] in ['flat', 'gdx']]
        matrices_dumpy(seq_inf_pair, 'SeqInf', 'Information', work_dir)
        return seq_inf_pair

    @staticmethod
    @time_used('Homology Searching')
    def sequences_homology_search(igp, sm, wd, e, st):
        """
        :param igp: input genome directory
        :param sm: homology searching method
        :param wd: working directory
        :param e: threads for e-value
        :param st: numbers of threads for searching
        :return: None
        """
        input_genomes = os.listdir(igp)
        if sm == 'mmseqs':
            db_s = [os.path.join(igp, genomes) for genomes in os.listdir(igp)]
        else:
            db_s = make_blast_db(input_genomes, igp, sm, wd)
        blast_results_path = os.path.join(wd, 'BlastResults')
        run_blast_search_parallel(
            input_genomes, igp, sm, db_s, blast_results_path, e, st)

    @staticmethod
    @time_used('SSN Building')
    def ssn_construction(md, si, wd, nt):
        """
        :param md: define the most distance in homology
        :param si: sequences information object
        :param wd: working directory
        :param nt: number of threads for network building
        :return: ssn_r sequences similarity network
        """
        matrices_dir = os.path.join(wd, "matrices")
        brp = os.path.join(wd, "BlastResults")
        get_bh_matrix_parallel(brp, md, si, matrices_dir, nt)
        get_connected_matrix_parallel(si, md, matrices_dir, nt)
        ssn_r = build_ssn_parallel(matrices_dir, si, nt, wd)
        return ssn_r

    @staticmethod
    @time_used('HNN Construction')
    def ConstructedHnCC(ssn, dm, w, nt, gc, wd):
        """
        :param ssn: sequences similarity network
        :param dm: community detection method
        :param w: attribution for edges weight
        :param nt: numbers of thread for community detection
        :param gc: gamma coefficient
        :param wd: working directory
        :return: hierarchical network of genes families
        """
        hnn = HHN(ssn)
        connected_components = hnn.splitConnectedComponents()
        turn = 1
        node_list_all = []
        edge_list_all = []
        genes_ids_list_all = []
        genome_ids_list_all = []
        genes_num_list_all = []
        genome_num_list_all = []
        modularity_list_all = []
        community_all, continued = SplitTasks1(connected_components)
        ssn.vs['index'] = ssn.vs.indices
        while community_all:
            my_queues = mp.Manager().Queue()
            if len(community_all) <= int(nt):
                parallel_num = len(community_all)
            else:
                parallel_num = nt
            pro = mp.Pool(processes=int(parallel_num))
            for cc in community_all:
                if continued is None:
                    pro.apply_async(
                        func=get_connected_of_ccs, args=(my_queues, ssn, cc, dm, w, gc))
                else:
                    pro.apply_async(
                        func=GetConnectedOfCCsV2, args=(my_queues, ssn, cc, dm, w, gc))
            connects, second_list = GetConnectedList(my_queues, len(community_all), turn)
            if second_list:
                max_cc = sorted(second_list, key=lambda X: len(X[1]), reverse=True)[0][1]
                print('Max Nodes Numbers of CC %d ' % len(max_cc))
                if len(max_cc) < 5000:
                    community_all, continued = SplitTasks2(second_list)
                else:
                    community_all, continued = SplitTasks1(second_list)
            else:
                community_all = second_list
            pro.close()
            pro.join()
            node_list, edge_list, genes_ids_list, genome_ids_list, genes_num_list, genome_num_list, modularity_list = connects
            node_list_all.extend(node_list)
            edge_list_all.extend(edge_list)
            genes_ids_list_all.extend(genes_ids_list)
            genome_ids_list_all.extend(genome_ids_list)
            genes_num_list_all.extend(genes_num_list)
            genome_num_list_all.extend(genome_num_list)
            modularity_list_all.extend(modularity_list)
            turn += 1
        hhm = igraph.Graph()
        attributes_dic = dict(
            geneIDs=genes_ids_list_all, genomeIDs=genome_ids_list_all, genesNum=genes_num_list_all,
            genomesNum=genome_num_list_all,
            modularity=modularity_list_all)
        hhm.add_vertices(node_list_all, attributes=attributes_dic)
        hhm.add_edges(edge_list_all)
        hhm.write_gml(os.path.join(wd, 'hm.gml'))
        return hhm

    @staticmethod
    @time_used('HNN analysis and OGs generation')
    def hnn_analysis(ssn, hg, si, op, r, nt, wd):
        """
        :param ssn: sequences similarity network
        :param hg: hierarchical graph
        :param si: sequences information
        :param op: species overlap
        :param r: refined
        :param nt: numbers of thread
        :param wd: working directory
        """
        hm = HHN(ssn, hg)
        hhng = hm.GetEvolutionEvents(op, wd)
        refind_dir = os.path.join(wd, 'refined')
        if r:
            hhnR = RefinedNodesEvents(ssn, hhng, si, refind_dir, nt)
            hhnR.write_gml(os.path.join(wd, 'hmmR.gml'), ids="no-such-attribute")
            hhngO = hhnR.copy()
            hhng = hhnR
        else:
            hhng = hhng
            hhng.write_gml(os.path.join(wd, 'hmmO.gml'))
            hhngO = hhng.copy()
        hhng.vs['Event'] = list(map(str, hhng.vs['Event']))
        hhngO.vs['Event'] = list(map(str, hhng.vs['Event']))
        hhngF, OGDic = ExtractOGSorted(hhng, ssn, si, 'I', 0, 0)
        node_recode_dic, OGForPairwise = WriteOGFiles(si, hhngO, OGDic, wd)
        OutputGraph = os.path.join(wd, 'hmm.gml')
        output_final_graph(si, hhngF, OutputGraph, node_recode_dic)
        GetOrthologsFromOGsPP(si, OGForPairwise, wd, nt)


def processes():
    pro_l = np.array(['recode_genome_seq', 'sequences_homology_search', 'ssn_construction', 'ConstructedHnCC',
                      'hnn_analysis'])
    pro_args = [('InputDir', 'WorkingDir', 'GenomeExtension'),
                ('recodeDir', 'seqSearchingMethod', 'WorkingDir', 'E_values', 'threads'),
                ('MostDistance', 'return_obj[0]', 'WorkingDir', 'NetworkThreads'),
                ('return_obj[2]', 'detectionMethod', 'weight', 'NetworkThreads', 'GammaCo', 'WorkingDir'),
                ('return_obj[2]', 'return_obj[3]', 'return_obj[0]', 'Overlaps', 'refined', 'NetworkThreads',
                 'WorkingDir')]
    return pro_l, pro_args


def continued_method(pro_list, args_list, ssn_f, hhn_f, seq_f):
    """
    :param pro_list:
    :param args_list:
    :param ssn_f:
    :param hhn_f:
    :param seq_f:
    :return:
    """
    continued_rank = list(map(os.path.exists, [seq_f, hhn_f, ssn_f]))
    if continued_rank[0] and continued_rank[1]:
        pro_list[[0, 2, 3]] = ['read_generate_file'] * 3
        pro_list[1] = None
        args_list[:4] = ['seq_inf_file', None, 'ssn_file', 'hhn_file']
    elif continued_rank[0] and continued_rank[2]:
        pro_list[1] = None
        pro_list[[0, 2]] = ['read_generate_file'] * 2
        args_list[:3] = ['seq_inf_file', None, 'ssn_file']
    return pro_list, args_list


def main(arguments):
    InputDir = arguments['input_dir']
    GenomeExtension = arguments['extension']
    continuedMethod = arguments['continued']
    E_values = arguments['e_values']
    threads = arguments['search_threads']
    seqSearchingMethod = arguments['search_method']
    outDir = arguments['output_dir']
    weight = arguments['weight_types']
    if weight not in ['NBS', 'Sim', 'BS']:
        print('Error! Attribute does not exist :%s\nCurrent version only provide Sim ,BS , NBS' % weight)
        sys.exit(0)
    NetworkThreads = arguments['graph_threads']
    detectionMethod = arguments['community_detect_method']
    MostDistance = arguments['distance']
    Overlaps = arguments['species_overlap']
    WorkingDir = os.path.join(outDir, 'WorkingDirectory')
    recodeDir = os.path.join(WorkingDir, 'GenomeSeq')
    GammaCo = arguments['gamma_coefficient']
    parameter = ' '.join(sys.argv)
    refined = arguments['refined']
    call_pro_list, call_args_list = processes()
    if continuedMethod:
        ssn_file = os.path.join(WorkingDir, 'ssn.gml')
        hhn_file = os.path.join(WorkingDir, 'hm.gml')
        seq_inf_file = os.path.join(WorkingDir, 'SeqInfInformation.pic')
        call_pro_list, call_args_list = continued_method(call_pro_list, call_args_list, ssn_file, hhn_file,
                                                         seq_inf_file)
    process_call = Call()
    process_call.make_dirs(working_dir=WorkingDir)
    return_obj = []
    for pro, arg in zip(call_pro_list, call_args_list):
        if pro is not None and pro != 'None':
            if isinstance(arg, tuple):
                pro_arg = tuple(map(eval, arg))
                pro_return_obj = methodcaller(pro, *pro_arg)(process_call)
            else:
                pro_arg = eval(arg)
                pro_return_obj = methodcaller(pro, pro_arg)(process_call)
        else:
            pro_return_obj = None
        return_obj.append(pro_return_obj)


if __name__ == '__main__':
    parameters = get_parameters()
    main(parameters)
