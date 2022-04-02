import os
import sys
import subprocess
import multiprocessing as mp
import time
import argparse


def get_parameters():
    parse = argparse.ArgumentParser()
    parse.add_argument('-q', required=True, action='store', help='Query sequences directory')
    parse.add_argument('-o', required=True, action='store', help='blast out directory')
    parse.add_argument('-t', required=False, action='store', help='thread numbers default: 8', default='8')
    parse.add_argument('--db', required=False, action='store', help='database path default: /home/licongjian/COGs-new',
                       default='/home/licongjian/COGs-new')
    parse.add_argument('-e', required=False, action='store', help='blastp E-value default:0.001', default='0.001')
    parse.add_argument('-c', required=False, action='store', help='similarity threads default:20', default='20')
    parse.add_argument('-s', required=False, action='store', help='coverage threads default:40', default='40')
    arguments = parse.parse_args()
    return arguments


def input_file(file_name):
    with open(file_name) as f:
        data = f.readlines()
    return data


def output_file(file_name, content):
    with open(file_name, 'w') as f:
        f.writelines(content)


def ReadFile(FileName):
    with open(FileName) as f:
        data = f.read()
    return data


def ReadFastaFile(Fasta):
    seqDic = {}
    seqString = ''
    for line in ReadFile(Fasta).split('\n'):
        if line.startswith('>'):
            seqString += '\n' + line + '\n'
        else:
            seqString += line.strip()
    seqList = seqString.strip('\n').split('\n')
    for index in range(int(len(seqList) / 2)):
        seqDic[seqList[2 * index]] = seqList[2 * index + 1]
    return seqDic


def processing(ntd=0, ttn=100, scale=0.4, message='', symbol='='):
    """
    Print process bar on screen
    :param ntd: Number of Tasks Done
    :param ttn: Total Task Number
    :param scale: arrows ('>') length scaling
    :param message: extra information
    :param symbol: symbol to print
    :return:
    """
    h = int(ntd * 100 / ttn)  # percentage of number task done
    i = int(h * scale)  # length of arrow ('>')
    j = int(100 * scale - i)  # length of spaces
    arrows = symbol * i + ' ' * j
    sys.stdout.write('\r%s' % message + arrows + '[%s%%]' % h)
    sys.stdout.flush()


def make_search_db(references_seq, db_dir, threads):
    references_seq_prefix = os.path.abspath(references_seq).split('/')[-1].split('.')[0]
    db = os.path.join(db_dir, 'db%s' % references_seq_prefix)
    MakeDBCmd = ' '.join(
        ['diamond', 'makedb', '--in', references_seq, '--db', db, '--threads', str(threads), '>', '/dev/null'])
    make_db = subprocess.Popen(MakeDBCmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    make_db.wait()
    return db


def search_sequences(queue, query_file, db_name, out_put_dir, E):
    BlastFilePath = os.path.join(out_put_dir, '%s.out' % os.path.abspath(query_file).split('/')[-1])
    command = ' '.join([
        'diamond', 'blastp', '--more-sensitive', '-p', '1', '-q', query_file, '-d', '%s.dmnd' % db_name,
        '--evalue', str(E), '-f', '6', '--out', BlastFilePath, '--quiet'])
    search = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    search.wait()
    if search.returncode == 0:
        os.system('echo END >> %s' % BlastFilePath)
    queue.put(search.returncode)


def recall_results(queues, done_task, total_task):
    while True:
        Time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
        processing(done_task, message='[%s]Sequences Searching' % Time)
        statue = queues.get()
        if statue == 0:
            pass
        done_task += 1
        if done_task >= total_task:
            break


def search_sequences_pp(done_numbers, total_numbers, query_genes_list, db_name, searching_threads, output_dir, e_value):
    process = mp.Pool(int(searching_threads))
    my_queue = mp.Manager().Queue()
    for query in query_genes_list:
        process.apply_async(func=search_sequences, args=(my_queue, query, db_name, output_dir, e_value))
    recall_results(my_queue, done_numbers, total_numbers)
    process.close()
    process.join()


def recode_cog_db(cog_db_dir, db_dir, seqInf):
    recode_db_file = ''
    for cog_file in os.listdir(cog_db_dir):
        cog = ReadFastaFile(os.path.join(cog_db_dir, cog_file))
        for locus, seq in cog.items():
            recode_locus = '%s_%s' % (locus.strip(), cog_file)
            seqInf[recode_locus.strip('>').split(' ')[0]] = seq
            recode_db_file += '%s\n%s\n' % (recode_locus, seq)
    print(os.path.join(db_dir, 'cog_db.fasta'))
    output_file(os.path.join(db_dir, 'cog_db.fasta'), recode_db_file)


def handle_search_results(search_out_dir, sequence_inf, sim_threads, cov_threads):
    cov_threads = cov_threads * 100 if cov_threads <= 1 else cov_threads
    sim_threads = sim_threads * 100 if sim_threads <= 1 else sim_threads
    search_dic = {}
    for search_out in os.listdir(search_out_dir):
        search_file = os.path.join(search_out_dir, search_out)
        for line in input_file(search_file):
            if 'END' not in line:
                query_gene = line.strip().split('\t')[0]
                refer_gene = line.strip().split('\t')[1]
                cog = refer_gene.split('_')[-1]
                sim = line.strip().split('\t')[2]
                Q_Align_length = abs(int(line.strip().split('\t')[6]) - int(line.strip().split('\t')[7]))
                R_Align_length = abs(int(line.strip().split('\t')[8]) - int(line.strip().split('\t')[9]))
                refer_coverage = R_Align_length / len(sequence_inf[refer_gene]) * 100
                queries_coverage = Q_Align_length / len(sequence_inf[query_gene]) * 100
                coverage = queries_coverage if queries_coverage < refer_coverage else refer_coverage
                line_inf = line.strip().split('\t')
                line_inf.append(cog)
                if float(sim) >= sim_threads and coverage >= cov_threads:
                    if query_gene not in search_dic.keys():
                        search_dic[query_gene] = [line_inf]
                    else:
                        search_dic[query_gene].append(line_inf)
    print(list(search_dic.values())[0])
    return search_dic


def annotation(search_results, output_dir):
    annotation_result = ''
    for gene, hit_genes in search_results.items():
        best_hit = sorted(hit_genes, key=lambda X: (float(X[2]), -float(X[-3]), float(X[-2])), reverse=True)[0]
        if gene == best_hit[0]:
            annotation_result += '\t'.join(best_hit) + '\n'
    output_file(os.path.join(output_dir, 'cog_annotation.txt'), annotation_result)


def Get_sequences_inf(genes_query_list):
    seq_inf = {}
    for file in genes_query_list:
        seq_fasta = ReadFastaFile(file)
        for locus, seq in seq_fasta.items():
            seq_inf[locus.split(' ')[0].strip('>')] = seq
    return seq_inf


def read_genome_files(genome_dir):
    genome_list = []
    for genome_file in os.listdir(genome_dir):
        if genome_file.split('.')[-1] in ['faa', 'fna', 'fasta', 'fas']:
            genome_list.append(os.path.join(genome_dir, genome_file))
    return genome_list


def check_point(genome_list, blast_out_dir):
    genome_done = 0
    genome_done_list = []
    for blast_file in os.listdir(blast_out_dir):
        if 'END' in input_file(os.path.join(blast_out_dir, blast_file))[-1]:
            genome = blast_file.split('/')[-1].split('.out')[0]
            genome_done_list.append(genome)
            genome_done += 1
    genome_continued = []
    for genome_file in genome_list:
        genome_ = genome_file.split('/')[-1]
        if genome_ not in genome_done_list:
            genome_continued.append(genome_file)
    return genome_done, genome_continued


if __name__ == '__main__':
    args = get_parameters()
    query_genome = args.q
    refer_path = args.db
    output_directory = args.o
    pro_threads = int(args.t)
    cov_ = float(args.c)
    sim_ = float(args.s)
    e_ = float(args.e)
    working_dir = os.path.join(output_directory, 'working_dir')
    blast_dir = os.path.join(working_dir, 'blast_out')
    for dir_ in [output_directory, working_dir, blast_dir]:
        try:
            os.mkdir(dir_)
        except FileExistsError:
            pass
    genome_files = read_genome_files(query_genome)
    seq_information = Get_sequences_inf(genome_files)
    recode_cog_db(refer_path, working_dir, seq_information)
    if os.listdir(blast_dir):
        genome_done_numbers, genome_continued_list = check_point(genome_files, blast_dir)
    else:
        genome_done_numbers, genome_continued_list = 0, genome_files
    if genome_continued_list:
        db_seq = os.path.join(working_dir, 'cog_db.fasta')
        database = make_search_db(db_seq, working_dir, pro_threads)
        search_sequences_pp(genome_done_numbers, len(genome_files), genome_continued_list, database, pro_threads,
                            blast_dir, e_)
    search_results_dic = handle_search_results(blast_dir, seq_information, sim_, cov_)
    annotation(search_results_dic, output_directory)
