import igraph
import sys
import os
import pickle as pic
import shutil
import numpy as np
from collections import Counter


def MatricesLoad(matrixType, name, path):
    with open(os.path.join(path, '%s%s.pic' % (matrixType, name)), 'rb') as picFile:
        M = pic.load(picFile)
    return M


def input_file(inputName):
    with open(inputName) as f:
        data = f.read()
    return data


def ReadFile(FileName):
    with open(FileName) as file:
        data = file.readlines()
    return data


def OutputFile(fileName, content):
    with open(fileName, 'w') as f:
        f.writelines(content)


def Copy_og_file(from_, to_):
    if os.path.isfile(from_):
        shutil.copy(from_, to_)
    else:
        print('%s File Not Found' % from_)
        sys.exit(0)


def transfer_cc_to_trees(cc_subgraph, cc_numbers, vs_list, genome_numbers, genes_numbers, edges, edge_weight,
                         cog_list, cog_c_list, cog_gene_list, pathway_list, mapping):
    try:
        rooted_nodes = cc_subgraph.vs.select(_degree=2)[0]
    except IndexError:
        cc_subgraph.write_gml('test.gml')
        print(cc_subgraph, cc_numbers)
    rooted_nodes = cc_subgraph.vs.select(_degree=2)[0]
    rooted_nodes_index = rooted_nodes.index
    vs_tuple = cc_subgraph.bfs(vid=rooted_nodes_index)
    vs_list.append(str(cc_numbers))
    genome_numbers.append(rooted_nodes['genomesNum'])
    genes_numbers.append(rooted_nodes['genesNum'])
    cog_list.append("None")
    cog_c_list.append("None")
    cog_gene_list.append("None")
    pathway_list.append("None")
    for vs_index in vs_tuple[0]:
        if vs_index != rooted_nodes_index:
            nodes = cc_subgraph.vs[vs_index]
            if nodes['name'] == 'None' or nodes['name'] is None or int(nodes['genesNum']) == 1:
                pass
            else:
                shortest_path = rooted_nodes.get_shortest_paths(nodes)[0]
                shortest_path_length = len(shortest_path[1:])
                edge_weight.append(int(shortest_path_length))
                vs_list.append(nodes['name'])
                try:
                    cog_ID, cog_C, cog_gene, cog_pathway = mapping[nodes['name']]
                except KeyError:
                    cog_ID, cog_C, cog_gene, cog_pathway = ['None'] * 4
                cog_list.append(cog_ID)
                cog_c_list.append(cog_C)
                cog_gene_list.append(cog_gene)
                pathway_list.append(cog_pathway)
                edges.append([cc_numbers, nodes['name']])
                genome_numbers.append(nodes['genomesNum'])
                genes_numbers.append(nodes['genesNum'])
    graph_list = [vs_list, genome_numbers, genes_numbers, edges, edge_weight]
    return graph_list


def Get_cc_og_inf(cc_graph, cc_numbers, output_Dir, og_dir, numbers_threads):
    og_numbers = 0
    core_og_numbers = 0
    single_og_numbers = 0
    nodes = cc_graph.vs.select(name_ne="None", genesNum_gt=1)
    cc_dir = os.path.join(output_Dir, cc_numbers)
    try:
        os.mkdir(cc_dir)
    except FileExistsError:
        pass
    cc_inf = ''
    og_seq = ''
    for node in nodes:
        if 'OG' in node['name']:
            og_file_from = os.path.join(og_dir, '%s.fasta' % node['name'])
            og_seq += input_file(og_file_from) + '\n'
            og_file_to = os.path.join(cc_dir, '%s.fasta' % node['name'])
            Copy_og_file(og_file_from, og_file_to)
            og_numbers += 1
            cc_inf += '%s\t%d\t%d\n' % (node['name'], int(node['genesNum']), int(node['genomesNum']))
            if int(node['genomesNum']) >= numbers_threads:
                core_og_numbers += 1
                if int(node['genesNum']) == numbers_threads:
                    single_og_numbers += 1
        else:
            print(node['name'])
    OutputFile(os.path.join(cc_dir, '%s_og_inf.txt' % cc_numbers), cc_inf)
    return og_numbers, core_og_numbers, single_og_numbers, og_seq


def pangenome_structure(genomeNum, whole_genome_num, genesNum):
    proportion = genomeNum / whole_genome_num
    if 0 <= proportion < 0.15:
        pan = 'cloud'
    elif 0.15 <= proportion < 0.95:
        pan = 'shell'
    elif 0.95 <= proportion < 0.99:
        pan = 'soft core'
    elif 0.99 <= proportion < 1:
        pan = 'core'
    else:
        if genomeNum == genesNum:
            pan = 'SOG'
        else:
            pan = 'HOG'
    return pan


def Statistic_graph_cc(whole_graph, output_dir, ogs_dir, genome_numbers, mapping_inf, ortho_seq_dir):
    plot_graph = igraph.Graph()
    while_graph_gml = igraph.read(whole_graph)
    og_seq_for_annotation = ''
    name_list = []
    genome_list = []
    gene_list = []
    edge_list = []
    edge_weight = []
    cog_id = []
    cog_cate = []
    cog_gene = []
    cog_path = []
    for num, sub_cc in enumerate(while_graph_gml.components()):
        subgraph = while_graph_gml.subgraph(sub_cc)
        CC_Numbers = '{}{:0>5d}'.format('CC', num)
        if len(sub_cc) > 1:
            transfer_cc_to_trees(subgraph, CC_Numbers, name_list, genome_list, gene_list, edge_list, edge_weight,
                                 cog_id, cog_cate, cog_gene, cog_path, mapping_inf)
        else:
            og_node = subgraph.vs[0]
            if og_node['genesNum'] > 1:
                try:
                    cog_single, cong_c_s, cog_g_s, cog_p_s = mapping_inf[og_node['name']]
                except KeyError:
                    cog_single, cong_c_s, cog_g_s, cog_p_s = ["None"] * 4
                cog_id.append(cog_single)
                cog_cate.append(cong_c_s)
                cog_gene.append(cog_g_s)
                cog_path.append(cog_p_s)
                name_list.append(og_node['name'])
                genome_list.append(og_node['genomesNum'])
                gene_list.append(og_node['genesNum'])
            else:
                pass
    nodes_attribution = dict(genomesNum=genome_list, genesNum=gene_list, COG_ID=cog_id, COG_Cate=cog_cate,
                             COG_Genes=cog_gene, COG_Pathway=cog_path)
    edges_attribution = dict(rooted_distances=edge_weight)
    plot_graph.add_vertices(name_list, attributes=nodes_attribution)
    plot_graph.add_edges(edge_list, attributes=edges_attribution)
    OutputFile(os.path.join(output_dir, 'sequences_for_annotation.fasta'), og_seq_for_annotation)
    split_plot_graph(plot_graph, output_dir, ortho_seq_dir)
    plot_graph.write_gml(os.path.join(output_dir, 'radiation_graph.gml'))


def CalshannoEnt(data):
    """
    :param data: label list
    :return: shnno ent
    """
    label_stat = Counter(data)
    shanno_ent_list = []
    for label, numbers in label_stat.items():
        prob = numbers / len(data)
        shanno = prob * np.log2(prob)
        shanno_ent_list.append(shanno)
    shanno_ent = -1 * np.sum(shanno_ent_list)
    return shanno_ent


def CC_ShannoEnt(cc_graph):
    labels_list = cc_graph.vs['COG_Cate']
    shanno_ent_cc = CalshannoEnt(labels_list)
    return shanno_ent_cc


def find_ogs(cc_graphs, orthology_dir, sog_dir, find_type):
    if find_type == 'sog':
        core_genes = cc_graphs.vs.select(genesNum_eq=109, genomesNum_eq=109)
    else:
        core_genes = cc_graphs.vs.select(genesNum_gt=109, genomesNum_eq=109)
    og_name_list = []
    rooted_sum_numbers = 0
    if len(core_genes) > 0:
        single_genes_nodes = [index for index in core_genes.indices if 'CC' not in cc_graphs.vs[index]['name']]
        sog_numbers = len(single_genes_nodes)
        if len(single_genes_nodes) > 1:
            edge_select = cc_graphs.es.select(_target=single_genes_nodes)
            rooted_list = edge_select['rooted_distances']
            rooted_sum = np.sum(rooted_list)
            rooted_sum_numbers = rooted_sum
            if rooted_sum >= 200:
                min_rooted_es = sorted(edge_select, key=lambda X:X['rooted_distances'])[0]
                node_name = min_rooted_es.target_vertex['name']
                og_name_list.append(node_name)
            else:
                for e_og in edge_select:
                    node_name = e_og.target_vertex['name']
                    og_name_list.append(node_name)
        elif len(single_genes_nodes) == 1:
            node_name = cc_graphs.vs[single_genes_nodes[0]]['name']
            og_name_list.append(node_name)
        else:
            pass
    else:
        sog_numbers = 0
    for og in og_name_list:
        og_file = os.path.join(orthology_dir, '%s.fasta' % og)
        to = os.path.join(sog_dir, '%s.fasta' % og)
        Copy_og_file(og_file, to)
    return len(og_name_list), sog_numbers, rooted_sum_numbers


def split_plot_graph(plot_whole_graph, out_dir, ortho_dir):
    single_10 = []
    small_down_10 = []
    medium_10_50 = []
    large_50_200 = []
    ultra_above = []
    sub_link = []
    sub_unlike = []
    shnno = ''
    sog = ''
    sog_dir = os.path.join(out_dir, 'sog')
    hog_dir = os.path.join(out_dir, 'hog')
    try:
        os.mkdir(sog_dir)
    except FileExistsError:
        pass
    try:
        os.mkdir(hog_dir)
    except FileExistsError:
        pass
    for cc in plot_whole_graph.components():
        sub = plot_whole_graph.subgraph(cc)
        sog_select_og_numbers, sog_og_numbers, sog_root_sum = find_ogs(sub, ortho_dir, sog_dir, 'sog')
        hog_select_og_numbers, hog_og_numbers, hog_root_sum = find_ogs(sub, ortho_dir, hog_dir, 'hog')
        if sog_og_numbers > 0 or hog_og_numbers > 0:
            sog += '%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (len(cc), sog_og_numbers, sog_select_og_numbers, hog_og_numbers, hog_select_og_numbers, sog_root_sum, hog_root_sum )
        if len(cc) > 2:
            og_name = [names for names in sub.vs['name'] if 'CC' in names][0]
            shannos = CC_ShannoEnt(sub)
            vcount = sub.vcount()
            shnno += '%s\t%f\t%d\n' % (og_name, shannos, vcount)
        for c in cc:
            if len(cc) <= 2:
                sub_unlike.append(c)
            else:
                sub_link.append(c)
            if len(cc) <= 2:
                single_10.append(c)
            elif 2 < len(cc) <= 10:
                small_down_10.append(c)
            elif 10 < len(cc) <= 50:
                medium_10_50.append(c)
            elif 50 < len(cc) <= 200:
                large_50_200.append(c)
            else:
                ultra_above.append(c)
    single_graph = plot_whole_graph.subgraph(single_10)
    small_down_graph = plot_whole_graph.subgraph(small_down_10)
    medium_graph = plot_whole_graph.subgraph(medium_10_50)
    large_graph = plot_whole_graph.subgraph(large_50_200)
    ultra_graph = plot_whole_graph.subgraph(ultra_above)
    subgraph_link = plot_whole_graph.subgraph(sub_link)
    subgraph_unlink = plot_whole_graph.subgraph(sub_unlike)
    subgraph_link.write_gml(os.path.join(out_dir, 'radiation_subgraph_link.gml'))
    subgraph_unlink.write_gml(os.path.join(out_dir, 'radiation_subgraph_unlink.gml'))
    OutputFile(os.path.join(out_dir, 'shanno_stat.txt'), shnno)
    OutputFile(os.path.join(out_dir, 'sog_select.txt'), sog)
    for n, g in enumerate([single_graph, small_down_graph, medium_graph, large_graph, ultra_graph]):
        g.write_gml(os.path.join(out_dir, 'radiation_subgraph%d.gml' % n))


def get_og_annotation_file(annotate_file):
    og_match = {}
    for line in ReadFile(annotate_file):
        og, cog_id, cog_c, genes, pathway, genesNum = line.strip().split('\t')
        og_name = og.split('.')[0]
        if cog_c != "None":
            cog_c = cog_c[0]
        else:
            cog_c = cog_c
        og_match[og_name] = [cog_id, cog_c, genes, pathway]
    return og_match


if __name__ == '__main__':
    whole_graph_file = sys.argv[1:][0]
    orthogroups_dir = sys.argv[1:][1]
    Out = sys.argv[1:][2]
    genome_all_numbers = sys.argv[1:][3]
    annotation_file = sys.argv[1:][4]
    orthologous_sequences_dir = sys.argv[1:][5]
    cog_mapping = get_og_annotation_file(annotation_file)
    try:
        os.mkdir(Out)
    except FileExistsError:
        pass
    Statistic_graph_cc(whole_graph_file, Out, orthogroups_dir, int(genome_all_numbers), cog_mapping,
                       orthologous_sequences_dir)
