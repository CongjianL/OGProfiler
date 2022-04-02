import igraph
import sys


def OutputFile(fileName, content):
    with open(fileName, 'w') as f:
        f.writelines(content)


def ReadFile(FileName):
    with open(FileName) as file:
        data = file.readlines()
    return data


def FindNodes(hhnName, ssnName, statName):
    hhn = igraph.read(hhnName)
    ssn = igraph.read(ssnName)
    WeightStat = ''
    InsperateNodes = hhn.vs.select(Event='None', genesNum_gt=1)
    print(len(InsperateNodes))
    for node in InsperateNodes:
        nodesIDs = node['geneIDs'].strip().split(' ')
        sub = ssn.subgraph(nodesIDs)
        for es in sub.es:
            WeightStat += '%s\t%s\t%s\n' % (es.source_vertex['name'],  es.target_vertex['name'], str(es['Sim']))
    OutputFile(statName, WeightStat)


if __name__ == '__main__':
    HHN = sys.argv[1:][0]
    SSN = sys.argv[1:][1]
    Out = sys.argv[1:][2]
    FindNodes(HHN, SSN, Out)


