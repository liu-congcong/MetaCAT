import os
from datetime import datetime
from math import ceil, log10

from .fasta import readFastaFile


def main(parameters):
    if not parameters.output.endswith('.mapping'):
        parameters.output = parameters.output + '.mapping'
    clusterWidth = ceil(log10(len(parameters.fasta) + 1))
    clusterPrefix = os.path.basename(parameters.output)[ : -8]
    x = list()
    i = 1
    openFile = open(parameters.output, 'w')
    openFile.write('Sequence ID\tCluster ID\n')
    for file in parameters.fasta:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading \"{file}\".', flush = True)
        clusterSize = 0
        for sequenceID, sequence in readFastaFile(file):
            x.append(f'{sequenceID}\t{clusterPrefix}.{i:0{clusterWidth}d}\n')
            clusterSize += len(sequence)
        if clusterSize >= parameters.min_cluster:
            openFile.writelines(x)
            i += 1
        x.clear()
    openFile.close()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
