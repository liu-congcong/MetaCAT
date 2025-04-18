import os
from datetime import datetime
from math import ceil, log10

from .fasta import read_fasta_file as readFastaFile


def main(parameters):
    if not parameters.output.endswith('.mapping'):
        parameters.output = parameters.output + '.mapping'
    clusterWidth = ceil(log10(len(parameters.fasta) + 1))
    clusterPrefix = os.path.basename(parameters.output)[ : -8]
    openFile = open(parameters.output, 'w')
    openFile.write('Sequence ID\tCluster ID\n')
    for i, file in enumerate(parameters.fasta, start = 1):
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading \"{file}\".', flush = True)
        for sequence, _ in readFastaFile(file):
            openFile.write(f'{sequence}\t{clusterPrefix}.{i:0{clusterWidth}d}\n')
    openFile.close()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
