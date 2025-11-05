import gzip
import os
from datetime import datetime

from .checkm2 import readCheckm2File
from .fasta import readFastaFile


def isGzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def readMappingFile(file):
    clusterID2sequenceIDs = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, mode = 'rt')
    assert openFile.readline().rstrip('\n') in ('Sequence ID\tCluster ID', 'Sequence ID\tGenome ID\tLength'), f'\"{file}\" must have a single header line: Sequence ID<tab>Cluster ID.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        clusterID2sequenceIDs.setdefault(lines[1], list()).append(lines[0])
    openFile.close()
    return clusterID2sequenceIDs


def writeFastaFiles(sequenceID2sequence, sequenceID2length, clusterID2sequenceIDs, minSize, output, lineLength = 100):
    for clusterID, sequenceIDs in clusterID2sequenceIDs.items():
        if sum(sequenceID2length[sequenceID] for sequenceID in sequenceIDs) >= minSize:
            file = os.path.join(output, f'{clusterID}.fasta')
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing to {file}.', flush = True)
            openFile = open(file, 'w')
            for sequenceID in sequenceIDs:
                sequence = sequenceID2sequence[sequenceID]
                openFile.write('>' + sequenceID + '\n')
                filePointer = 0
                while openFile.write(sequence[filePointer : filePointer + lineLength]):
                    openFile.write('\n')
                    filePointer += lineLength
            openFile.close()
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading assembly file.', flush = True)
    sequenceID2sequence = dict()
    sequenceID2length = dict()
    for sequenceID, sequence in readFastaFile(parameters.fasta):
        sequenceID2sequence[sequenceID] = sequence
        sequenceID2length[sequenceID] = len(sequence)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading mapping file.', flush = True)
    clusterID2sequenceIDs_ = readMappingFile(parameters.mapping)
    clusterID2sequenceIDs = dict()
    if parameters.checkm2 is not None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading CheckM2 file.', flush = True)
        clusters = readCheckm2File(parameters.checkm2, parameters.contamination, parameters.completeness)
        for i, j in clusterID2sequenceIDs_.items():
            if i in clusters:
                clusterID2sequenceIDs[i] = j
    else:
        clusterID2sequenceIDs = clusterID2sequenceIDs_
    writeFastaFiles(sequenceID2sequence, sequenceID2length, clusterID2sequenceIDs, parameters.min_cluster, parameters.output)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
