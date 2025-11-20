import gzip
import os
from ctypes import c_int64
from datetime import datetime
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value
from uuid import uuid4

import numpy
from threadpoolctl import threadpool_limits

from .bam import decodeRead, indexBam, readIndices
from .processbar import ProcessBar


def isGzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def readAnnotationFile(file):
    cluster2i = dict()
    classifications = list()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Cluster ID\tClassification', f'\"{file}\" is not a valid annotation file.'
    for i, line in enumerate(openFile):
        lines = line.rstrip('\n').split('\t')
        cluster2i[lines[0]] = i
        classifications.append(lines[1])
    openFile.close()
    return (cluster2i, classifications)


def readMappingFile(file, cluster2i):
    sequence2i = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Sequence ID\tCluster ID', f'\"{file}\" is not a valid mapping file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        sequence2i[lines[0]] = cluster2i[lines[1]]
    openFile.close()
    return sequence2i


def readAbundanceFile(file, classification2i, individuals):
    individual2i = dict((j, i) for i, j in enumerate(individuals))
    x = numpy.zeros(shape = (len(classification2i), len(individuals)), dtype = numpy.float32)
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    headers = openFile.readline().rstrip('\n').split('\t')
    assert headers[0] == 'Abundance', f'\"{file}\" is not a valid abundance file.'
    y = [individual2i[i] for i in headers[1 : ]]
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[0] in classification2i:
            x[classification2i[lines[0]], y] = lines[1 : ]
    openFile.close()
    return x == 0


def createProcesses(inputFiles, threads, mapq, identity, coverage, nReads, n):
    queue = Queue(threads)
    processes = list()
    outputFiles = list()
    processBar = ProcessBar(n)
    n = Value(c_int64, 0)
    for i in range(threads):
        outputFiles.append(uuid4().hex)
        processes.append(Process(target = countACGT, args = (queue, inputFiles, outputFiles[-1], mapq, identity, coverage, nReads, n, processBar)))
        processes[-1].start()
    return (queue, processes, outputFiles, n)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None, None, None, None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


@threadpool_limits.wrap(limits = 1)
def countACGT(queue, inputFiles, outputFile, mapq, identity, coverage, nReads, n, processBar):
    nFiles = len(inputFiles)
    tempIndividualACGT = list()
    openFile = open(outputFile, 'w', buffering = 10485760)
    while True:
        sequence, classification, length, fileOffsets, dataOffsets, dataSizes = queue.get()
        if sequence is None:
            break
        x = numpy.zeros(shape = (length, nFiles + 1, 5), dtype = numpy.float32)
        y = numpy.arange(length, dtype = numpy.int64)
        for i, (inputFile, fileOffset, dataOffset, dataSize) in enumerate(zip(inputFiles, fileOffsets, dataOffsets, dataSizes), start = 1):
            if dataSize:
                for read, ungappedRegions in decodeRead(inputFile, fileOffset, dataOffset, dataSize, mapq, identity):
                    # read: 01234 - ACGT? #
                    for sequenceStart, sequenceEnd, readStart, readEnd in ungappedRegions:
                        numpy.add.at(x[ : , i, : ], (y[sequenceStart : sequenceEnd], read[readStart : readEnd]), 1)
        x[ : , : , 4] = 0 # remove ? #
        numpy.einsum('ijk->ik', x[ : , 1 : ], out = x[ : , 0])
        numpy.putmask(x, x < nReads, 0)
        individualN = numpy.sum(x, axis = 2)  # (positions, individuals + 1) #
        majorAllele = numpy.argmax(x[ : , 0, : 4], axis = 1) # (positions, ), major allele has been determined #
        nIndividuals = numpy.count_nonzero(individualN[ : , 1 : ], axis = 1) # (positions, ) #
        y = nIndividuals / (x.shape[1] - 1) >= coverage # (positions, ) #
        individualN += 10 * numpy.finfo(numpy.float32).eps
        x /= individualN[ : , : , None]
        numpy.putmask(y, numpy.count_nonzero(x[numpy.arange(x.shape[0]), 1 : , majorAllele] > 0.05, axis = 1) / (nIndividuals + 10 * numpy.finfo(numpy.float32).eps) > 0.99, False)
        numpy.round(x, 4, out = x)
        for i in numpy.flatnonzero(y):
            populationACGT = '/'.join(str(j) for j in x[i, 0, : 4])
            for xi in x[i, 1 : ]:
                if numpy.max(xi) > 0:
                    tempIndividualACGT.append('/'.join(str(j) for j in xi[ : 4]))
                else:
                    tempIndividualACGT.append('NA')
            tempIndividualACGT_ = '\t'.join(tempIndividualACGT)
            openFile.write(f'{classification}\t{sequence}\t{i + 1}\t{populationACGT}\t{tempIndividualACGT_}\n')
            tempIndividualACGT.clear()
        n.acquire()
        n.value += length
        processBar.plot(n.value)
        n.release()
    openFile.close()
    return None


def writeFile(individuals, inputFiles, outputFile):
    processBar = ProcessBar(len(inputFiles))
    open4w = open(outputFile, 'wb', buffering = 10485760)
    open4w.write(f'Classification\tChromosome\tPosition\tPopulation\t{individuals}\n'.encode('utf8'))
    for i, inputFile in enumerate(inputFiles, start = 1):
        open4r = open(inputFile, 'rb', buffering = 10485760)
        while open4w.write(open4r.read(10485760)):
            pass
        open4r.close()
        os.remove(inputFile)
        processBar.plot(i)
    open4w.close()
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Indexing all bam files.', flush = True)
    indexBam(parameters.bam, parameters.threads_index)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading all index files.', flush = True)
    sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndices(parameters.bam, parameters.threads_index)

    individuals = [os.path.splitext(os.path.basename(i))[0] for i in parameters.bam]

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading annotation file.', flush = True)
    cluster2i, classifications = readAnnotationFile(parameters.annotation)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading mapping file.', flush = True)
    sequence2i = readMappingFile(parameters.mapping, cluster2i)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Calling variants.', flush = True)
    queue, processes, tempfiles, n = createProcesses(parameters.bam, parameters.threads_call, parameters.mapq, parameters.identity, parameters.coverage, parameters.depth, numpy.sum(lengths))
    for sequence, length, fileOffset, dataOffset, dataSize in zip(sequences, lengths, fileOffsets, dataOffsets, dataSizes):
        i = sequence2i[sequence]
        queue.put((sequence, classifications[i], length, fileOffset, dataOffset, dataSize))
    freeProcesses(queue, processes)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing to \"{parameters.output}\".', flush = True)
    writeFile('\t'.join(individuals), tempfiles, parameters.output)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
