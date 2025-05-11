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


def createProcesses(inputFiles, threads, mapq, identity, coverage, nReads, n):
    queue = Queue(threads)
    processes = list()
    outputFiles = list()
    processBar = ProcessBar(n)
    n = Value(c_int64, 0)
    with threadpool_limits(limits = 1):
        for i in range(threads):
            outputFiles.append(uuid4().hex)
            processes.append(Process(target = countACGT, args = (queue, inputFiles, outputFiles[-1], mapq, identity, coverage, nReads, n, processBar)))
            processes[-1].start()
    return (queue, processes, outputFiles, n)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None, None, None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def countACGT(queue, inputFiles, outputFile, mapq, identity, coverage, nReads, n, processBar):
    nFiles = len(inputFiles)
    tempIndividualACGT = list()
    openFile = open(outputFile, 'w')
    while True:
        sequence, length, fileOffsets, dataOffsets, dataSizes = queue.get()
        if sequence is None:
            break
        x = numpy.zeros(shape = (length, nFiles + 1, 5), dtype = numpy.float32)
        for i, (inputFile, fileOffset, dataOffset, dataSize) in enumerate(zip(inputFiles, fileOffsets, dataOffsets, dataSizes), start = 1):
            for read, ungappedRegions in decodeRead(inputFile, fileOffset, dataOffset, dataSize, mapq, identity):
                # read: 01234 - ACGT? #
                for sequenceStart, sequenceEnd, readStart, readEnd in ungappedRegions:
                    x[numpy.arange(sequenceStart, sequenceEnd, dtype = numpy.int64), i, read[readStart : readEnd]] += 1
        for position, acgt in callVariants(x[ : , : , : 4], coverage, nReads):
            for i in acgt[1 : ]:
                if numpy.any(i > 0):
                    tempIndividualACGT.append('/'.join(i.astype(numpy.str_).tolist()))
                else:
                    tempIndividualACGT.append('NA')
            openFile.write('\t'.join([sequence, str(position), '/'.join(acgt[0].astype(numpy.str_).tolist()), '\t'.join(tempIndividualACGT)]) + '\n')
            tempIndividualACGT.clear()
        n.acquire()
        n.value += length
        processBar.plot(n.value)
        n.release()
    openFile.close()
    return None


def callVariants(x, coverage, n):
    x[x < n] = 0 # (positions, individuals + 1, 4) #
    numpy.einsum('ijk->ik', x[ : , 1 : ], out = x[ : , 0]) # x[ : , 0]: (positions, 4) #
    individualN = numpy.sum(x, axis = 2) # (positions, individuals + 1) #
    majorAllele = numpy.argmax(x[ : , 0], axis = 1) # (positions, ), major allele has been determined #
    nIndividuals = numpy.count_nonzero(individualN[ : , 1 : ], axis = 1) # (positions, ) #
    mask = nIndividuals / x.shape[1] >= coverage # positions #
    individualN += 10 * numpy.finfo(numpy.float32).eps
    x /= individualN[ : , : , None]
    mask[numpy.count_nonzero(x[numpy.arange(x.shape[0]), 1 : , majorAllele] > 0.05, axis = 1) / (nIndividuals + 10 * numpy.finfo(numpy.float32).eps) > 0.99] = False # snp has been determined #
    numpy.round(x, 4, out = x)
    for i in numpy.flatnonzero(mask):
        yield (i + 1, x[i]) # x[i]: (individuals + 1, 4) #
    return None


def writeFile(individuals, inputFiles, outputFile):
    individuals = [os.path.splitext(os.path.basename(i))[0] for i in individuals]
    open4w = open(outputFile, 'wb')
    open4w.write(('\t'.join(['Chromosome', 'Position', 'Population'] + individuals) + '\n').encode('utf8'))
    for inputFile in inputFiles:
        open4r = open(inputFile, 'rb')
        while open4w.write(open4r.read(10485760)):
            pass
        open4r.close()
        os.remove(inputFile)
    open4w.close()
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Indexing all bam files.', flush = True)
    indexBam(parameters.input, parameters.threads_index)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading all index files.', flush = True)
    sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndices(parameters.input, parameters.threads_index)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Calling variants.', flush = True)
    queue, processes, tempfiles, n = createProcesses(parameters.input, parameters.threads_call, parameters.mapq, parameters.identity, parameters.coverage, parameters.depth, numpy.sum(lengths))
    for sequence, length, fileOffset, dataOffset, dataSize in zip(sequences, lengths, fileOffsets, dataOffsets, dataSizes):
        queue.put((sequence, length, fileOffset, dataOffset, dataSize))
    freeProcesses(queue, processes)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing to \"{parameters.output}\".', flush = True)
    writeFile(parameters.input, tempfiles, parameters.output)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
