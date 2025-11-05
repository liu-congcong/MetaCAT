from ctypes import c_float, c_int32
from math import log
from multiprocessing import JoinableQueue, Process
from multiprocessing.sharedctypes import RawArray

import numpy
from scipy.sparse import csr_array
from threadpoolctl import threadpool_limits

from .coverage_model import model as coverageModel
from .kmer_frequency_model import model as kmerFrequencyModel

log0 = -1e10


@threadpool_limits.wrap(limits = 1)
def workerProcess(
    processQ, nodes, n1, m1, n2, m2,
    neighborIndices, neighborAffinities,
    nModelW, modelW,
    flag,
    kmerFreq, kmerFreqSquaredRowSum, length, kmerFreqDim, minKmerFreqP,
    covMean, covVar, covDim, minCovP
):
    # shared data #
    nodes = numpy.ndarray(shape = n1, dtype = numpy.int32, buffer = nodes) # 0 <= nodes < n1 #
    neighborX1 = numpy.ndarray(shape = (n1, 1, m1), dtype = numpy.int32, buffer = neighborIndices)
    neighborY1 = numpy.ndarray(shape = (n1, 1, m1), dtype = numpy.float32, buffer = neighborAffinities)
    neighborX2 = numpy.ndarray(shape = (n2, nModelW, m2), dtype = numpy.int32, buffer = neighborIndices)
    neighborY2 = numpy.ndarray(shape = (n2, nModelW, m2), dtype = numpy.float32, buffer = neighborAffinities)
    modelW = numpy.ndarray(shape = (nModelW, 2), dtype = numpy.float32, buffer = modelW)
    kmerFreq = numpy.ndarray(shape = (n1, kmerFreqDim), dtype = numpy.float32, buffer = kmerFreq)
    kmerFreqSquaredRowSum = numpy.ndarray(shape = n1, dtype = numpy.float32, buffer = kmerFreqSquaredRowSum)
    length = numpy.ndarray(shape = (6, n1), dtype = numpy.float32, buffer = length)
    covMean = numpy.ndarray(shape = (n1, covDim), dtype = numpy.float32, buffer = covMean)
    covVar = numpy.ndarray(shape = (n1, covDim), dtype = numpy.float32, buffer = covVar)
    token_ = -1

    while True:
        I = processQ.get()
        if I is None:
            processQ.task_done()
            break
        if token_ != flag[1]:
            token_ = flag[1]
            nodes_ = nodes[ : flag[0]]
            kmerFreq_ = kmerFreq[nodes_]
            kmerFreqSquaredRowSum_ = kmerFreqSquaredRowSum[nodes_]
            length_ = length[ : , nodes_]
            covMean_ = covMean[nodes_]
            covVar_ = covVar[nodes_]
            nNeighbors1 = min(m1, flag[0])
            nNeighbors2 = min(m2, flag[0])
            temp = numpy.empty(shape = (5, flag[0]), dtype = numpy.float32)
        for i in range(I[0], I[1]):
            kmerFrequencyModel(temp, i, kmerFreq_, kmerFreqSquaredRowSum_, length_)
            coverageModel(temp[1 : ], i, covMean_, covVar_)
            numpy.putmask(temp[0], temp[0] < minKmerFreqP, log0) # kmerFreq model #
            numpy.putmask(temp[1], temp[1] < minCovP, log0) # cov model #
            temp[ : 2, i] = log0
            if flag[2] != -1: # neighborX1, neighborY1 #
                numpy.dot(modelW[flag[2]], temp[ : 2], out = temp[2])
                neighbors = numpy.argpartition(temp[2], kth = -nNeighbors1)[-nNeighbors1 : ]
                neighborX1[i, 0, : nNeighbors1] = neighbors
                neighborY1[i, 0, : nNeighbors1] = temp[2, neighbors]
            else: # neighborX2, neighborY2 #
                for j in range(nModelW):
                    numpy.dot(modelW[j], temp[ : 2], out = temp[2])
                    neighbors = numpy.argpartition(temp[2], kth = -nNeighbors2)[-nNeighbors2 : ]
                    neighborX2[i, j, : nNeighbors2] = neighbors
                    neighborY2[i, j, : nNeighbors2] = temp[2, neighbors]
        processQ.task_done()
    return None


def createProcesses(
    kmerFreq, kmerFreqSquaredRowSum, length, kmerFreqDim, minKmerFreqP,
    covMean, covVar, covDim, minCovP,
    nModelW, modelW, n1, m1, n2, m2, t
):
    '''
    nModelW: number of model weights
    modelW: model weights
    n1: number of sequences
    m1: number of neighbors
    n2: number of seed sequences
    m2: number of neighbors
    '''

    sharedNeighborX = RawArray(c_int32, max(n1 * m1, n2 * m2 * nModelW))
    sharedNeighborY = RawArray(c_float, max(n1 * m1, n2 * m2 * nModelW))
    neighborX1 = numpy.ndarray(shape = (n1, 1, m1), dtype = numpy.int32, buffer = sharedNeighborX)
    neighborY1 = numpy.ndarray(shape = (n1, 1, m1), dtype = numpy.float32, buffer = sharedNeighborY)
    neighborX2 = numpy.ndarray(shape = (n2, nModelW, m2), dtype = numpy.int32, buffer = sharedNeighborX)
    neighborY2 = numpy.ndarray(shape = (n2, nModelW, m2), dtype = numpy.float32, buffer = sharedNeighborY)
    nodes = RawArray(c_int32, n1)
    flag = RawArray(c_int32, 3) # nodesSize, token, modelFlag #
    processQ = JoinableQueue(t)
    processes = list()
    for i in range(t):
        processes.append(
            Process(
                target = workerProcess,
                args = (
                    processQ, nodes, n1, m1, n2, m2,
                    sharedNeighborX, sharedNeighborY,
                    nModelW, modelW,
                    flag,
                    kmerFreq, kmerFreqSquaredRowSum, length, kmerFreqDim, log(minKmerFreqP) if minKmerFreqP else log0,
                    covMean, covVar, covDim, log(minCovP) if minCovP else log0
                )
            )
        )
        processes[-1].start()
    return (processQ, nodes, flag, neighborX1, neighborY1, neighborX2, neighborY2, processes)


def freeProcesses(processQ, processes):
    for process in processes:
        processQ.put(None)
    processQ.close()
    processQ.join()
    for process in processes:
        process.join()
        process.close()
    return None


def createAffinity(processQ, nodes, affinities, container, flag, sequences, nNeighbors, nModelWeights, batchSize):
    container[ : sequences.size] = sequences
    flag[0] = sequences.size
    flag[1] += 1
    if nModelWeights:
        flag[2] = -1
    else:
        nModelWeights = 1
    for i in range(0, sequences.size, batchSize):
        processQ.put((i, min(i + batchSize, sequences.size)))
    processQ.join()

    for i in range(nModelWeights):
        x = csr_array(
            (
                numpy.exp(affinities[ : sequences.size, i, : nNeighbors].flat),
                nodes[ : sequences.size, i, : nNeighbors].flatten(),
                numpy.arange(sequences.size + 1) * nNeighbors
            ),
            shape = (sequences.size, sequences.size), dtype = numpy.float32, copy = False
        )
        yield x
    return None