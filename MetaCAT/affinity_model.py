from ctypes import c_float, c_int32
from math import log
from multiprocessing import JoinableQueue, Process
from multiprocessing.sharedctypes import RawArray, RawValue

import numpy
from scipy.sparse import csr_array
from threadpoolctl import threadpool_limits

from .coverage_model import model as coverageModel
from .kmer_frequency_model import model as kmerFrequencyModel

log0 = -1e10


def workerProcess(
    processQ, indices, nIndices, token, n1, n2, m,
    neighborIndices, neighborAffinities,
    nModelW, modelW, modelFlag,
    kmerFreq, kmerFreqSquaredRowSum, length, kmerFreqDim, minKmerFreqP,
    covMean, covVar, covDim, minCovP
):
    # shared data #
    indices = numpy.ndarray(shape = n1, dtype = numpy.int32, buffer = indices)
    neighborIndices1 = numpy.ndarray(shape = (n1, 1, m), dtype = numpy.int32, buffer = neighborIndices)
    neighborAffinities1 = numpy.ndarray(shape = (n1, 1, m), dtype = numpy.float32, buffer = neighborAffinities)
    neighborIndices2 = numpy.ndarray(shape = (n2, nModelW, m), dtype = numpy.int32, buffer = neighborIndices)
    neighborAffinities2 = numpy.ndarray(shape = (n2, nModelW, m), dtype = numpy.float32, buffer = neighborAffinities)
    modelW = numpy.ndarray(shape = (nModelW, 2), dtype = numpy.float32, buffer = modelW)
    kmerFreq = numpy.ndarray(shape = (n1, kmerFreqDim), dtype = numpy.float32, buffer = kmerFreq)
    kmerFreqSquaredRowSum = numpy.ndarray(shape = n1, dtype = numpy.float32, buffer = kmerFreqSquaredRowSum)
    length = numpy.ndarray(shape = (6, n1), dtype = numpy.float32, buffer = length)
    covMean = numpy.ndarray(shape = (n1, covDim), dtype = numpy.float32, buffer = covMean)
    covVar = numpy.ndarray(shape = (n1, covDim), dtype = numpy.float32, buffer = covVar)
    token_ = -1

    while True:
        i = processQ.get()
        if i is None:
            processQ.task_done()
            break
        if token_ != token.value:
            token_ = token.value
            indices_ = indices[ : nIndices.value]
            kmerFreq_ = kmerFreq[indices_]
            kmerFreqSquaredRowSum_ = kmerFreqSquaredRowSum[indices_]
            length_ = length[ : , indices_]
            covMean_ = covMean[indices_]
            covVar_ = covVar[indices_]
            nNeighbors = min(m, nIndices.value)
            temp = numpy.empty(shape = (5, nIndices.value), dtype = numpy.float32)

        kmerFrequencyModel(temp, i, kmerFreq_, kmerFreqSquaredRowSum_, length_)
        coverageModel(temp[1 : ], i, covMean_, covVar_)

        temp[0][temp[0] < minKmerFreqP] = log0 # kmerFreq model #
        temp[1][temp[1] < minCovP] = log0 # cov model #
        temp[ : 2, i] = log0

        if modelFlag.value != -1: # neighborIndices1, neighborAffinities1 #
            numpy.dot(modelW[modelFlag.value], temp[ : 2], out = temp[2])
            neighbors = numpy.argpartition(temp[2], kth = -nNeighbors)[-nNeighbors : ]
            neighbors = neighbors[numpy.argsort(temp[2, neighbors])[ : : -1]] # sort from high to low #
            neighborIndices1[i, 0, : nNeighbors] = neighbors
            neighborAffinities1[i, 0, : nNeighbors] = temp[2, neighbors]
        else: # neighborIndices2, neighborAffinities2 #
            for j in range(nModelW):
                numpy.dot(modelW[j], temp[ : 2], out = temp[2])
                neighbors = numpy.argpartition(temp[2], kth = -nNeighbors)[-nNeighbors : ]
                neighbors = neighbors[numpy.argsort(temp[2, neighbors])[ : : -1]] # sort from high to low #
                neighborIndices2[i, j, : nNeighbors] = neighbors
                neighborAffinities2[i, j, : nNeighbors] = temp[2, neighbors]
        processQ.task_done()
    return None


def createProcesses(
    kmerFreq, kmerFreqSquaredRowSum, length, kmerFreqDim, minKmerFreqP,
    covMean, covVar, covDim, minCovP,
    nModelW, modelW, n1, n2, m, t
):
    '''
    nModelW: number of model weights
    modelW: model weights
    n1: number of sequences
    n2: number of seed sequences
    m: number of neighbors
    '''

    sharedNeighborIndices = RawArray(c_int32, m * max(n1, n2 * nModelW))
    sharedNeighborAffinities = RawArray(c_float, m * max(n1, n2 * nModelW))
    neighborIndices1 = numpy.ndarray(shape = (n1, 1, m), dtype = numpy.int32, buffer = sharedNeighborIndices)
    neighborAffinities1 = numpy.ndarray(shape = (n1, 1, m), dtype = numpy.float32, buffer = sharedNeighborAffinities)
    neighborIndices2 = numpy.ndarray(shape = (n2, nModelW, m), dtype = numpy.int32, buffer = sharedNeighborIndices)
    neighborAffinities2 = numpy.ndarray(shape = (n2, nModelW, m), dtype = numpy.float32, buffer = sharedNeighborAffinities)
    modelFlag = RawValue(c_int32, -1)
    container = RawArray(c_int32, n1)
    containerSize = RawValue(c_int32, 0)
    token = RawValue(c_int32, 0)
    processQ = JoinableQueue(t)
    processes = list()
    with threadpool_limits(limits = 1):
        for i in range(t):
            processes.append(
                Process(
                    target = workerProcess,
                    args = (
                        processQ, container, containerSize, token, n1, n2, m,
                        sharedNeighborIndices, sharedNeighborAffinities,
                        nModelW, modelW, modelFlag,
                        kmerFreq, kmerFreqSquaredRowSum, length, kmerFreqDim, log(minKmerFreqP) if minKmerFreqP else log0,
                        covMean, covVar, covDim, log(minCovP) if minCovP else log0
                    )
                )
            )
            processes[-1].start()
    return (processQ, container, containerSize, token, modelFlag, neighborIndices1, neighborAffinities1, neighborIndices2, neighborAffinities2, processes)


def freeProcesses(processQ, processes):
    for process in processes:
        processQ.put(None)
    processQ.close()
    processQ.join()
    for process in processes:
        process.join()
        process.close()
    return None


def createAffinity(processQ, neighborIndices, neighborAffinities, containerSize, container, nNeighbors, x, flag):
    containerSize.value = x.size
    container[ : x.size] = x
    nNeighbors = min(nNeighbors, x.size)
    # put data into queue #
    for i in range(x.size):
        processQ.put(i)
    processQ.join()
    for i in range(flag):
        y = csr_array(
            (
                numpy.exp(neighborAffinities[ : x.size, i, : nNeighbors].flat),
                neighborIndices[ : x.size, i, : nNeighbors].flatten(),
                numpy.arange(x.size + 1) * nNeighbors
            ),
            shape = (x.size, x.size),
            dtype = numpy.float32,
            copy = False
        )
        yield y
    return None
