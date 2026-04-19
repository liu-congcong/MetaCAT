from ctypes import c_float, c_int32, c_int8
from multiprocessing import JoinableQueue, Process
from multiprocessing.sharedctypes import RawArray
from operator import itemgetter

import numpy
from threadpoolctl import threadpool_limits

from .marker import markerHashPartition, nArcMarkers, nBacMarkers, nMarkers


def get_score(x, flag, w1 = 0.5, w2 = 0.5):
    '''
    x: (nMarkers, )
    flag: 0 sort by score, 1 sort by recall
    '''
    arcMarkers = x[markerHashPartition[0] : markerHashPartition[1]]
    uArcMarkers = numpy.count_nonzero(arcMarkers) # unique arc markers #
    rArcMarkers = numpy.sum(arcMarkers) - uArcMarkers # redundant arc markers #
    arcRecall = uArcMarkers / nArcMarkers
    arcScore = arcRecall - w1 * numpy.sum(arcMarkers > 1) / (uArcMarkers + 10 * numpy.finfo(numpy.float32).eps) - w2 * rArcMarkers / nArcMarkers

    bacMarkers = x[markerHashPartition[2] : markerHashPartition[3]]
    uBacMarkers = numpy.count_nonzero(bacMarkers) # unique bac markers #
    rBacMarkers = numpy.sum(bacMarkers) - uBacMarkers # redundant bac markers #
    bacRecall = uBacMarkers / nBacMarkers
    bacScore = bacRecall - w1 * numpy.sum(bacMarkers > 1) / (uBacMarkers + 10 * numpy.finfo(numpy.float32).eps) - w2 * rBacMarkers / nBacMarkers
    return max((arcScore, arcRecall), (bacScore, bacRecall), key = itemgetter(flag))


def evaluate_cluster(markerTable):
    '''
    Return:
        score, recall, k, mappingSets
    '''
    score, recall = get_score(numpy.sum(markerTable, axis = 0), 1)
    k2weight = dict()
    mappings = set()
    for mapping in markerTable.T:
        mapping = numpy.flatnonzero(mapping)
        if mapping.size:
            k2weight[mapping.size] = k2weight.get(mapping.size, 0) + 1
            mappings.add(tuple(mapping))
    k2weight[0] = 0
    maxWeight = max(k2weight.values())
    k = max(k for k, weight in k2weight.items() if weight >= maxWeight * 0.5)
    if k == 1 and sum(k2weight.values()) - k2weight[1] > 5:
        k = 2
    return (score, recall, k, [numpy.asarray(i, dtype = numpy.int32) for i in sorted(mappings)])


@threadpool_limits.wrap(limits = 1)
def workerProcess(processQ, sequences, markerMatrix, x, y1, y2, xn, yn, m):
    '''
    processQ:
    sequences: total seeds
    markerMatrix: (nTotalSeeds, nMarkers), int8
    x: lpPartitions, (nModelWeights * nMarkers, nTotalSeeds)
    y1: (nModelWeights, nTotalSeeds), >0 for high quality, =0 for unprocessed, <0 for low quality
    y2: (nModelWeights, 2)
    xn: nModelWeights * nMarkers
    yn: nModelWeights
    m: nTotalSeeds
    '''
    sequences = numpy.ndarray(shape = m, dtype = numpy.int32, buffer = sequences)
    markerMatrix = numpy.ndarray(shape = (m, nMarkers), dtype = numpy.int8, buffer = markerMatrix) # (nTotalSeeds, nMarkers) #
    x = numpy.ndarray(shape = (xn, m), dtype = numpy.int32, buffer = x) # (nModelWeights * nMarkers, nTotalSeeds) #
    y1 = numpy.ndarray(shape = (yn, m), dtype = numpy.int32, buffer = y1) # (nModelWeights, nTotalSeeds) #
    y2 = numpy.ndarray(shape = (yn, 2), dtype = numpy.float32, buffer = y2) # (nModelWeights, 2) #
    while True:
        sequences_, xi, xj, yi = processQ.get()
        '''
        sequences_: array of seeds
        xi: partitions start
        xj: partitions end
        xi - xj: groups of seeds
        yi: modelWeight index
        '''
        if sequences_ is None:
            processQ.task_done()
            break
        markerMatrix_ = markerMatrix[numpy.flatnonzero(numpy.isin(sequences, sequences_))] # (nSequences_, mMarkers) #
        partitions = x[xi : xj, : sequences_.size] + 1 # partitions += 1, partitions >= 0, unclustered sequences = 0 #
        # fill scoreIJ #
        scoreIJ = dict() # i: group index, j: j-th seed in seeds #
        for i in range(markerMatrix_.shape[0]):
            scoreIJ[(-1, i)] = get_score(markerMatrix_[i], 0)[0]
        for i, partition in enumerate(partitions):
            for j in numpy.unique(partition):
                if j:
                    scoreIJ[(i, j)] = get_score(numpy.sum(markerMatrix_[partition == j], axis = 0), 0)[0]
        y1[yi, : sequences_.size] = 0
        y1i = 1
        y2[yi] = 0
        while scoreIJ:
            (i, j), score = max(scoreIJ.items(), key = itemgetter(1))
            if i == -1:
                mask = numpy.array([j], dtype = numpy.int32)
            else:
                mask = numpy.flatnonzero(partitions[i] == j)
            for k in mask:
                del scoreIJ[(-1, k)]
            if score >= 0.5:
                y1[yi, mask] = y1i
                y1i += 1
                y2[yi, 0] += score
            elif score >= 0.1:
                y1[yi, mask] = -y1i
                y1i += 1
                y2[yi, 1] += score
            partitions[ : , mask] *= -1
            for i, partition in enumerate(partitions):
                for j in numpy.unique(partition[mask]):
                    if j:
                        mask_ = partition == -j
                        if numpy.any(mask_):
                            scoreIJ[(i, -j)] = get_score(numpy.sum(markerMatrix_[mask_], axis = 0), 0)[0]
                        else:
                            del scoreIJ[(i, -j)]
            partitions[ : , mask] = 0
        processQ.task_done()
    return None


def createProcesses(sequences, markerTable, x, xn, yn):
    '''
    sequences: sorted total seeds, [int, ]
    markerTable: (nTotalSeeds, nMarkers), int8
    x: lpPartitions, (nModelWeights * nMarkers, nSeeds)
    xn: lpPartitionN, nModelWeights * nMarkers, int
    yn: nModelWeights, int
    '''
    m = len(sequences)
    sharedOptimalPartitions = RawArray(c_int32, yn * m)
    optimalPartitions = numpy.ndarray(
        shape = (yn, m),
        dtype = numpy.int32,
        buffer = sharedOptimalPartitions
    )
    sharedScores = RawArray(c_float, yn * 2)
    scores = numpy.ndarray(shape = (yn, 2), dtype = numpy.float32, buffer = sharedScores)
    processQ = JoinableQueue(yn)
    processes = list()
    for i in range(yn):
        processes.append(
            Process(
                target = workerProcess,
                args = (
                    processQ, sequences, markerTable, x,
                    sharedOptimalPartitions, sharedScores, xn, yn, m
                )
            )
        )
        processes[-1].start()
    return (processQ, optimalPartitions, scores, processes)


def freeProcesses(processQ, processes):
    for process in processes:
        processQ.put((None, None, None, None))
    processQ.close()
    processQ.join()
    for process in processes:
        process.join()
        process.close()
    return None


def score_partitions(processQ, sequences, blockSize, nBlocks):
    '''
    nBlocks: nModelWeights
    sequences: array of seeds
    '''
    for i in range(nBlocks):
        processQ.put((sequences, i * blockSize, (i + 1) * blockSize, i))
    processQ.join()
    return None


def createSeedMarker(sequence2markers):
    sequences = sorted(sequence2markers.keys())
    n = len(sequences)
    x = RawArray(c_int32, n)
    X = numpy.ndarray(shape = n, dtype = numpy.int32, buffer = x)
    y = RawArray(c_int8, n * nMarkers)
    Y = numpy.ndarray(shape = (n, nMarkers), dtype = numpy.int8, buffer = y)
    for i, sequence in enumerate(sequences):
        x[i] = sequence
        Y[i, sequence2markers[sequence]] = 1
    return (x, X, y, Y)
