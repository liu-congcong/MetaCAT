from ctypes import c_float, c_int32
from math import ceil
from multiprocessing import JoinableQueue, Process
from multiprocessing.sharedctypes import RawArray
from operator import itemgetter

import numpy
from threadpoolctl import threadpool_limits

from .marker import markerHashPartition, nArcMarkers, nBacMarkers, nMarkers

highQualityScore = 0.5
minScore = 0.1


def get_score(x, flag, w1 = 0.5, w2 = 0.5):
    '''
    flag: 0 sort by score, 1 sort by recall
    '''
    arcMarkers = x[markerHashPartition[0] : markerHashPartition[1]]
    uArcMarkers = numpy.count_nonzero(arcMarkers) # unique arc markers #
    rArcMarkers = numpy.sum(arcMarkers) - uArcMarkers # redundant arc markers #
    arcRecall = uArcMarkers / nArcMarkers
    arcScore = arcRecall - w1 * numpy.sum(arcMarkers > 1) / (uArcMarkers + 1e-10) - w2 * rArcMarkers / nArcMarkers

    bacMarkers = x[markerHashPartition[2] : markerHashPartition[3]]
    uBacMarkers = numpy.count_nonzero(bacMarkers) # unique bac markers #
    rBacMarkers = numpy.sum(bacMarkers) - uBacMarkers # redundant bac markers #
    bacRecall = uBacMarkers / nBacMarkers
    bacScore = bacRecall - w1 * numpy.sum(bacMarkers > 1) / (uBacMarkers + 1e-10) - w2 * rBacMarkers / nBacMarkers
    return max((arcScore, arcRecall), (bacScore, bacRecall), key = itemgetter(flag))


def get_score_(x, flag):
    arcMarkers = x[markerHashPartition[0] : markerHashPartition[1]]
    uArcMarkers = numpy.count_nonzero(arcMarkers)
    arcRecall = uArcMarkers / nArcMarkers
    arcPrecision = uArcMarkers / (numpy.sum(arcMarkers) + 1e-10)
    arcF1 = 2.0 * arcRecall * arcPrecision / (arcRecall + arcPrecision + 1e-10)

    bacMarkers = x[markerHashPartition[2] : markerHashPartition[3]]
    uBacMarkers = numpy.count_nonzero(bacMarkers)
    bacRecall = uBacMarkers / nBacMarkers
    bacPrecision = uBacMarkers / (numpy.sum(bacMarkers) + 1e-10)
    bacF1 = 2.0 * bacRecall * bacPrecision / (bacRecall + bacPrecision + 1e-10)
    return max((arcF1, arcRecall), (bacF1, bacRecall), key = itemgetter(flag))


def evaluate_sequences(sequences, seed2markers, nSeeds):
    '''
    Parameters:
        sequences: the list of sequences, need to be sorted by length in reverse order.
        seed2markers: the hash of seed - markers pairs, sorted markers.
        nSeeds: max number of seeds.
    Return:
        score, recall, k, seeds, mappingSets
    '''
    x = numpy.zeros(shape = nMarkers, dtype = numpy.int32)
    seeds = list()
    marker2mappings = dict()
    i = 0
    for sequence in sequences:
        if sequence in seed2markers:
            x[seed2markers[sequence]] += 1
            for marker in seed2markers[sequence]:
                marker2mappings.setdefault(marker, list()).append(i)
            seeds.append(sequence)
            i += 1
            if i > nSeeds:
                break
    mappingSets = set()
    k2weight = dict()
    if marker2mappings:
        for marker, mappings in marker2mappings.items():
            # mappings must be a sorted tuple, not an array #
            k = len(mappings)
            mappingSets.add(tuple(mappings))
            k2weight[k] = k2weight.get(k, 0) + 1
        maxWeight = max(k2weight.values())
        k = max(k for k, weight in k2weight.items() if weight >= maxWeight * 0.5)
        if k == 1 and sum(k2weight.values()) - k2weight[1] > 5:
            k = 2
    else:
        k = 0
    return (*get_score(x, 1), k, numpy.asarray(seeds, dtype = numpy.int64), [numpy.asarray(i, dtype = numpy.int64) for i in sorted(mappingSets)])


def workerProcess(processQ, sequences, markerMatrix, x, y1, y2, xn, yn, m):
    '''
    processQ:
    sequences: total seeds
    markerMatrix: (nTotalSeeds, nMarkers)
    x: lpPartitions, (nModelWeights * nMarkers, nTotalSeeds)
    y1: (nModelWeights, nTotalSeeds)
    y2: (nModelWeights, )
    xn: nModelWeights * nMarkers
    yn: nModelWeights
    m: nTotalSeeds

    fill y1 and y2
    '''
    sequences = numpy.ndarray(shape = m, dtype = numpy.int32, buffer = sequences)
    markerMatrix = numpy.ndarray(shape = (m, nMarkers), dtype = numpy.int32, buffer = markerMatrix)
    x = numpy.ndarray(shape = (xn, m), dtype = numpy.int32, buffer = x)
    y1 = numpy.ndarray(shape = (yn, m), dtype = numpy.int32, buffer = y1)
    while True:
        sequences_, xi, xj, yi = processQ.get()
        if sequences_ is None:
            processQ.task_done()
            break
        markerMatrix_ = markerMatrix[numpy.flatnonzero(numpy.isin(sequences, sequences_))] # (nSequences_, mMarkers) #
        partitions = x[xi : xj, : sequences_.size] + 1 # partitions += 1, partitions >= 0, unclustered sequences = 0 #
        y1[yi, : sequences_.size] = -1
        scoreIJ = dict()
        for i in range(markerMatrix_.shape[0]):
            scoreIJ[(-1, i)] = get_score(markerMatrix_[i], 0)[0]
        for i, partition in enumerate(partitions):
            for j in numpy.unique(partition):
                if j:
                    scoreIJ[(i, j)] = get_score(numpy.sum(markerMatrix_[partition == j], axis = 0), 0)[0]
        y2[yi] = 0
        while scoreIJ:
            (i, j), score = max(scoreIJ.items(), key = itemgetter(1))
            if i == -1:
                mask = numpy.array([j], dtype = numpy.int64)
            else:
                mask = numpy.flatnonzero(partitions[i] == j)
            for k in mask:
                del scoreIJ[(-1, k)]
            if score >= minScore:
                y1[yi, mask] = mask[0]
            if score >= highQualityScore:
                y2[yi] += score
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


def workerProcessLow(processQ, sequences, markerMatrix, x, y1, y2, xn, yn, m):
    '''
    processQ:
    sequences: total seeds
    markerMatrix: (nTotalSeeds, nMarkers)
    x: lpPartitions, (nModelWeights * nMarkers, nTotalSeeds)
    y1: (nModelWeights, nTotalSeeds)
    y2: (nModelWeights, )
    xn: nModelWeights * nMarkers
    yn: nModelWeights
    m: nTotalSeeds

    fill y1 and y2
    '''
    sequences = numpy.ndarray(shape = m, dtype = numpy.int32, buffer = sequences)
    markerMatrix = numpy.ndarray(shape = (m, nMarkers), dtype = numpy.int32, buffer = markerMatrix)
    x = numpy.ndarray(shape = (xn, m), dtype = numpy.int32, buffer = x)
    y1 = numpy.ndarray(shape = (yn, m), dtype = numpy.int32, buffer = y1)
    scoresDataType = numpy.dtype([('flag', numpy.bool_), ('marker', numpy.int32), ('label', numpy.int64), ('score', numpy.float32)])
    temp = list()
    while True:
        sequences_, xi, xj, yi = processQ.get()
        if sequences_ is None:
            processQ.task_done()
            break
        markerMatrix_ = markerMatrix[numpy.flatnonzero(numpy.isin(sequences, sequences_))] # (nSequences_, mMarkers) #
        X = x[xi : xj, : sequences_.size] + 1 # partitions += 1 #
        n = 0
        for xi in X:
            uniqueXi = numpy.unique(xi)
            uniqueXiOffset = 0 if uniqueXi[0] else 1
            temp.append(uniqueXi[uniqueXiOffset : ])
            n += uniqueXi.size - uniqueXiOffset
        scores = numpy.empty(shape = n + sequences_.size, dtype = scoresDataType)
        for i in range(sequences_.size):
            scores[i] = (True, -1, i, get_score(markerMatrix_[i], 0)[0])
        n = sequences_.size
        for i, (xi, uniqueXi) in enumerate(zip(X, temp)):
            for xj in uniqueXi:
                scores[n] = (True, i, xj, get_score(numpy.sum(markerMatrix_[xi == xj], axis = 0), 0)[0])
                n += 1
        temp.clear()
        y1[yi, : sequences_.size] = -1
        y1label = 0
        while numpy.any(scores['flag']):
            mask = numpy.flatnonzero(scores['flag'] == True)
            i = mask[numpy.argmax(scores['score'][mask])] # scores[i]: max score #
            scoreI = scores[i] # flag, marker, label, score #
            mask = numpy.flatnonzero(X[scoreI[1]] == scoreI[2])
            scoreI[0] = False
            if scoreI[3] >= minScore:
                y1[yi, mask] = y1label
                y1label += 1
            if scoreI[3] >= highQualityScore:
                y2[yi] += scoreI[3]
            X[ : , mask] = -1
            for i, xi in enumerate(X):
                if i != scoreI[1]:
                    uniqueXi = numpy.unique(xi[mask])
                    for xj in uniqueXi:
                        j = numpy.flatnonzero((scores['marker'] == i) & (scores['label'] == xj))
                        scores['score'][j] = get_score(numpy.sum(markerMatrix_[xi == xj], axis = 0), 0)[0]
        processQ.task_done()
    return None


def createProcesses(sequences, markerTable, x, xn, yn, t):
    '''
    sequences: sorted total seeds, [int, ]
    markerTable: (nTotalSeeds, nMarkers)
    x: lpPartitions, (nModelWeights * nMarkers, nSeeds)
    xn: lpPartitionN, nModelWeights * nMarkers, int
    yn: nModelWeights, int
    t: threads, int
    '''
    m = len(sequences)

    # optimal partitions #
    sharedOptimalPartitions = RawArray(c_int32, yn * m)
    optimalPartitions = numpy.ndarray(
        shape = (yn, m),
        dtype = numpy.int32,
        buffer = sharedOptimalPartitions
    )

    # optimal partition scores #
    sharedScores = RawArray(c_float, yn)
    scores = numpy.ndarray(shape = yn, dtype = numpy.float32, buffer = sharedScores)

    processQ = JoinableQueue(yn)
    processes = list()
    with threadpool_limits(limits = ceil(t / yn)):
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
    '''
    for i in range(nBlocks):
        processQ.put((sequences, i * blockSize, (i + 1) * blockSize, i))
    processQ.join()
    return None


def createMarkerTable(sequence2markers):
    sequences = sorted(sequence2markers.keys())
    n = len(sequences)
    x = RawArray(c_int32, n * nMarkers)
    y = numpy.ndarray(shape = (n, nMarkers), dtype = numpy.int32, buffer = x)
    for i, sequence in enumerate(sequences):
        y[i, sequence2markers[sequence]] = 1
    return x
