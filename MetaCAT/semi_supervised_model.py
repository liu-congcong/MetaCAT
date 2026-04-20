from math import ceil
import numpy

from .affinity_model import createAffinity
from .label_propagation import annealingLabelPropagation
from .score_model import score_partitions


def seedModel(
    affinityProcessQ, affinityContainer, affinityFlag, affinityNeighborIndices, affinityNeighborAffinities,
    lpProcessQ,
    scoreProcessQ, scorePartitions, scoreScores,
    seeds, mappingSets, nNeighbors, nModelWeights, nThreads
):
    '''
    affinityFlag: containerSize, token, modelFlag
    scorePartitions: (nModelWeights, nSeeds), >0 for high quality, =0 for unprocessed, <0 for low quality
    scoreScores: (nModelWeights, 2)
    seeds: array
    mappingSets: list of array
    '''
    nNeighbors = min(seeds.size, nNeighbors)
    nPartitions = 0
    probabilities = list()
    for x in createAffinity(affinityProcessQ, affinityNeighborIndices, affinityNeighborAffinities, affinityContainer, affinityFlag, seeds, nNeighbors, nModelWeights, ceil(0.2 * seeds.size / nThreads)):
        probabilities.append(numpy.percentile(numpy.min(x.data.reshape(seeds.size, nNeighbors), axis = 1), 20))
        for mappingSet in mappingSets:
            if mappingSet.size > 1 and mappingSet.size < seeds.size:
                lpProcessQ.put((nPartitions, x, mappingSet))
                nPartitions += 1
    lpProcessQ.join()

    # score partitions, fill scorePartitions and scoreScores #
    score_partitions(scoreProcessQ, seeds, int(nPartitions / nModelWeights), nModelWeights)
    i = numpy.argmax(scoreScores, axis = 0)
    if scoreScores[i[0], 0] > 0:
        affinityFlag[2] = i[0]
    else:
        affinityFlag[2] = i[1]
    return (scorePartitions[affinityFlag[2], : seeds.size], max(0.8, probabilities[affinityFlag[2]]))


def sequenceModel(
    affinityProcessQ, affinityContainer, affinityFlag, affinityNeighborIndices, affinityNeighborAffinities,
    sequences, seedSequences, seedLabels, skip0, nNeighbors, probability, nThreads
):
    '''
    affinityFlag: containerSize, token, modelFlag
    seedLabels: >0 for high quality, =0 for unprocessed, <0 for low quality
    skip0: skip 0 if True
    '''
    nNeighbors = min(sequences.size, nNeighbors)
    mappings = numpy.flatnonzero(numpy.isin(sequences, seedSequences, assume_unique = True))
    y = numpy.full_like(sequences, -1, dtype = numpy.int32)
    for i in numpy.unique(seedLabels):
        if i or (not skip0):
            mask = seedLabels == i
            y[mappings[mask]] = numpy.min(seedSequences[mask])
    uniqueY = numpy.unique(y)

    if uniqueY[0] == -1 and uniqueY.size > 2: # unlabeled and labeled data #
        # affinity matrix #
        x = next(
            createAffinity(
                affinityProcessQ, affinityNeighborIndices, affinityNeighborAffinities, affinityContainer, affinityFlag,
                sequences, nNeighbors, 0, ceil(0.2 * sequences.size / nThreads)
            )
        )
        x.data[x.data < probability] = 0
        y[ : ] = annealingLabelPropagation(x, y, ((0.97, 10, 0.97), (0.90, 10, 0.90), (0.85, 30, 0.7)))
        mask = y == -1
        if numpy.any(mask):
            y[mask] = numpy.min(sequences[mask])
    elif uniqueY[0] == -1 or uniqueY.size == 1:
        y[ : ] = -1 - numpy.min(sequences)
    else:
        for i in uniqueY:
            mask = y == i
            y[mask] = numpy.min(sequences[mask])
    return y
