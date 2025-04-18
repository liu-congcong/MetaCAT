import numpy

from .affinity_model import createAffinity
from .label_propagation import annealingLabelPropagation, calculatePartition
from .score_model import score_partitions


def seedModel(
    affinityProcessQ, affinityContainer, affinityContainerSize, affinityToken, affinityNeighborIndices, affinityNeighborAffinities, affinityModelFlag,
    lpProcessQ,
    scoreProcessQ, scorePartitions, scoreScores,
    seeds, mappingSets, nNeighbors, nModelWeights
):
    '''
    scorePartitions: (nModelWeights, nSeeds)
    scoreScores: (nModelWeights, )
    seeds: array
    mappingSets: list of array
    '''
    # affinity matrix #
    affinityContainerSize.value = seeds.size
    affinityContainer[ : affinityContainerSize.value] = seeds
    affinityModelFlag.value = -1
    x = createAffinity(
        affinityProcessQ, affinityNeighborIndices, affinityNeighborAffinities, affinityContainerSize,
        affinityContainer, nNeighbors, seeds, nModelWeights
    ) # generator #
    nPartitions = calculatePartition(lpProcessQ, x, mappingSets) # nModelWeights * #mappingSets #
    affinityToken.value += 1
    # score partitions, fill scorePartitions and scoreScores #
    score_partitions(scoreProcessQ, seeds, int(nPartitions / nModelWeights), nModelWeights)
    i = numpy.argmax(scoreScores)
    affinityModelFlag.value = i
    return (scorePartitions[i, : seeds.size])


def sequenceModel(
    affinityProcessQ, affinityContainer, affinityContainerSize, affinityToken, affinityNeighborIndices, affinityNeighborAffinities,
    sequences, seedSequences, seedLabels, nNeighbors
):
    mappings = numpy.flatnonzero(numpy.isin(sequences, seedSequences, assume_unique = True))
    y = numpy.full_like(sequences, -1, dtype = numpy.int32)
    for i in numpy.unique(seedLabels):
        mask = seedLabels == i
        y[mappings[mask]] = numpy.min(seedSequences[mask])
    uniqueY = numpy.unique(y)
    if uniqueY[0] == -1 and uniqueY.size > 2: # unlabeled and labeled data #
        # affinity matrix #
        affinityContainerSize.value = sequences.size
        affinityContainer[ : affinityContainerSize.value] = sequences
        x = createAffinity(
            affinityProcessQ, affinityNeighborIndices, affinityNeighborAffinities, affinityContainerSize,
            affinityContainer, nNeighbors, sequences, 1
        ) # x: generator #

        y = annealingLabelPropagation(next(x), y, ((0.97, 10, 0.97), (0.9, 10, 0.9), (0.7, 50, 0.7)))
        affinityToken.value += 1

        mask = y == -1
        if numpy.any(mask):
            y[mask] = numpy.min(sequences[mask])

    else:
        y = numpy.full_like(sequences, -1 - numpy.min(sequences), dtype = numpy.int32)
    return y
