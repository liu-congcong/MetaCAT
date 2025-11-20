import gzip
import os
import sys
from ctypes import c_float
from datetime import datetime
from math import ceil, log, log10
from multiprocessing.sharedctypes import RawArray
from operator import itemgetter

import numpy
from sklearn.decomposition import PCA

from .affinity_model import createProcesses as createAffinityProcesses
from .affinity_model import freeProcesses as freeAffinityProcesses
from .fasta import readFastaFile
from .kmer import count_kmers
from .kmer_frequency_model import parseW as parseSequenceLength
from .label_propagation import createProcesses as createLPProcesses
from .label_propagation import freeProcesses as freeLPProcesses
from .marker import markerHash, nMarkers
from .processbar import ProcessBar
from .score_model import createProcesses as createScoreProcesses
from .score_model import createSeedMarker, evaluate_cluster
from .score_model import freeProcesses as freeScoreProcess
from .semi_supervised_model import seedModel as clusterSeeds
from .semi_supervised_model import sequenceModel as clusterSequences


def is_gzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def read_seed_file(file):
    sequence2markers = dict()
    if is_gzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, mode = 'rt')
    openFile.readline()
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        sequence2markers[lines[0]] = [int(i) for i in lines[1].split(' ')]
    openFile.close()
    return sequence2markers


def determine_samples(file):
    if is_gzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, mode = 'rt')
    n = int(openFile.readline().count('\t') * 0.5)
    openFile.close()
    return n


def read_coverage_file(file, seqID2index):
    if is_gzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, mode = 'rt')
    n = len(seqID2index)
    m = int(openFile.readline().count('\t') * 0.5)
    sharedCovMean = RawArray(c_float, n * m) # float32 #
    covMean = numpy.ndarray(shape = (n, m), dtype = numpy.float32, buffer = sharedCovMean)
    sharedCovVar = RawArray(c_float, n * m) # float32 #
    covVar = numpy.ndarray(shape = (n, m), dtype = numpy.float32, buffer = sharedCovVar)
    covMeanIndices = tuple(range(1, 2 * m + 1, 2))
    covVarIndices = tuple(range(2, 2 * m + 1, 2))
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[0] in seqID2index:
            i = seqID2index[lines[0]]
            covMean[i] = [lines[j] for j in covMeanIndices]
            covVar[i] = [lines[j] for j in covVarIndices]
    openFile.close()
    covVar /= min(m, 10) # A larger value of m (i.e., a smaller value of variance) enhances difference in the means of the two sequences. #
    covVar += 1e-3
    return (covMean, sharedCovMean, sharedCovVar)


def parseModelWeight(w1, w2, n):
    w1.sort(reverse = True)
    w2.sort()
    if n > 1:
        w1.insert(0, 1 / (log(n + 1) + 1))
        w2.insert(0, 1 / (1 / log(n + 1) + 1))
    w = RawArray(c_float, len(w1) + len(w2))
    for i, (W1, W2) in enumerate(zip(w1, w2)):
        w[2 * i] = W1
        w[2 * i + 1] = W2
    return (i + 1, w)


def output_clusters(sequenceStructs, labels, filePrefix, noFastas, lineLength = 100):
    uniqueLabels = numpy.unique(labels)
    if uniqueLabels[0] == - labels.size - 1:
        uniqueLabels = uniqueLabels[1 : ]
    labelWidth = ceil(log10(uniqueLabels.size + 1))
    clusterPrefix = os.path.basename(filePrefix)
    openMappingFile = open(filePrefix + '.mapping', 'w')
    openMappingFile.write('Sequence ID\tCluster ID\n')
    for labelID, label in enumerate(uniqueLabels, start = 1):
        sequences = numpy.flatnonzero(labels == label)
        if noFastas:
            for sequence in sequences:
                openMappingFile.write(f'{sequenceStructs[sequence][0]}\t{clusterPrefix}.{labelID:0{labelWidth}d}\n')
        else:
            openFastaFile = open(f'{filePrefix}.{labelID:0{labelWidth}d}.fasta', 'w')
            for sequence in sequences:
                openMappingFile.write(f'{sequenceStructs[sequence][0]}\t{clusterPrefix}.{labelID:0{labelWidth}d}\n')
                openFastaFile.write('>' + sequenceStructs[sequence][0] + '\n')
                SEQUENCE = sequenceStructs[sequence][1]
                filePointer = 0
                while openFastaFile.write(SEQUENCE[filePointer : filePointer + lineLength]):
                    openFastaFile.write('\n')
                    filePointer += lineLength
            openFastaFile.close()
    openMappingFile.close()
    return None


def readGroundTruthFile(file, sequenceID2sequence): # only for debug #
    i = 0
    genome2i = dict()
    genome = numpy.empty(shape = len(sequenceID2sequence), dtype = numpy.int32)
    if is_gzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, mode = 'rt')
    openFile.readline() # Sequence ID, Genome ID, Length #
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[0] in sequenceID2sequence:
            if lines[1] not in genome2i:
                genome2i[lines[1]] = i
                i += 1
            genome[sequenceID2sequence[lines[0]]] = genome2i[lines[1]]
    openFile.close()
    return genome


def main(parameters):
    if parameters.swdpgmm_engine == 'auto':
        try:
            from .swdpgmm_gpu import SWDPGMM
        except Exception:
            from .swdpgmm_cpu import SWDPGMM
    elif parameters.swdpgmm_engine == 'gpu':
        from .swdpgmm_gpu import SWDPGMM
    else:
        from .swdpgmm_cpu import SWDPGMM

    # determine number of samples #
    nCov = determine_samples(parameters.coverage)
    if parameters.min_sequence_length is None:
        parameters.min_sequence_length = 1000 if nCov == 1 else 500

    # load seed file #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading seeds file.', flush = True)
    sequenceID2markers = read_seed_file(parameters.seed) # sorted markers #
    nSeeds = len(sequenceID2markers)
    if nSeeds == 0:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Empty seeds file.', flush = True)
        sys.exit()

    # read fasta file #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading fasta file.', flush = True)
    sequenceStructs = list()
    for sequenceID, SEQUENCE in readFastaFile(parameters.fasta):
        sequenceLength = len(SEQUENCE)
        if (sequenceLength >= parameters.min_sequence_length) or (sequenceID in sequenceID2markers):
            sequenceStructs.append((sequenceID, SEQUENCE, sequenceLength))
    sequenceStructs.sort(key = itemgetter(2), reverse = True)
    nSequences = len(sequenceStructs)

    sequenceID2sequence = dict()
    length = numpy.empty(shape = nSequences, dtype = numpy.float32)
    for i, (sequenceID, SEQUENCE, sequenceLength) in enumerate(sequenceStructs):
        sequenceID2sequence[sequenceID] = i
        length[i] = sequenceLength

    sequence2markers = dict()
    for sequenceID, markers in sequenceID2markers.items():
        sequence2markers[sequenceID2sequence[sequenceID]] = markers
    del sequenceID2markers

    # read coverage file #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading coverage file.', flush = True)
    covMean, sharedCovMean, sharedCovVar = read_coverage_file(parameters.coverage, sequenceID2sequence)
    covMean = numpy.log1p(covMean)
    covMean -= numpy.min(covMean, axis = 0)
    covMean /= numpy.max(covMean, axis = 0) + 10 * numpy.finfo(numpy.float32).eps # covMean >= 0 and covMean <= 1 #

    # calculate kmer frequency #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Counting kmers of all sequences.', flush = True)
    kmerFreq, sharedKmerFreq = count_kmers(nSequences, sequenceStructs, parameters.threads)
    sharedKmerFreqSquaredRowSum = RawArray(c_float, numpy.einsum('ij,ij->i', kmerFreq, kmerFreq))

    pca = PCA(n_components = None, whiten = False, svd_solver = 'covariance_eigh', random_state = parameters.random_number)
    pca.fit(kmerFreq)
    k = min(
        max(
            numpy.sum(numpy.cumsum(pca.explained_variance_ratio_) < parameters.min_pca_variance_ratio) + 1,
            parameters.min_pca_components
        ),
        numpy.count_nonzero(pca.explained_variance_ratio_ >= 1e-5) # max value of k #
    )
    kmerFreq = pca.transform(kmerFreq)[ : , : k]

    # concatenate coverage and kmer frequency #
    x = numpy.concatenate((kmerFreq, covMean), axis = 1)
    y = numpy.empty(shape = nSequences, dtype = numpy.int32)
    del covMean
    del kmerFreq
    nModelWeights, modelWeights = parseModelWeight(parameters.kmer_frequence_weight, parameters.coverage_weight, nCov)

    # start processes for affinity model #
    parsedSequenceLength = parseSequenceLength(length)
    (
        affinityProcessQ, affinityContainer, affinityFlag,
        affinityNeighborIndices1, affinityNeighborAffinities1, affinityNeighborIndices2, affinityNeighborAffinities2,
        affinityProcesses
    ) = createAffinityProcesses(
        sharedKmerFreq, sharedKmerFreqSquaredRowSum, parsedSequenceLength, 136, parameters.min_kmer_frequence_probability,
        sharedCovMean, sharedCovVar, nCov, parameters.min_coverage_probability,
        nModelWeights, modelWeights, nSequences, parameters.sequence_neighbors, nSeeds, parameters.seed_neighbors, parameters.threads
    )

    # start processes for lp model #
    lpPartitionN = nModelWeights * nMarkers
    lpProcessQ, lpPartitions, lpProcesses = createLPProcesses(lpPartitionN, nSeeds, parameters.threads)

    # start process for score model #
    sharedSeedSequences, seedSequences, sharedSeedMarkers, seedMarkers = createSeedMarker(sequence2markers)
    scoreProcessQ, scorePartitions, scoreScores, scoreProcesses = createScoreProcesses(
        sharedSeedSequences, sharedSeedMarkers, lpPartitions, lpPartitionN, nModelWeights
    )

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Clustering sequences.', flush = True)
    # initialize the labels #
    labels = numpy.zeros(shape = nSequences, dtype = numpy.int32)
    labelMask = numpy.empty(shape = nSequences, dtype = numpy.bool_)
    processedNSequences = 0
    processBar = ProcessBar(nSequences)

    while processedNSequences < nSequences:
        for label in numpy.unique(labels[labels >= 0]):
            numpy.equal(labels, label, out = labelMask)
            sequences = numpy.flatnonzero(labelMask)
            seedSequencesMask = numpy.isin(seedSequences, sequences, assume_unique = True)
            if sequences.size == nSequences: # the 1st cluster for seeds #
                seedSequencesMask[length[seedSequences[seedSequencesMask]] < (1500 if nCov == 1 else 2500)] = False
            seeds = seedSequences[seedSequencesMask][ : parameters.max_seeds]
            score, recall, nClusters, seedSets = evaluate_cluster(seedMarkers[seedSequencesMask][ : parameters.max_seeds])
            if recall < 0.2 or numpy.sum(length[sequences]) < parameters.min_cluster: # clusters with recall < 0.2 or size < min_cluster will be removed #
                labels[sequences] = -nSequences - 1
            elif nClusters <= 1 or score >= 0.95:
                labels[sequences] = -sequences[0] - 1
            else: # nClusters > 1 and recall >= 0.2 #
                seedLabels, neighborProbability = clusterSeeds(
                    affinityProcessQ, affinityContainer, affinityFlag,
                    affinityNeighborIndices2, affinityNeighborAffinities2,
                    lpProcessQ, scoreProcessQ, scorePartitions, scoreScores,
                    seeds, seedSets, parameters.seed_neighbors, nModelWeights, parameters.threads
                ) # >0 for high quality, =0 for unprocessed, <0 for low quality #
                flag = True
                if ((nCov > 1 and nClusters >= parameters.min_swdpgmm_clusters) or (sequences.size >= parameters.min_swdpgmm_sequences) or (seeds.size == parameters.max_seeds)) and numpy.any(seedLabels):
                    y[ : ] = 0
                    for i in numpy.unique(seedLabels):
                        y[seeds[seedLabels == i]] = i
                    swdpgmm = SWDPGMM(minW = 2, maxIterations = parameters.max_swdpgmm_iterations)
                    z = swdpgmm.main(x[sequences], length[sequences] * 1e-3, y[sequences])
                    uniqueZ = numpy.unique(z)
                    if uniqueZ.size > 1:
                        for z_ in uniqueZ:
                            sequences_ = sequences[z == z_]
                            labels[sequences_] = sequences_[0]
                        flag = False
                if flag:
                    labels[sequences] = clusterSequences(
                        affinityProcessQ, affinityContainer, affinityFlag,
                        affinityNeighborIndices1, affinityNeighborAffinities1,
                        sequences, seeds, seedLabels, nCov == 1 and sequences.size == nSequences,
                        parameters.sequence_neighbors, neighborProbability, parameters.threads
                    )
        processedNSequences = numpy.count_nonzero(labels < 0)
        processBar.plot(processedNSequences)
    freeAffinityProcesses(affinityProcessQ, affinityProcesses)
    freeLPProcesses(lpProcessQ, lpProcesses)
    freeScoreProcess(scoreProcessQ, scoreProcesses)
    output_clusters(sequenceStructs, labels, parameters.output, parameters.no_fasta_clusters)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
