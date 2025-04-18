import gzip
import os
import sys
from ctypes import c_float, c_int32
from datetime import datetime
from math import ceil, log10
from multiprocessing.sharedctypes import RawArray
from operator import itemgetter

import numpy
from sklearn.decomposition import PCA

from .affinity_model import createProcesses as createAffinityProcesses
from .affinity_model import freeProcesses as freeAffinityProcesses
from .fasta import read_fasta_file
from .kmer import count_kmers
from .kmer_frequency_model import parseW as parseSequenceLength
from .label_propagation import createProcesses as createLPProcesses
from .label_propagation import freeProcesses as freeLPProcesses
from .marker import markerHash, nMarkers
from .processbar import ProcessBar
from .score_model import createMarkerTable
from .score_model import createProcesses as createScoreProcesses
from .score_model import evaluate_sequences
from .score_model import freeProcesses as freeScoreProcess
from .semi_supervised_model import seedModel as clusterSeeds
from .semi_supervised_model import sequenceModel as clusterSequences

maxNModelWeights = 5


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
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[0] in markerHash:
            marker = markerHash[lines[0]]
            for sequenceID in lines[1 : ]:
                sequence2markers.setdefault(sequenceID, list()).append(marker)
    openFile.close()

    for sequence, markers in sequence2markers.items():
        sequence2markers[sequence] = numpy.sort(markers)
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


def parseModelWeight(ws1, ws2):
    if len(ws1) > maxNModelWeights or len(ws2) > maxNModelWeights:
        ws1 = [0.9, 0.7, 0.5, 0.3, 0.1]
        ws2 = [0.1, 0.3, 0.5, 0.7, 0.9]
    ws1.sort(reverse = True)
    ws2.sort()
    w = RawArray(c_float, len(ws1) + len(ws2))
    for i, (w1, w2) in enumerate(zip(ws1, ws2)):
        w[2 * i] = w1
        w[2 * i + 1] = w2
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


def read_mapping_file(input_file): # only for debug #
    sequence_id2genome_id = dict()
    genome_id2genome = dict()
    genome_sizes = list()
    genome = 0
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        if not line.startswith('@'):
            lines = line.rstrip('\n').split('\t')
            # sequence genome length #
            sequence_id2genome_id[lines[0]] = lines[1]
            if lines[1] not in genome_id2genome:
                genome_id2genome[lines[1]] = genome
                genome_sizes.append(0)
                genome += 1
            genome_ = genome_id2genome[lines[1]]
            genome_sizes[genome_] += int(lines[2]) if len(lines) >= 3 else 0
    open_file.close()
    genome_sizes = numpy.asarray(genome_sizes, dtype = numpy.float32)
    return (genome_id2genome, sequence_id2genome_id, genome_sizes)


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
    randomGenerator = numpy.random.default_rng(parameters.random_number)
    nCov = determine_samples(parameters.coverage)
    if parameters.min_sequence_length is None:
        parameters.min_sequence_length = 1500 if nCov == 1 else 500

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
    for sequenceID, SEQUENCE in read_fasta_file(parameters.fasta):
        length_ = len(SEQUENCE)
        if length_ >= parameters.min_sequence_length or sequenceID in sequenceID2markers:
            sequenceStructs.append((sequenceID, SEQUENCE, length_))
    sequenceStructs.sort(key = itemgetter(2), reverse = True)
    nSequences = len(sequenceStructs)

    sequenceID2sequence = dict()
    length = numpy.empty(shape = nSequences, dtype = numpy.float32)
    for i, (sequenceID, SEQUENCE, length_) in enumerate(sequenceStructs):
        sequenceID2sequence[sequenceID] = i
        length[i] = length_

    sequence2markers = dict()
    for sequenceID, markers in sequenceID2markers.items():
        sequence2markers[sequenceID2sequence[sequenceID]] = markers
    del sequenceID2markers

    # sequence length #
    length = numpy.array(length, dtype = numpy.float32)

    # read coverage file #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading coverage file.', flush = True)
    covMean, sharedCovMean, sharedCovVar = read_coverage_file(parameters.coverage, sequenceID2sequence)

    # calculate kmer frequency #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Counting kmers of all sequences.', flush = True)
    kmerFreq, sharedKmerFreq = count_kmers(nSequences, sequenceStructs, parameters.threads)
    '''
    if os.access(f'{os.path.basename(parameters.fasta)}.{parameters.min_sequence_length}.kmers.npy', os.R_OK):
        sharedKmerFreq = RawArray(c_float, numpy.load(f'{os.path.basename(parameters.fasta)}.{parameters.min_sequence_length}.kmers.npy').flat)
        kmerFreq = numpy.ndarray(shape = (nSequences, 136), dtype = numpy.float32, buffer = sharedKmerFreq)
    else:
        kmerFreq, sharedKmerFreq = count_kmers(nSequences, sequenceStructs, parameters.threads)
        numpy.save(f'{os.path.basename(parameters.fasta)}.{parameters.min_sequence_length}.kmers.npy', kmerFreq)
    '''

    nModelWeights, modelWeights = parseModelWeight(parameters.kmer_frequence_weight, parameters.coverage_weight)

    # start processes for affinity model #
    parsedSequenceLength = parseSequenceLength(length)
    sharedKmerFreqSquaredRowSum = RawArray(c_float, numpy.einsum('ij,ij->i', kmerFreq, kmerFreq))
    (
        affinityProcessQ, affinityContainer, affinityContainerSize, affinityToken, affinityModelFlag,
        affinityNeighborIndices1, affinityNeighborAffinities1, affinityNeighborIndices2, affinityNeighborAffinities2,
        affinityProcesses
    ) = createAffinityProcesses(
        sharedKmerFreq, sharedKmerFreqSquaredRowSum, parsedSequenceLength, 136, parameters.min_kmer_frequence_probability,
        sharedCovMean, sharedCovVar, nCov, parameters.min_coverage_probability,
        nModelWeights, modelWeights, nSequences, nSeeds, parameters.neighbors, parameters.threads
    )

    # start processes for lp model #
    lpPartitionN = nModelWeights * nMarkers
    lpProcessQ, lpPartitions, lpProcesses = createLPProcesses(lpPartitionN, nSeeds, parameters.threads)

    # start process for score model #
    markerTable = createMarkerTable(sequence2markers) # (nSeeds, nMarkers) #
    sharedSeedSequences = RawArray(c_int32, sorted(sequence2markers.keys()))
    scoreProcessQ, scorePartitions, scoreScores, scoreProcesses = createScoreProcesses(
        sharedSeedSequences, markerTable, lpPartitions, lpPartitionN, nModelWeights, parameters.threads
    )

    covMean = numpy.log(covMean + 1) # covMean >= 0 #
    covMean -= numpy.min(covMean, axis = 0)
    covMean /= numpy.max(covMean, axis = 0) + 10 * numpy.finfo(numpy.float32).eps # covMean >= 0 and covMean <= 1 #
    pca = PCA(n_components = None, whiten = False, random_state = randomGenerator.integers(numpy.iinfo(numpy.int32).max))
    pca.fit(kmerFreq)
    k = max(numpy.sum(numpy.cumsum(pca.explained_variance_ratio_) < parameters.min_pca_variance_ratio), parameters.min_pca_components)
    x = numpy.concatenate((pca.transform(kmerFreq)[ : , : k], covMean), axis = 1)
    y = numpy.empty(shape = nSequences, dtype = numpy.int32)

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
            score, recall, nClusters, seeds, seedSets = evaluate_sequences(
                sequences[length >= (4000 if nCov == 1 else 2500)] if sequences.size == nSequences else sequences, # sorted sequences #
                sequence2markers, # sorted markers #
                parameters.max_seeds
            )
            if recall < 0.2 or numpy.sum(length[sequences]) < parameters.min_cluster:
                labels[sequences] = -nSequences - 1
            elif nClusters <= 1 or score >= 0.95: # all clusters with the number of genomes is one will be output #
                labels[sequences] = -sequences[0] - 1
            else: # nClusters > 1 and recall >= 0.2 #
                seedLabels = clusterSeeds(
                    affinityProcessQ, affinityContainer, affinityContainerSize, affinityToken,
                    affinityNeighborIndices2, affinityNeighborAffinities2, affinityModelFlag,
                    lpProcessQ, scoreProcessQ, scorePartitions, scoreScores,
                    seeds, seedSets, parameters.neighbors, nModelWeights
                )
                if (nClusters >= parameters.min_swdpgmm_clusters) or (sequences.size >= parameters.min_swdpgmm_sequences) or (seeds.size == parameters.max_seeds):
                    y[ : ] = -1
                    for i in numpy.unique(seedLabels):
                        y[seeds[seedLabels == i]] = i
                    swdpgmm = SWDPGMM(minW = 2, maxIterations = parameters.max_swdpgmm_iterations)
                    z = swdpgmm.main(x[sequences], length[sequences] * 1e-3, y[sequences])
                    uniqueZ = numpy.unique(z)
                    if uniqueZ.size > 1:
                        for z_ in uniqueZ:
                            sequences_ = sequences[z == z_]
                            labels[sequences_] = sequences_[0]
                        continue
                labels[sequences] = clusterSequences(
                    affinityProcessQ, affinityContainer, affinityContainerSize, affinityToken,
                    affinityNeighborIndices1, affinityNeighborAffinities1,
                    sequences, seeds, seedLabels, parameters.neighbors
                )
        processedNSequences = numpy.count_nonzero(labels < 0)
        processBar.plot(processedNSequences)
    freeAffinityProcesses(affinityProcessQ, affinityProcesses)
    freeLPProcesses(lpProcessQ, lpProcesses)
    freeScoreProcess(scoreProcessQ, scoreProcesses)
    output_clusters(sequenceStructs, labels, parameters.output, parameters.no_fasta_clusters)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
