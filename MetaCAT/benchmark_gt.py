import gzip
import os

import numpy
from matplotlib.pyplot import close, subplots

from .colors import parseRGBA


def isGzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def readGTFile(file):
    sequence2index = dict()
    sequenceLengths = list()
    genomeLengths = list()
    genomes = list()
    nSequences = 0
    nGenomes = 0
    genome2index = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Sequence ID\tGenome ID\tLength', f'\"{file}\" must have a single header line: Sequence ID<tab>Genome ID<tab>Length.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[1] not in genome2index:
            genome2index[lines[1]] = nGenomes
            nGenomes += 1
        sequence2index[lines[0]] = nSequences
        nSequences += 1
        sequenceLengths.append(int(lines[2]))
        genomes.append(genome2index[lines[1]])
    openFile.close()
    sequenceLengths = numpy.asarray(sequenceLengths, dtype = numpy.float64)
    genomes = numpy.asarray(genomes, dtype = numpy.int64)
    genomeLengths = numpy.empty(shape = nGenomes, dtype = numpy.float64)
    for genome in range(nGenomes):
        genomeLengths[genome] = numpy.sum(sequenceLengths, where = genomes == genome)
    return (sequence2index, sequenceLengths, genomes, genomeLengths, nGenomes)


def readMappingFile(sequence2index, file):
    cluster2sequences = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Sequence ID\tCluster ID', f'\"{file}\" must have a single header line: Sequence ID<tab>Cluster ID.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        sequence = sequence2index[lines[0]]
        cluster2sequences.setdefault(lines[1], list()).append(sequence)
    openFile.close()
    return list(cluster2sequences.values())


def calculateScores(sequenceLengths, genomeLengths, genomes, nGenomes, clusterSequences):
    nClusters = len(clusterSequences)
    scores = numpy.empty(shape = (nClusters, 2), dtype = numpy.float64)
    benchmarks = numpy.zeros(shape = (nClusters, nGenomes), dtype = numpy.float64) # nClusters, nGenomes #
    for i, sequences in enumerate(clusterSequences):
        for sequence in sequences:
            benchmarks[i, genomes[sequence]] += sequenceLengths[sequence]
    clusterLength = numpy.sum(benchmarks, axis = 1)
    i = numpy.argmax(benchmarks / genomeLengths, axis = 1)
    _ = benchmarks[numpy.arange(nClusters), i]
    numpy.divide(_, genomeLengths[i], out = scores[ : , 0])
    numpy.divide(_, clusterLength, out = scores[ : , 1])
    return scores


def scatterPlot(dataset, programs, program2scores, output, color):
    nPrograms = len(program2scores)
    figure, axes = subplots(nrows = 1, ncols = nPrograms, figsize = (0.2 + 1.1 * nPrograms, 1.58), layout = 'constrained', squeeze = False)
    for subplot, program in zip(axes[0], programs):
        scores = program2scores[program]
        subplot.scatter(scores[ : , 1], scores[ : , 0], s = 10, color = color, marker = 'o', edgecolors = 'none', rasterized = True)
        subplot.set_xticks([0, 0.5, 1])
        subplot.set_xlim(-0.05, 1.05)
        subplot.set_yticks([0, 0.5, 1])
        subplot.set_ylim(-0.05, 1.05)
        subplot.tick_params(labelsize = 8)
        #subplot.xaxis.set_minor_locator(AutoMinorLocator(2))
        #subplot.yaxis.set_minor_locator(AutoMinorLocator(2))
        #subplot.grid(True, which = 'both', linestyle = '--', linewidth = 0.5)
        subplot.set_title(program, fontdict = {'fontsize': 9})
        #subplot.set_xlabel('Precision', fontdict = {'fontsize': 9})
        #subplot.set_ylabel('Recall', fontdict = {'fontsize': 9})
    figure.suptitle(dataset, size = 10)
    figure.supxlabel('Precision', size = 9)
    figure.supylabel('Recall', size = 9)
    figure.savefig(output, dpi = 600, bbox_inches = 'tight')
    close()
    return None


def barPlot(precisionSet, recallSet, dataset, programs, program2scores, output, colors):
    nPrecisionSet = len(precisionSet)
    nRecallSet = len(recallSet)
    figure, axes = subplots(
        nrows = 1,
        ncols = nPrecisionSet,
        sharey = True,
        figsize = (3.5 * nPrecisionSet, 0.6 + 0.2 * len(programs)),
        layout = 'constrained',
        squeeze = False
    )
    x = numpy.zeros(shape = 2 * nRecallSet, dtype = numpy.int64)
    for subplot, precision in zip(axes.flat, precisionSet):
        for program in programs:
            scores = program2scores[program]
            precisionMask = scores[ : , 1] >= precision
            for i, recall in enumerate(recallSet):
                x[i] = numpy.count_nonzero((scores[ : , 0] >= recall) & precisionMask)
                yield (program, precision, recall, x[i])
            x[nRecallSet + 1 : ] = x[ : nRecallSet - 1]
            x[ : nRecallSet] -= x[nRecallSet : ]
            subplot.barh(program, x[ : nRecallSet], left = x[nRecallSet : ], alpha = None, height = 0.6, color = colors, linewidth = 0, label = recallSet)
        subplot.set_title(dataset, fontdict = {'fontsize': 10})
        subplot.set_xlabel(f'Number of genomes (Precision ≥ {precision:.2f})', fontdict = {'fontsize': 9})
        subplot.tick_params(labelsize = 8)
        #subplot.xaxis.set_minor_locator(AutoMinorLocator(2))
        #subplot.grid(True, which = 'both', axis = 'x', linestyle = '--', linewidth = 0.5)
    figure.legend(
        recallSet,
        ncol = nRecallSet,
        title = 'Recall',
        bbox_to_anchor = (0.5, 1.0),
        loc = 'lower center',
        fontsize = 8,
        title_fontsize = 9
    )
    figure.savefig(output, dpi = 300, bbox_inches = 'tight')
    close()
    return None


def main(parameters):
    if parameters.dataset is None:
        parameters.dataset = os.path.splitext(os.path.basename(parameters.ground_truth))[0]
    if parameters.label is None:
        parameters.label = list()
        for i in parameters.mapping:
            parameters.label.append(os.path.basename(i))
    parameters.precision.sort(reverse = True)
    parameters.recall.sort(reverse = True)
    sequence2index, sequenceLengths, genomes, genomeLengths, nGenomes = readGTFile(parameters.ground_truth)
    label2scores = dict()
    for mappingFile, label in zip(parameters.mapping, parameters.label):
        sequenceSets = readMappingFile(sequence2index, mappingFile)
        label2scores[label] = calculateScores(sequenceLengths, genomeLengths, genomes, nGenomes, sequenceSets)
    colors = parseRGBA(parameters.color, len(parameters.recall), parameters.color_min, parameters.color_max)
    scatterPlot(parameters.dataset, parameters.label, label2scores, parameters.output + '.quality.pdf', (*colors[0][ : 3], 0.1))
    program2clusters = dict()
    for program, precision, recall, x in barPlot(parameters.precision, parameters.recall, parameters.dataset, parameters.label[ : : -1], label2scores, parameters.output + '.pdf', colors):
        program2clusters.setdefault(program, list()).append(str(x))
    openFile = open(parameters.output + '.tsv', 'w')
    header = ['Dataset', 'Program']
    for precision in parameters.precision:
        for recall in parameters.recall:
            header.append(f'Precision: {precision:.2f} Recall: {recall:.2f}')
    openFile.write('\t'.join(header) + '\n')
    for program, clusters in program2clusters.items():
        clusters = '\t'.join(clusters)
        openFile.write(f'{parameters.dataset}\t{program}\t{clusters}\n')
    openFile.close()
    return None
