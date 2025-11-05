import gzip
from itertools import combinations
from math import ceil, log10, sqrt

import numpy
from matplotlib.pyplot import close, subplots
from scipy.stats import false_discovery_control, mannwhitneyu

zeroAbundance = 1e-5
log10ZeroAbundance = log10(zeroAbundance)


def isGzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def readGroupFile(file):
    individual2group = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'ID\tGroup', f'\"{file}\" must have a single header line: ID<tab>Group.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        individual2group[lines[0]] = lines[1]
    openFile.close()
    return individual2group


def readAbundanceFile(file, classifications):
    nodes = list()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    header = openFile.readline().rstrip('\n')
    assert header.startswith('Abundance'), f'\"{file}\" must have a single header line: Abundance ... .'
    individuals = header.split('\t')[1 : ]
    abundance = list()
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[0].rsplit(';', maxsplit = 1)[-1][0] in classifications:
            nodes.append(lines[0])
            abundance.append(lines[1 : ])
    openFile.close()
    abundance = numpy.asarray(abundance, dtype = numpy.float64)
    numpy.clip(abundance, zeroAbundance, 1, out = abundance)
    numpy.log10(abundance, out = abundance)
    return (individuals, nodes, abundance)


def test(x1, x2, coverage, alpha, multipleTest):
    '''
    x1: (classifications, individuals)
    x2: (classifications, individuals)
    '''
    mask1 = numpy.count_nonzero(x1 > log10ZeroAbundance, axis = 1) >= coverage * x1.shape[1]
    mask2 = numpy.count_nonzero(x2 > log10ZeroAbundance, axis = 1) >= coverage * x2.shape[1]
    mask3 = numpy.flatnonzero(mask1 & mask2)
    p = numpy.full(shape = x1.shape[0], fill_value = numpy.nan, dtype = numpy.float64)
    p[mask3] = mannwhitneyu(x1[mask3], x2[mask3], axis = 1).pvalue
    if multipleTest == 'bonferroni':
        index = numpy.flatnonzero(p <= alpha / mask3.size)
    elif multipleTest == 'benjamini-hochberg':
        p[mask3] = false_discovery_control(p[mask3], method = 'bh')
        index = numpy.flatnonzero(p <= alpha)
    else:
        index = numpy.flatnonzero(p <= alpha)
    return (p, index[numpy.argsort(p[index])])


def boxPlot(nodes, rows, columns, x1, x2, p, group1, group2, output):
    '''
    nodes: list()
    x1: (classifications, individuals)
    x2: (classifications, individuals)
    p: (classifications, )
    '''
    #randomGenerator = numpy.random.default_rng(0)
    figure, axes = subplots(
        nrows = rows,
        ncols = columns,
        figsize = (2.05 * columns, 1.15 * rows),
        layout = 'constrained',
        squeeze = False
    )
    for subplot in axes.flat[len(nodes) : ]:
        subplot.remove()
    for i, (node, subplot) in enumerate(zip(nodes, axes.flat)):
        subplot.boxplot(
            [x1[i], x2[i]],
            orientation = 'horizontal',
            notch = True,
            bootstrap = 1000,
            medianprops = dict(color = 'red'),
            widths = 0.6,
            patch_artist = False,
            tick_labels = [group1, group2],
            showfliers = False,
            flierprops = dict(marker = '.', markerfacecolor = 'gray', markersize = 1, markeredgecolor = 'none')
        )
        '''
        subplot.scatter(
            x1[i], 0.75 + 0.5 * randomGenerator.random(x1[i].size, dtype = numpy.float64),
            s = 1,
            c = 'gray',
            marker = '.'
        )
        subplot.scatter(
            x2[i], 1.75 + 0.5 * randomGenerator.random(x2[i].size, dtype = numpy.float64),
            s = 1,
            c = 'gray',
            marker = '.'
        )
        '''
        subplot.set_title(f'{node.rsplit(";", maxsplit = 1)[-1][3 : ]}\n$p$: {p[i]:.2e}', fontdict = {'fontsize': 9})
        subplot.set_xlabel('$log_{10}{(RA)}$', fontdict = {'fontsize': 9})
        subplot.tick_params(labelsize = 8, pad = 0.01)
    figure.savefig(output, dpi = 600, bbox_inches = 'tight')
    close(figure)
    return None


def main(parameters):
    individual2group = readGroupFile(parameters.group)
    individuals, nodes, x = readAbundanceFile(parameters.abundance, set(parameters.classification))
    group2indices = dict()
    for i, individual in enumerate(individuals):
        if individual in individual2group:
            group2indices.setdefault(individual2group[individual], list()).append(i)
    for group1, group2 in combinations(sorted(group2indices.keys()), 2):
        if (parameters.comparison is not None) and ((group1 not in parameters.comparison) or (group2 not in parameters.comparison)):
            continue
        x1 = x[ : , group2indices[group1]] # (classifications, individuals) #
        x2 = x[ : , group2indices[group2]] # (classifications, individuals) #
        p, index = test(x1, x2, parameters.coverage, parameters.alpha, parameters.multiple_test)
        if index.size:
            columns = parameters.columns
            rows = parameters.rows
            if rows is None and columns is None:
                columns = ceil(sqrt(index.size))
                rows = ceil(index.size / columns)
            elif rows is None:
                rows = ceil(index.size / columns)
            elif columns is None:
                columns = ceil(index.size / rows)
            else:
                assert rows * columns >= index.size, 'Size of grid must be greater than or equal to the number of plots.'
            boxPlot(
                [nodes[i] for i in index], rows, columns, x1[index], x2[index], p[index],
                group1, group2, f'{parameters.output}.{group1}-{group2}.pdf'
            )
        for i, j, k, l in zip(nodes, p, x1, x2):
            print(i, j, group1, numpy.median(k), group2, numpy.median(l), sep = '\t', flush = True)
