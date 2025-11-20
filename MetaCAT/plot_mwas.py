from datetime import datetime
import gzip
from math import floor, log10

import numpy
from matplotlib.pyplot import close, subplots

from .colors import generateColors


def isGzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def readVariantAnnotationFile(file):
    x = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Classification\tChromosome\tPosition\tAnnotation', f'\"{file}\" is not a valid variant annotation file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        x[(lines[1], int(lines[2]))] = lines[3]
    openFile.close()
    return x


def readMwasFile(file):
    classification2chromosome2xyz = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().startswith('Classification\tChromosome\tPosition\tAllele\tBeta\tSE\tP\tSignificance'), f'\"{file}\" is not a valid mwas file.'
    n = 0
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        flag = 2 if lines[7] == '-' else 1
        classification2chromosome2xyz.setdefault(lines[0], dict()).setdefault(lines[1], list()).append((int(lines[2]), float(lines[6]), flag))
        n += 1
    openFile.close()
    return (n, classification2chromosome2xyz)


def qq(y, output, width, height):
    y = numpy.sort(y)
    x = numpy.arange(y.size, 0, -1, dtype = numpy.float64)
    x -= 0.5
    x /= x.size
    numpy.log10(x, out = x)
    numpy.negative(x, out = x)

    xMin = x[0]
    xMax = x[-1] + 0.02 * (x[-1] - x[0])
    yMax = y[-1] + 0.5

    figure, axes = subplots(figsize = (width, height), squeeze = True)
    axes.scatter(x, y, c = '#666699', s = 5, alpha = 1, marker = 'o', edgecolors = 'none', rasterized = True)
    axes.set_xlim(xMin, xMax)
    axes.set_ylim(0, yMax)
    axes.plot([xMin, xMax], [xMin, xMax], color = '#888888', linewidth = 1, linestyle = '--')
    axes.set_xlabel(r'Expected $-log_{10}{(p)}$', fontdict = {'fontsize': 9})
    axes.set_xticks(numpy.arange(floor(xMax) + 1), numpy.arange(floor(xMax) + 1), size = 8)
    axes.set_ylabel(r'Observed $-log_{10}{(p)}$', fontdict = {'fontsize': 9})
    axes.set_yticks(numpy.arange(0, floor(yMax) + 1, 2), numpy.arange(0, floor(yMax) + 1, 2), size = 8)
    figure.savefig(output, dpi = 600, bbox_inches = 'tight')
    close()
    return None


def readAnnotationFile(file):
    x = list()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Cluster ID\tClassification', f'\"{file}\" is not a valid annotation file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        x.append((lines[1], lines[0]))
    openFile.close()
    x.sort()
    return x


def manhattan(xyz, variantAnnotations, offsets, alpha, output, width, height):
    figure, axes = subplots(2, 1, height_ratios = [19, 1], figsize = (width, height), squeeze = False)
    xy = xyz[0 : 2, xyz[2] == 2]
    axes[0, 0].scatter(*xy, c = '#CCCCCC', s = 8, alpha = 1, marker = 'o', rasterized = True)
    xy = xyz[0 : 2, xyz[2] == 1]
    axes[0, 0].scatter(*xy, c = '#888888', s = 8, alpha = 1, marker = 'o', rasterized = True)
    axes[0, 0].axhline(y = alpha, color = '#DD4A4A', linewidth = 1, linestyle = '--')
    xy = xyz[0 : 2, xyz[2] == 0]
    for i, j in zip(xy.T, variantAnnotations):
        axes[0, 0].annotate(
            j, xy = i, xytext = (i[0] + 5000000, i[1] + 0.9),
            arrowprops = dict(arrowstyle = '-', color = '#DD4A4A', shrinkA = 2, shrinkB = 2),
            size = 6,
            color = '#DD4A4A',
            horizontalalignment = 'left',
            verticalalignment = 'bottom'
        )
    axes[0, 0].scatter(*xy, c = '#DD4A4A', s = 8, alpha = 1, marker = 'o', rasterized = True)
    axes[0, 0].xaxis.set_visible(False)
    axes[0, 0].set_xlim(-0.01 * offsets[-1, 1], 1.01 * offsets[-1, 1])
    axes[0, 0].set_ylim(0, numpy.max(xyz[1]) + 0.5)
    axes[0, 0].set_ylabel('$-log_{10}{(p)}$', fontdict = {'fontsize': 9})
    axes[0, 0].tick_params(labelsize = 8)

    axes[1, 0].barh(0.1, offsets[ : , 1] - offsets[ : , 0], height = 0.16, left = offsets[ : , 0], color = generateColors(offsets.shape[0]), edgecolor = '#111111', linewidth = 0.05)
    axes[1, 0].set_xticks([])
    axes[1, 0].yaxis.set_visible(False)
    axes[1, 0].spines[ : ].set_visible(False)
    axes[1, 0].set_ylim(0, 0.2)
    axes[1, 0].set_xlim(-0.01 * offsets[-1, 1], 1.01 * offsets[-1, 1])
    axes[1, 0].set_xlabel('Species', fontdict = {'fontsize': 9})
    figure.tight_layout(h_pad = 0.1)
    figure.savefig(output, dpi = 600, bbox_inches = 'tight')
    close(figure)
    return None


def main(parameters):
    if parameters.variant_annotation is None:
        sequencePosition2annotation = dict()
    else:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading variant annotation file.', flush = True)
        sequencePosition2annotation = readVariantAnnotationFile(parameters.variant_annotation)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading mwas file.', flush = True)
    n, classification2chromosome2xyz = readMwasFile(parameters.mwas)

    offsets = numpy.empty(shape = (len(classification2chromosome2xyz), 2), dtype = numpy.float64)
    offset = 0
    xyz = numpy.empty(shape = (3, n), dtype = numpy.float64)
    variantAnnotations = list()
    I = 0
    for i, classification in enumerate(sorted(classification2chromosome2xyz.keys())):
        offsets[i, 0] = offset
        chromosome2xyz = classification2chromosome2xyz[classification]
        for chromosome in sorted(chromosome2xyz.keys()):
            chromosome2xyz[chromosome].sort()
            for j in chromosome2xyz[chromosome]:
                if (chromosome, j[0]) in sequencePosition2annotation:
                    variantAnnotations.append(sequencePosition2annotation[(chromosome, j[0])])
                    xyz[ : , I] = [offset + j[0], j[1], 0]
                else:
                    xyz[ : , I] = [offset + j[0], j[1], j[2]]
                I += 1
            offset += j[0]
        offsets[i, 1] = offset
    numpy.clip(xyz[1], 10 ** -parameters.max_negative_log_p, 1, out = xyz[1])
    numpy.log10(xyz[1], out = xyz[1])
    numpy.negative(xyz[1], out = xyz[1])

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Creating qq plot.', flush = True)
    qq(xyz[1], parameters.output + '.qq.pdf', parameters.qq_width, parameters.qq_height)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Creating manhattan plot.', flush = True)
    manhattan(xyz, variantAnnotations, offsets, -log10(parameters.alpha / n), parameters.output + '.manhattan.pdf', parameters.manhattan_width, parameters.manhattan_height)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
