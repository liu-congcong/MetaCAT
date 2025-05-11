from collections.abc import Mapping, Sequence
import gzip
import os
from ctypes import c_double, c_float, c_int64
from datetime import datetime
from functools import partial
from multiprocessing import Lock, Process, Queue
from multiprocessing.sharedctypes import RawArray, Value
from uuid import uuid4

import numpy
from numpy.typing import NDArray
from sklearn.decomposition import PCA
from statsmodels.genmod.families.family import Binomial, Gaussian
from statsmodels.genmod.generalized_linear_model import GLM
from threadpoolctl import threadpool_limits

from .processbar import ProcessBar

mapping1 = ('A', 'C', 'G', 'T')
mapping2 = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def isGzipped(file):
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()
    return magicCode == b'\x1f\x8b'


def readPhenotypeFile(file, phenotypeColumn, id2index):
    phenotypeColumn = phenotypeColumn - 1
    n = len(id2index)
    flag = numpy.zeros(shape = n, dtype = numpy.bool_)
    phenotypes = numpy.empty(shape = n, dtype = numpy.float32)
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().startswith('ID'), f'\"{file}\" is not a valid phenotype file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if (lines[0] in id2index) and (lines[phenotypeColumn].lower() != 'na'):
            flag[id2index[lines[0]]] = True
            phenotypes[id2index[lines[0]]] = float(lines[phenotypeColumn])
    openFile.close()
    return (flag, phenotypes)


def readAnnotationFile(file: str) -> tuple[Sequence[str], Mapping[str, int]]:
    classifications = list()
    classification2index = dict()
    cluster2classification = dict()
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
    for i, j in enumerate(x):
        assert j[0] not in classification2index, f'\"{j[0]}\" is not unique.'
        classifications.append(j[0])
        cluster2classification[j[1]] = i
    return (classifications, cluster2classification)


def readMappingFile(file: str, cluster2classification: Mapping[str, int]) -> Mapping[str, int]:
    sequence2classification = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().rstrip('\n') == 'Sequence ID\tCluster ID', f'\"{file}\" is not a valid mapping file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        sequence2classification[lines[0]] = cluster2classification[lines[1]]
    openFile.close()
    return sequence2classification


def readAbundanceFile(file: str, classifications: Sequence[str], individuals: Sequence[str]) -> NDArray[numpy.float32]:
    n = len(classifications)
    m = len(individuals)
    X = RawArray(c_float, n * m) # offset + n * m #
    x = numpy.ndarray(shape = (n, m), dtype = numpy.float32, buffer = X)
    classification2index = dict((j, i) for i, j in enumerate(classifications))
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    headers = openFile.readline().rstrip('\n').split('\t')
    assert headers[0] == 'Abundance', f'\"{file}\" is not a valid abundance file.'
    id2index = dict((j, i) for i, j in enumerate(headers[1 : ], start = 1))
    indices = [id2index[i] for i in individuals]
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[0] in classification2index:
            x[classification2index[lines[0]]] = [float(lines[i]) for i in indices]
    openFile.close()
    return X


def readCovariateFile(file, id2index):
    n = len(id2index)
    flag = numpy.zeros(shape = n, dtype = numpy.bool_)
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    line = openFile.readline()
    assert line.startswith('ID'), f'\"{file}\" is not a valid covariate file.'
    m = line.count('\t')
    x = numpy.empty(shape = (n, m), dtype = numpy.float32)
    for line in openFile:
        id, covariates = line.rstrip('\n').split('\t', maxsplit = 1)
        if (id in id2index) and ('na' not in covariates.lower()):
            flag[id2index[id]] = True
            x[id2index[id]] = covariates.split('\t')
    openFile.close()
    return (flag, x)


def createADM(x):
    y = x.T @ x
    y *= -2
    xx = numpy.einsum('ij,ij->j', x, x, dtype = numpy.float64)
    y += xx
    y += xx[ : , None]
    y.flat[ : : xx.size + 1] = 0
    y **= 0.5
    return y


def readMWASFile(file):
    chrPos2alleleBeta = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().startswith('Classification\tChromosome\tPosition\tAllele\tBeta\tSE\tP\tSignificance'), f'\"{file}\" is not a valid mwas file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        beta = float(lines[4])
        if lines[7] == '+':
            chrPos2alleleBeta[(lines[1], lines[2])] = (mapping2[lines[3]], beta)
    openFile.close()
    return chrPos2alleleBeta


def testProcess(queue, file, n, m, phenotypes, covariates, abundances, variantColumns, chrPos2alleleBeta):
    m0 = numpy.count_nonzero(phenotypes == 0)
    m1 = numpy.count_nonzero(phenotypes == 1)
    if m0 + m1 == phenotypes.size:
        flag = 0
        Glm = partial(GLM, family = Binomial())
    else:
        flag = 1
        Glm = partial(GLM, family = Gaussian())
    m0 *= m
    m1 *= m
    m *= phenotypes.size

    N = 0
    x = numpy.hstack((numpy.ones(shape = (phenotypes.size, 2), dtype = numpy.float32), covariates))
    mask = numpy.empty(shape = phenotypes.size, dtype = numpy.bool_)
    acgt = numpy.empty(shape = 4, dtype = numpy.float32)
    openFile = open(file, 'w')
    # sequence position major beta se p gt1 gt2 ... #
    numpy.seterr(all = 'raise')
    while True:
        abundanceOffset, sequence, position, populationGenotype, individualGenotypes = queue.get()
        if sequence is None:
            break
        if chrPos2alleleBeta is None:
            acgt[ : ] = populationGenotype.split('/')
            allele = numpy.argmax(acgt)
        elif (sequence, position) in chrPos2alleleBeta:
            allele = chrPos2alleleBeta[(sequence, position)][0]
        else:
            continue
        # fill genotypes #
        individualGenotypes = individualGenotypes.split('\t')
        for i, variantColumn in enumerate(variantColumns):
            variant = individualGenotypes[variantColumn]
            if (variant != 'NA') and abundances[abundanceOffset + variantColumn]:
                mask[i] = True
                x[i, 1] = float(variant.split('/')[allele])
            else:
                mask[i] = False
        # mwas #
        phenotypes_ = phenotypes[mask]
        if (flag == 0 and numpy.count_nonzero(phenotypes_ == 0) >= m0 and numpy.count_nonzero(phenotypes_ == 1) >= m1) or (flag == 1 and numpy.count_nonzero(mask) >= m):
            x[mask, 1] -= numpy.mean(x[mask, 1])
            if numpy.any(x[mask, 1]):
                glm = Glm(phenotypes_, x[mask])
                try:
                    glmResult = glm.fit(maxiter = 9999)
                    openFile.write(f'{sequence}\t{position}\t{mapping1[allele]}\t{glmResult.params[1]}\t{glmResult.bse[1]}\t{glmResult.pvalues[1]}\n')
                except Exception:
                    openFile.write(f'{sequence}\t{position}\t{mapping1[allele]}\t0\t0\t1\n')
            else:
                openFile.write(f'{sequence}\t{position}\t{mapping1[allele]}\t0\t0\t1\n')
            N += 1
    openFile.close()
    n.acquire()
    n.value += N
    n.release()
    return None


def createTestProcesses(threads, phenotypes, covariates, abundances, variantColumns, m, chrPos2alleleBeta):
    queue = Queue(threads)
    processes = list()
    files = list()
    n = Value(c_int64, 0)
    with threadpool_limits(limits = 1):
        for i in range(threads):
            files.append(uuid4().hex)
            processes.append(
                Process(
                    target = testProcess,
                    args = (queue, files[-1], n, m, phenotypes, covariates, abundances, variantColumns, chrPos2alleleBeta)
                )
            )
            processes[-1].start()
    return (queue, processes, files, n)


def gdmProcess(queue, lock, x, y, n):
    x = numpy.ndarray(shape = (n, n), dtype = numpy.float32, buffer = x)
    y = numpy.ndarray(shape = (n, n), dtype = numpy.float32, buffer = y)
    matrix4 = numpy.empty(shape = 4, dtype = numpy.float32)
    matrixN1 = numpy.empty(shape = n, dtype = numpy.bool_)
    matrixN2 = numpy.empty(shape = n, dtype = numpy.float32)
    matrixNN1 = numpy.empty(shape = (n, n), dtype = numpy.bool_)
    matrixNN2 = numpy.empty(shape = (n, n), dtype = numpy.float32)
    X = numpy.zeros(shape = (n, n), dtype = numpy.float32)
    Y = numpy.zeros(shape = (n, n), dtype = numpy.float32)
    while True:
        populationGenotype, individualGenotypes = queue.get()
        if populationGenotype is None:
            break
        matrix4[ : ] = populationGenotype.split('/')
        majorAllele = numpy.argmax(matrix4)
        for i, individual in enumerate(individualGenotypes.split('\t')):
            if individual == 'NA':
                matrixN1[i] = False
            else:
                matrixN1[i] = True
                matrixN2[i] = float(individual.split('/')[majorAllele])
        numpy.bitwise_and(matrixN1, matrixN1[ : , None], out = matrixNN1)
        numpy.subtract(matrixN2, matrixN2[ : , None], out = matrixNN2, where = matrixNN1)
        numpy.square(matrixNN2, out = matrixNN2, where = matrixNN1)
        numpy.add(X, matrixNN2, out = X, where = matrixNN1)
        Y += matrixNN1
    lock.acquire()
    x += X
    y += Y
    lock.release()
    return None


def createGDMProcesses(threads, n):
    x = RawArray(c_double, n * n)
    X = numpy.ndarray(shape = (n, n), dtype = numpy.float32, buffer = x)
    y = RawArray(c_double, n * n)
    Y = numpy.ndarray(shape = (n, n), dtype = numpy.float32, buffer = y)
    Y += 10 * numpy.finfo(numpy.float32).eps
    queue = Queue(threads)
    lock = Lock()
    processes = list()
    with threadpool_limits(limits = 1):
        for i in range(threads):
            processes.append(Process(target = gdmProcess, args = (queue, lock, x, y, n)))
            processes[-1].start()
    return (queue, processes, lock, X, Y)


def freeProcesses(queue, processes, n):
    none = (None, ) * n
    for process in processes:
        queue.put(none)
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def readVariantFileHeader(file: str) -> tuple[Sequence[str], Mapping[str, int], NDArray[numpy.int64]]:
    id2index = dict()
    x = list()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    line = openFile.readline().rstrip('\n')
    assert line.startswith('Chromosome\tPosition\tPopulation'), f'\"{file}\" is not a valid variant file.'
    individuals = line.split('\t')[3 : ]
    for i, j in enumerate(individuals):
        id2index[j] = i
        x.append(i)
    openFile.close()
    return (individuals, id2index, numpy.asarray(x, dtype = numpy.int64))


def readVariantFile(file):
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    line = openFile.readline().rstrip('\n')
    assert line.startswith('Chromosome\tPosition\tPopulation'), f'\"{file}\" is not a valid variant file.'
    for line in openFile:
        yield line.rstrip('\n').split('\t', maxsplit = 3)
    openFile.close()
    return None


def normalizeX(x):
    x = x - numpy.mean(x, axis = 0)
    y = numpy.einsum('ij,ij->j', x, x)
    y /= x.shape[0] - 1
    y **= 0.5
    y += 10 * numpy.finfo(numpy.float32).eps
    x /= y
    return x


def bonferroniCorrection(alpha, classifications, sequence2classification, inputFiles, outputFile):
    open4w = open(outputFile, 'w')
    open4w.write('Classification\tChromosome\tPosition\tAllele\tBeta\tSE\tP\tSignificance\n')
    processBar = ProcessBar(len(inputFiles))
    for i, inputFile in enumerate(inputFiles, start = 1):
        open4r = open(inputFile, 'r')
        # Chromosome Position Allele Beta SE P ... #
        for line in open4r:
            lines = line.rstrip('\n').split('\t')
            lines.append('+' if float(lines[5]) <= alpha else '-')
            lines.insert(0, classifications[sequence2classification[lines[0]]])
            open4w.write('\t'.join(lines) + '\n')
        open4r.close()
        os.remove(inputFile)
        processBar.plot(i)
    open4w.close()
    return None


def write2file(alpha, classifications, sequence2classification, inputFiles, outputFile, chrPos2alleleBeta):
    open4w = open(outputFile, 'w')
    open4w.write('Classification\tChromosome\tPosition\tAllele\tBeta\tSE\tP\tSignificance\tDirection\n')
    for inputFile in inputFiles:
        open4r = open(inputFile, 'r')
        # Chromosome Position Allele Beta SE P #
        for line in open4r:
            lines = line.rstrip('\n').split('\t')
            Beta = chrPos2alleleBeta.get((lines[0], lines[1]), (None, 0))[1]
            beta = float(lines[3])
            if not Beta:
                direction = '.'
            elif (beta > 0 and Beta > 0) or (beta < 0 and Beta < 0):
                direction = '+'
            else:
                direction = '-'
            lines.append('+' if float(lines[5]) <= alpha else '-')
            lines.append(direction)
            lines.insert(0, classifications[sequence2classification[lines[0]]])
            open4w.write('\t'.join(lines) + '\n')
        open4r.close()
        os.remove(inputFile)
    open4w.close()
    return None


def main(parameters):
    individuals, id2index, variantColumns = readVariantFileHeader(parameters.variant)
    nIndividuals = len(individuals)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading phenotype file.', flush = True)
    flag, phenotypes = readPhenotypeFile(parameters.phenotype, parameters.phenotype_column, id2index)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading annotation file.', flush = True)
    classifications, cluster2classification = readAnnotationFile(parameters.annotation)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading mapping file.', flush = True)
    sequence2classification = readMappingFile(parameters.mapping, cluster2classification)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading abundance file.', flush = True)
    abundances = readAbundanceFile(parameters.abundance, classifications, individuals)

    if parameters.covariate is not None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading covariate file.', flush = True)
        flag_, covariates = readCovariateFile(parameters.covariate, id2index)
        flag &= flag_

    if parameters.mwas is not None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading mwas file.', flush = True)
        chrPos2alleleBeta = readMWASFile(parameters.mwas)
    else:
        chrPos2alleleBeta = None

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Creating gdm.', flush = True)
    queue, processes, lock, gdm, normalizer = createGDMProcesses(parameters.threads, len(id2index))
    N = 0
    for i in readVariantFile(parameters.variant):
        # i: [Chromosome, Position, Population, ID<tab>ID] #
        queue.put((i[2], i[3]))
        N += 1
        if not N % 10000:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> {N} variants have been loaded.', end = '\r', flush = True)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> {N} variants have been loaded.', flush = True)
    freeProcesses(queue, processes, 2)
    gdm /= normalizer
    gdm **= 0.5

    indices = numpy.flatnonzero(flag)
    individuals = [individuals[i] for i in indices]
    phenotypes = phenotypes[indices]
    variantColumns = variantColumns[indices]
    pca = PCA(n_components = min(parameters.pca_components, len(individuals) - 1), whiten = False, random_state = 0)
    if parameters.covariate is not None:
        covariates = numpy.hstack(
            (
                normalizeX(covariates[indices]),
                pca.fit_transform(gdm[indices][ : , indices])
            ),
            dtype = numpy.float32
        )
    else:
        covariates = pca.fit_transform(gdm[indices][ : , indices])

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Associating phenotypes with variants.', flush = True)
    queue, processes, tempfiles, n = createTestProcesses(parameters.threads, phenotypes, covariates, abundances, variantColumns, parameters.size, chrPos2alleleBeta)
    processBar = ProcessBar(N)
    for i, j in enumerate(readVariantFile(parameters.variant), start = 1):
        # i: [Chromosome, Position, Population, ID<tab>ID] #
        queue.put((sequence2classification[j[0]] * nIndividuals, j[0], j[1], j[2], j[3]))
        if not i % 1000 or i == N:
            processBar.plot(i)
    freeProcesses(queue, processes, 5)
    if n.value:
        if parameters.mwas:
            write2file(parameters.alpha, classifications, sequence2classification, tempfiles, parameters.output, chrPos2alleleBeta)
        else:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Performing the Bonferroni correction (alpha = {parameters.alpha}, n = {n.value}).', flush = True)
            bonferroniCorrection(parameters.alpha / n.value, classifications, sequence2classification, tempfiles, parameters.output)
    else:
        for i in tempfiles:
            os.remove(i)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
