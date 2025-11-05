import gzip
import os
from collections.abc import Mapping, Sequence
from ctypes import c_double, c_int64
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


def readPhenotypeFile(file, phenotypeColumn, individual2i):
    phenotypeColumn = phenotypeColumn - 1
    n = len(individual2i)
    flag = numpy.zeros(shape = n, dtype = numpy.bool_)
    phenotypes = numpy.empty(shape = n, dtype = numpy.float32)
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().startswith('ID'), f'\"{file}\" is not a valid phenotype file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if (lines[0] in individual2i) and (lines[phenotypeColumn].lower() != 'na'):
            flag[individual2i[lines[0]]] = True
            phenotypes[individual2i[lines[0]]] = float(lines[phenotypeColumn])
    openFile.close()
    return (flag, phenotypes)


def readCovariateFile(file, individual2i):
    n = len(individual2i)
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
        if (id in individual2i) and ('na' not in covariates.lower()):
            flag[individual2i[id]] = True
            x[individual2i[id]] = covariates.split('\t')
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
    chrPos2allele = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    assert openFile.readline().startswith('Classification\tChromosome\tPosition\tAllele\tBeta\tSE\tP\tSignificance'), f'\"{file}\" is not a valid mwas file.'
    for line in openFile:
        lines = line.rstrip('\n').split('\t')
        if lines[7] == '+':
            chrPos2allele[(lines[1], lines[2])] = mapping2[lines[3]]
    openFile.close()
    return chrPos2allele


@threadpool_limits.wrap(limits = 1)
def testProcess(queue, file, m, phenotypes, covariates, variantColumns, classification2abundance, chrPos2allele, n):
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
    nTests = 0
    x = numpy.hstack((numpy.ones(shape = (phenotypes.size, 3), dtype = numpy.float64), covariates))
    mask = numpy.empty(shape = phenotypes.size, dtype = numpy.bool_)
    acgt = numpy.empty(shape = 4, dtype = numpy.float64)
    openFile = open(file, 'w')
    # Classification, Chromosome, Position, Allele, Beta, SE, P #
    numpy.seterr(all = 'raise')
    while True:
        i = queue.get()
        if i is None:
            break
        I = i.split('\t')
        if chrPos2allele is None:
            acgt[ : ] = I[3].split('/')
            allele = numpy.argmax(acgt)
        elif (I[1], I[2]) in chrPos2allele:
            allele = chrPos2allele[(I[1], I[2])]
        else:
            continue
        # fill genotypes #
        for i, j in enumerate(zip(classification2abundance[I[0]], variantColumns)):
            if j[0] > 0 and I[j[1]] != 'NA':
                mask[i] = True
                x[i, 1] = float(I[j[1]].split('/')[allele])
                x[i, 2] = j[0]
            else:
                mask[i] = False
        # mwas #
        phenotypes_ = phenotypes[mask]
        if (flag == 0 and numpy.count_nonzero(phenotypes_ == 0) >= m0 and numpy.count_nonzero(phenotypes_ == 1) >= m1) or (flag == 1 and numpy.count_nonzero(mask) >= m):
            nTests += 1
            glm = Glm(phenotypes_, x[mask])
            try:
                glmResult = glm.fit(maxiter = 9999)
                openFile.write(f'{I[0]}\t{I[1]}\t{I[2]}\t{mapping1[allele]}\t{glmResult.params[1]}\t{glmResult.bse[1]}\t{glmResult.pvalues[1]}\n')
            except Exception:
                openFile.write(f'{I[0]}\t{I[1]}\t{I[2]}\t{mapping1[allele]}\t0\t0\t1\n')
    openFile.close()
    n.acquire()
    n.value += nTests
    n.release()
    return None


def createTestProcesses(threads, m, phenotypes, covariates, variantColumns, classification2abundance, chrPos2allele):
    queue = Queue(threads)
    processes = list()
    files = list()
    n = Value(c_int64, 0)
    for i in range(threads):
        files.append(uuid4().hex)
        processes.append(
            Process(
                target = testProcess,
                args = (queue, files[-1], m, phenotypes, covariates, variantColumns, classification2abundance, chrPos2allele, n)
            )
        )
        processes[-1].start()
    return (queue, processes, files, n)


@threadpool_limits.wrap(limits = 1)
def gdmProcess(queue, lock, x, y, n, classification2abundance):
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
        i = queue.get()
        if i is None:
            break
        I = i.split('\t')
        matrix4[ : ] = I[3].split('/')
        majorAllele = numpy.argmax(matrix4)
        for i, j in enumerate(zip(classification2abundance[I[0]], I[4 : ])):
            if j[0] > 0 and j[1] != 'NA':
                matrixN1[i] = True
                matrixN2[i] = float(j[1].split('/')[majorAllele])
            else:
                matrixN1[i] = False
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


def createGDMProcesses(threads, n, classification2abundance):
    x = RawArray(c_double, n * n)
    X = numpy.ndarray(shape = (n, n), dtype = numpy.float32, buffer = x)
    y = RawArray(c_double, n * n)
    Y = numpy.ndarray(shape = (n, n), dtype = numpy.float32, buffer = y)
    Y += 10 * numpy.finfo(numpy.float32).eps
    queue = Queue(threads)
    lock = Lock()
    processes = list()
    for i in range(threads):
        processes.append(Process(target = gdmProcess, args = (queue, lock, x, y, n, classification2abundance)))
        processes[-1].start()
    return (queue, processes, lock, X, Y)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put(None)
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def readVariantFileHeader(file: str) -> tuple[Sequence[str], Mapping[str, int], NDArray[numpy.int64]]:
    individual2i = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    line = openFile.readline().rstrip('\n')
    assert line.startswith('Classification\tChromosome\tPosition\tPopulation'), f'\"{file}\" is not a valid variant file.'
    individuals = line.split('\t')[4 : ]
    for i, j in enumerate(individuals):
        individual2i[j] = i
    openFile.close()
    return (individuals, individual2i)


def readVariantFile(file):
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    line = openFile.readline().rstrip('\n')
    assert line.startswith('Classification\tChromosome\tPosition\tPopulation'), f'\"{file}\" is not a valid variant file.'
    for line in openFile:
        yield line.rstrip('\n')
    openFile.close()
    return None


def readAbundanceFile(file, individuals):
    x = list()
    classification2i = dict()
    if isGzipped(file):
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, 'r')
    lines = openFile.readline().rstrip('\n').split('\t')
    individual2i = dict((j, i) for i, j in enumerate(lines[1 : ]))
    assert lines[0] == 'Abundance', f'\"{file}\" is not a valid abundance file.'
    y = [individual2i[i] for i in individuals]
    for i, line in enumerate(openFile):
        lines = line.rstrip('\n').split('\t')
        classification2i[lines[0]] = i
        x.append(lines[1 : ])
    openFile.close()
    x = numpy.asarray(x, dtype = numpy.float64)[ : , y]
    return (x, classification2i)


def normalizeX(x):
    x -= numpy.mean(x, axis = 0)
    y = numpy.einsum('ij,ij->j', x, x)
    y /= x.shape[0] - 1
    y **= 0.5
    y += 10 * numpy.finfo(numpy.float64).eps
    x /= y
    return None


def writeFile(inputFiles, outputFile, alpha):
    processBar = ProcessBar(len(inputFiles))
    open4w = open(outputFile, 'w', buffering = 10485760)
    open4w.write('Classification\tChromosome\tPosition\tAllele\tBeta\tSE\tP\tSignificance\n')
    for i, inputFile in enumerate(inputFiles, start = 1):
        open4r = open(inputFile, 'r')
        for line in open4r:
            line = line.rstrip('\n')
            significance = '-' if float(line.rsplit('\t', maxsplit = 1)[1]) > alpha else '+'
            open4w.write(f'{line}\t{significance}\n')
        open4r.close()
        os.remove(inputFile)
        processBar.plot(i)
    open4w.close()
    return None


def main(parameters):
    individuals, individual2i = readVariantFileHeader(parameters.variant)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading phenotype file.', flush = True)
    flag, phenotypes = readPhenotypeFile(parameters.phenotype, parameters.phenotype_column, individual2i)

    if parameters.covariate is not None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading covariate file.', flush = True)
        flag_, covariates = readCovariateFile(parameters.covariate, individual2i)
        flag &= flag_

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading abundance file.', flush = True)
    abundance, classification2i = readAbundanceFile(parameters.abundance, individuals)
    classification2abundance = dict()
    for classification, i in classification2i.items():
        classification2abundance[classification] = abundance[i]

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Creating gdm.', flush = True)
    queue, processes, lock, gdm, normalizer = createGDMProcesses(parameters.threads, len(individuals), classification2abundance)
    n = 0
    for i in readVariantFile(parameters.variant):
        # i: Classification, Chromosome, Position, Population, ID, ... #
        queue.put(i)
        n += 1
        if not n % 10000:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> {n} variants have been loaded.', end = '\r', flush = True)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> {n} variants have been loaded.', flush = True)
    freeProcesses(queue, processes)
    gdm /= normalizer
    gdm **= 0.5
    numpy.putmask(gdm, gdm == 0, 1)
    numpy.savez(parameters.output + '.gdm.npz', header = numpy.array(individuals), gdm = gdm, n = n)

    if parameters.mwas is not None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading mwas file.', flush = True)
        chrPos2allele = readMWASFile(parameters.mwas)
    else:
        chrPos2allele = None

    indices = numpy.flatnonzero(flag)
    individuals = [individuals[i] for i in indices]
    phenotypes = phenotypes[indices]
    abundance = abundance[ : , indices]

    classification2abundance.clear()
    for classification, i in classification2i.items():
        classification2abundance[classification] = abundance[i]

    pca = PCA(n_components = min(parameters.pca_components, len(individuals) - 1), whiten = False, random_state = 0)
    if parameters.covariate is not None:
        covariates = numpy.hstack(
            (
                covariates[indices],
                pca.fit_transform(gdm[indices][ : , indices])
            ),
            dtype = numpy.float32
        )
    else:
        covariates = pca.fit_transform(gdm[indices][ : , indices])
    normalizeX(covariates)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Associating phenotypes with variants.', flush = True)
    queue, processes, tempfiles, nTests = createTestProcesses(parameters.threads, parameters.size, phenotypes, covariates, indices + 4, classification2abundance, chrPos2allele)
    processBar = ProcessBar(n)
    for i, j in enumerate(readVariantFile(parameters.variant), start = 1):
        queue.put(j)
        if not i % 1000 or i == n:
            processBar.plot(i)
    freeProcesses(queue, processes)
    if parameters.mwas is None:
        parameters.alpha /= max(1, nTests.value)
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Bonferroni-corrected for {nTests.value} tests.', flush = True)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing to \"{parameters.output}\".', flush = True)
    writeFile(tempfiles, parameters.output, parameters.alpha)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
