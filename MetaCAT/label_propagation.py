from ctypes import c_int32
from multiprocessing import JoinableQueue, Process
from multiprocessing.sharedctypes import RawArray

import numpy
from scipy.sparse import csr_array
from threadpoolctl import threadpool_limits


def labelPropagation_standard(x, y, z, maxIterations = 100, tolerance = 1e-3):
    '''
    x: sparse affinity matrix (n, n)
    y: label vector (n, )
    z: label matrix (n, k)
    '''
    z[ : , : ] = 0
    uniqueY = numpy.unique(y)[1 : ]
    labeledY = y != -1
    for i in uniqueY:
        z[y == i, uniqueY == i] = 1
    labeledZ = z[labeledY]
    for i in range(maxIterations):
        tempZ = x @ z # (n, n) (n, k) -> (n, k) #
        tempZ /= numpy.sum(tempZ, axis = 1, keepdims = True) + 10 * numpy.finfo(x.dtype).eps
        tempZ[labeledY] = labeledZ
        z -= tempZ
        numpy.abs(z, out = z)
        if numpy.sum(z) < tolerance:
            z = tempZ
            break
        z = tempZ
    numpy.argmax(z, axis = 1, out = y)
    numpy.putmask(y, z[numpy.arange(y.size), y] == 0, -1)
    return None


def labelPropagation(x, y, z, maxIterations = 100, minP = 0, tolerance = 1e-3):
    '''
    x: sparse affinity matrix (n, n)
    y: label vector (n, )
    z: label matrix (n, k)
    '''
    z[ : , : ] = 0
    uniqueY = numpy.unique(y)[1 : ]
    labeledY = y != -1
    for i in uniqueY:
        z[y == i, uniqueY == i] = 1
    labeledZ = z[labeledY]
    for i in range(maxIterations):
        tempZ = x @ z # (n, n) (n, k) -> (n, k) #
        tempZ /= numpy.sum(tempZ, axis = 1, keepdims = True) + 10 * numpy.finfo(x.dtype).eps
        tempZ[labeledY] = labeledZ
        z -= tempZ
        numpy.abs(z, out = z)
        if numpy.sum(z) < tolerance:
            z = tempZ
            break
        z = tempZ
    numpy.argmax(z, axis = 1, out = y)
    p = z[numpy.arange(y.size), y]
    numpy.putmask(y, p <= minP, -1)
    return None


def annealingLabelPropagation(x, y, z = ((0.97, 10, 0), (0.9, 10, 0), (0.7, 30, 0))):
    '''
    x: sparse affinity matrix (n, n)
    y: label vector (n, ), >= -1
    z:
    '''
    uniqueY = numpy.unique(y)[1 : ]
    r = numpy.empty((x.shape[0], uniqueY.size), dtype = x.dtype) # (n, k) #
    for i, j, k in z:
        X = csr_array(
            (numpy.where(x.data >= i, x.data, 0), x.indices, x.indptr),
            shape = x.shape,
            dtype = numpy.float32
        )
        labelPropagation(X, y, r, j, k)
        if not numpy.any(y == -1):
            break
    Y = uniqueY[y]
    numpy.putmask(Y, y == -1, -1)
    return Y


@threadpool_limits.wrap(limits = 1)
def workerProcess(processQueue, partitions, n, m):
    partitions = numpy.ndarray(shape = (n, m), dtype = numpy.int32, buffer = partitions)
    Y = numpy.empty(shape = m, dtype = numpy.int32)
    while True:
        i, x, y = processQueue.get()
        if i is None:
            processQueue.task_done()
            break
        Y[ : x.shape[0]] = -1
        Y[y] = y
        annealingLabelPropagation(x, Y[ : x.shape[0]])
        partitions[i, : x.shape[0]] = Y[ : x.shape[0]]
        processQueue.task_done()
    return None


def createProcesses(n, m, t):
    '''
    m: number of total seed sequences
    '''
    partitions = RawArray(c_int32, n * m)
    processQueue = JoinableQueue(t)
    processes = list()
    for i in range(t):
        processes.append(Process(target = workerProcess, args = (processQueue, partitions, n, m)))
        processes[-1].start()
    return (processQueue, partitions, processes)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None, None))
    queue.close()
    queue.join()
    for process in processes:
        process.join()
        process.close()
    return None


def generateRandomData(n, k, random = 0):
    randomGenerator = numpy.random.default_rng(random)
    means = randomGenerator.random(size = (k, 2), dtype = numpy.float32) + numpy.arange(k)[ : , None] * 1
    covariance = numpy.identity(2, dtype = numpy.float32) * 0.1
    x = numpy.empty(shape = (k * (n + 1), 2), dtype = numpy.float32)
    y = numpy.full(shape = k * (n + 1), fill_value = -1, dtype = numpy.int32)
    for i in range(k):
        x[i * n : (i + 1) * n] = randomGenerator.multivariate_normal(
            mean = means[i], cov = covariance, size = n
        )
        x[n * k + i] = means[i]
        y[n * k + i] = i
    return (x, y)


def convert2csr(x):
    data = numpy.empty(shape = x.shape, dtype = numpy.float32)
    indices = numpy.empty(shape = x.shape, dtype = numpy.int32)
    indptr = numpy.arange(x.shape[0] + 1, dtype = numpy.int32) * (x.shape[0])
    for i, j in enumerate(x):
        k = numpy.argsort(j)
        data[i] = 1 - x[i][k]
        indices[i] = k
    y = csr_array((data.flatten(), indices.flatten(), indptr), shape = data.shape, dtype = numpy.float32)
    return y
