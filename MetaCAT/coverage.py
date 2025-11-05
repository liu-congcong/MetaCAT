import gzip
from ctypes import c_float, c_int64
from datetime import datetime
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import RawArray, Value

import numpy
from threadpoolctl import threadpool_limits

from .bam import getUngappedRegions, indexBam, readIndices
from .processbar import ProcessBar


def createProcesses(files, n, trim, mapq, identity, threads):
    m = len(files)
    x = RawArray(c_float, n * m * 2)
    X = numpy.ndarray(shape = (n, m, 2), dtype = numpy.float32, buffer = x)
    N = Value(c_int64, 0)
    queue = Queue(threads)
    processes = list()
    for i in range(threads):
        processes.append(Process(target = countProcess, args = (queue, x, n, m, files, trim, mapq, identity, N)))
        processes[-1].start()
    return (queue, processes, X, N)


@threadpool_limits.wrap(limits = 1)
def countProcess(queue, x, n, m, files, trim, mapq, identity, N):
    '''
    queue: queue
    x: (#sequences, #files, 2)
    n: #sequences
    m: #files
    trim: int
    mapq: int
    identity: float
    N: Value(c_int64, 0)
    '''
    processBar = ProcessBar(n)
    x = numpy.ndarray(shape = (n, m, 2), dtype = numpy.float32, buffer = x)
    while True:
        sequence, length, fileOffsets, dataOffsets, dataSizes = queue.get()
        if sequence is None:
            break
        y = numpy.zeros(shape = (length, m), dtype = numpy.int32)
        for i, (file, fileOffset, dataOffset, dataSize) in enumerate(zip(files, fileOffsets, dataOffsets, dataSizes)):
            if dataSize:
                for j in getUngappedRegions(file, fileOffset, dataOffset, dataSize, mapq, identity):
                    for region in j[2]:
                        y[region[0] : region[1], i] += 1
        numpy.mean(y[trim : length - trim], axis = 0, out = x[sequence, : , 0])
        numpy.var(y[trim : length - trim], axis = 0, out = x[sequence, : , 1], ddof = 1)
        N.acquire()
        N.value += 1
        processBar.plot(N.value)
        N.release()
    return None


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None, None, None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def writeFile(x, sequences, lengths, minLength, output):
    openFile = gzip.open(output, mode = 'wt', compresslevel = 9)
    openFile.write('\t'.join(['Sequence ID'] + (['Mean', 'Variance'] * x.shape[1])) + '\n')
    for sequence, length, y in zip(sequences, lengths, x):
        if length >= minLength:
            y = '\t'.join(str(i) for i in y.flatten())
            openFile.write(f'{sequence}\t{y}\n')
    openFile.close()
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Indexing all bam files.', flush = True)
    indexBam(parameters.bam, parameters.threads_index)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading all index files.', flush = True)
    sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndices(parameters.bam, parameters.threads_index)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Computing coverage.', flush = True)
    queue, processes, x, N = createProcesses(parameters.bam, len(sequences), parameters.trim, parameters.mapq, parameters.identity, parameters.threads_count)
    for i, (length, fileOffset, dataOffset, dataSize) in enumerate(zip(lengths, fileOffsets, dataOffsets, dataSizes)):
        if length >= parameters.min_sequence_length:
            queue.put((i, length, fileOffset, dataOffset, dataSize))
    freeProcesses(queue, processes)

    # write output #
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Writing to \"{parameters.output}\".', flush = True)
    writeFile(x, sequences, lengths, parameters.min_sequence_length, parameters.output)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
