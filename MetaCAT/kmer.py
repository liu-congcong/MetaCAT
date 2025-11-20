from ctypes import c_float, c_int64
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import RawArray, Value

import numpy

from .processbar import ProcessBar


def kmer_to_index(k):
    acgt = ('A', 'C', 'G', 'T')
    tgca = ('T', 'G', 'C', 'A')
    index = 0
    kmer2index = dict()
    for i in range(4 ** k):
        kmer_acgt = ''.join(acgt[i >> 2 * j & 3] for j in range(k))
        kmer_tgca = ''.join(tgca[i >> 2 * (k - 1 - j) & 3] for j in range(k))
        if kmer_acgt not in kmer2index:
            kmer2index[kmer_acgt] = index
            kmer2index[kmer_tgca] = index
            index += 1
    return kmer2index


def count_kmers_thread(processQueue, x, n, N):
    processBar = ProcessBar(N)
    kmer2index = kmer_to_index(4)
    while True:
        sequence, offset = processQueue.get()
        if sequence is None:
            break
        for i in range(len(sequence) - 3):
            kmer = sequence[i : i + 4]
            if kmer in kmer2index:
                j = offset + kmer2index[kmer]
                x[j] += 1
        n.acquire()
        n.value += 1
        processBar.plot(n.value)
        n.release()
    return None


def createProcesses(n, threads):
    x = RawArray(c_float, 136 * n)
    processes = list()
    queue = Queue(threads)
    N = Value(c_int64, 0)
    for i in range(threads):
        processes.append(Process(target = count_kmers_thread, args = (queue, x, N, n)))
        processes[-1].start()
    return (queue, processes, x, N)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def count_kmers(n, sequenceStructs, threads):
    queue, processes, sharedX, N = createProcesses(n, threads)
    # put all sequences to queue #
    for i, sequenceStruct in enumerate(sequenceStructs):
        queue.put((sequenceStruct[1], 136 * i))
    freeProcesses(queue, processes)
    # normalize kmers #
    x = numpy.ndarray(shape = (n, 136), dtype = numpy.float32, buffer = sharedX)
    x += 1 / 136
    x /= numpy.sum(x, axis = 1, keepdims = True)
    return (x, sharedX)
