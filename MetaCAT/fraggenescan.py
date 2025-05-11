import os
from ctypes import c_int64
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value
from subprocess import DEVNULL, run
from uuid import uuid4

from .fasta import split_fasta as splitFasta
from .processbar import ProcessBar


def fraggenescanProcess(queue, fraggenescan, type, n, processBar):
    while True:
        file = queue.get()
        if file is None:
            break
        completedProcess = run(
            [fraggenescan, '-w', '0', '-t', 'complete', '-s', file, '-o', file],
            stdout = DEVNULL, stderr = DEVNULL
        )
        assert not completedProcess.returncode, 'An error has occured while running FragGeneScan.'
        if os.access(file, os.F_OK):
            os.remove(file)
        if os.access(file + '.gff', os.F_OK):
            os.remove(file + '.gff')
        if os.access(file + '.out', os.F_OK):
            os.remove(file + '.out')
        if type == 'protein':
            os.rename(file + '.faa', file)
            if os.access(file + '.ffn', os.F_OK):
                os.remove(file + '.ffn')
        elif type == 'cds':
            os.rename(file + '.ffn', file)
            if os.access(file + '.faa', os.F_OK):
                os.remove(file + '.faa')
        n.acquire()
        n.value += 1
        processBar.plot(n.value)
        n.release()
    return None


def createProcesses(fraggenescan, type, n, threads):
    queue = Queue(threads)
    processes = list()
    N = Value(c_int64, 0)
    processBar = ProcessBar(n)
    for i in range(threads):
        processes.append(Process(target = fraggenescanProcess, args = (queue, fraggenescan, type, N, processBar)))
        processes[-1].start()
    return (queue, processes, N)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put(None)
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def runFraggenescan(fraggenescan, file, threads, type = 'protein'):
    '''
    type: protein | cds
    '''
    n = 3 * threads
    tempfiles = list()
    queue, processes, N = createProcesses(fraggenescan, type, n, threads)
    for tempfile in splitFasta(file, n):
        tempfiles.append(tempfile)
        queue.put(tempfile)
    freeProcesses(queue, processes)
    processBar = ProcessBar(1)
    processBar.plot(1)
    FILE = uuid4().hex
    open4w = open(FILE, 'wb')
    for tempfile in tempfiles:
        open4r = open(tempfile, 'rb')
        while open4w.write(open4r.read(10485760)):
            pass
        open4r.close()
        os.remove(tempfile)
    open4w.close()
    return FILE
