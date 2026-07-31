from ctypes import c_int64
from datetime import datetime
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value
from pathlib import Path
from subprocess import DEVNULL, run
from sys import exit, stderr

from .fasta import splitFastaFile
from .processbar import ProcessBar


class Fraggenescan:
    def __init__(self, fraggenescan, threads, type, temp):
        self.fraggenescan = fraggenescan
        self.threads = threads
        self.type = type
        self.temp = temp
        return None

    def workerProcess(self, queue, n):
        processBar = ProcessBar(3 * self.threads)
        while True:
            file = queue.get()
            if file is None:
                break
            completedProcess = run(
                [self.fraggenescan, '-w', '0', '-t', 'complete', '-s', file, '-o', file],
                stdout = DEVNULL, stderr = None
            )
            if completedProcess.returncode:
                print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to run FragGeneScan.', file = stderr, flush = True)
                exit(1)
            path = Path(file)
            path.unlink(missing_ok = True)
            path.with_suffix('.gff').unlink(missing_ok = True)
            path.with_suffix('.out').unlink(missing_ok = True)
            if self.type == 'protein':
                path.with_suffix('.faa').rename(file)
                path.with_suffix('.ffn').unlink(missing_ok = True)
            elif self.type == 'cds':
                path.with_suffix('.ffn').rename(file)
                path.with_suffix('.faa').unlink(missing_ok = True)
            with n.get_lock():
                n.value += 1
                processBar.plot(n.value)
        return None

    def createProcesses(self):
        queue = Queue(self.threads)
        processes = list()
        n = Value(c_int64, 0)
        for i in range(self.threads):
            processes.append(Process(target = self.workerProcess, args = (queue, n)))
            processes[-1].start()
        return (queue, processes, n)

    def freeProcesses(self, queue, processes):
        for process in processes:
            queue.put(None)
        queue.close()
        queue.join_thread()
        for process in processes:
            process.join()
            process.close()
        return self

    def main(self, inputFile, outputFile):
        tempfiles = list()
        queue, processes, n = self.createProcesses()
        for tempfile in splitFastaFile(inputFile, 3 * self.threads, self.temp):
            tempfiles.append(tempfile)
            queue.put(tempfile)
        self.freeProcesses(queue, processes)
        processBar = ProcessBar(1)
        processBar.plot(1)
        open4w = open(outputFile, 'wb', buffering = 10485760)
        for tempfile in tempfiles:
            open4r = open(tempfile, 'rb', buffering = 10485760)
            while open4w.write(open4r.read(10485760)):
                pass
            open4r.close()
            Path(tempfile).unlink()
        open4w.close()
        return self
