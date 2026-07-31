from datetime import datetime
from pathlib import Path
from shutil import rmtree
from subprocess import PIPE, STDOUT, run
from sys import exit, stderr
from uuid import uuid4

from .checkm2 import readCheckm2File


class Gtdbtk:
    def __init__(self, gtdbtk, threads, clusters, temp):
        self.gtdbtk = gtdbtk
        self.threads = threads
        self.clusters = clusters
        self.temp = Path(temp)
        return None

    def parseInputFiles(self, files):
        temp = self.temp.joinpath(uuid4().hex)
        openFile = temp.open('w')
        for i in files:
            cluster = Path(i).stem
            if self.clusters is None:
                openFile.write(f'{i}\t{cluster}\n')
            elif cluster in self.clusters:
                openFile.write(f'{i}\t{cluster}\n')
        openFile.close()
        return temp

    def classifyWf(self, inputFiles, outputFile):
        temp1 = self.parseInputFiles(inputFiles)
        temp2 = self.temp.joinpath(uuid4().hex)
        command = [self.gtdbtk, 'classify_wf', '--batchfile', str(temp1), '--out_dir', str(temp2), '--cpus', str(self.threads), '--force']
        completedProcess = run(command, stdout = None, stderr = None)
        if completedProcess.returncode:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to run GTDB-Tk.', file = stderr, flush = True)
            exit(1)
        temp1.unlink()

        flag = True
        open4w = open(outputFile, 'w')
        for i in ('bac120', 'ar53'):
            path = temp2.joinpath(f'gtdbtk.{i}.summary.tsv')
            if path.is_file():
                open4r = path.open('r')
                header = open4r.readline().rstrip('\n')
                if flag:
                    open4w.write(header + '\n')
                    flag = False
                for line in open4r:
                    open4w.write(line.rstrip('\n') + '\n')
                open4r.close()
        open4w.close()
        rmtree(str(temp2))
        return self


def main(parameters):
    if parameters.checkm2:
        clusters = readCheckm2File(parameters.checkm2, parameters.contamination, parameters.completeness)
    else:
        clusters = None
    gtdbtk = Gtdbtk(parameters.gtdbtk, parameters.threads, clusters, parameters.temp)
    gtdbtk.classifyWf(parameters.fasta, parameters.output)
    return None
