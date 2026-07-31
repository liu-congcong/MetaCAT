from datetime import datetime
from pathlib import Path
from subprocess import DEVNULL, run
from sys import exit, stderr
from uuid import uuid4


class Mash:
    def __init__(self, mash, d, p, k, s, temp):
        self.mash = mash
        self.d = d
        self.p = p
        self.k = k
        self.s = s
        self.temp = Path(temp)
        return None

    def parseInputFiles(self, files):
        temp = self.temp.joinpath(uuid4().hex)
        openFile = temp.open('w')
        for i in files:
            openFile.write(i + '\n')
        openFile.close()
        return temp

    def sketch(self, inputFile):
        outputFile = self.temp.joinpath(uuid4().hex)
        completedProcess = run(
            [self.mash, 'sketch', '-p', str(self.p), '-k', str(self.k), '-s', str(self.s), '-l', str(inputFile), '-o', str(outputFile)],
            stdout = DEVNULL, stderr = None
        )
        if completedProcess.returncode:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to run Mash.', file = stderr, flush = True)
            exit(1)
        inputFile.unlink()
        return outputFile.with_suffix('.msh')

    def dist(self, inputFile, outputFile):
        openFile = open(outputFile, 'wb')
        completedProcess = run(
            [self.mash, 'dist', '-p', str(self.p), '-d', str(self.d), str(inputFile), str(inputFile)],
            stdout = openFile, stderr = None
        )
        openFile.close()
        if completedProcess.returncode:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to run Mash.', file = stderr, flush = True)
            exit(1)
        inputFile.unlink()
        return self

    def main(self, inputFiles, outputFile):
        temp = self.parseInputFiles(inputFiles)
        temp = self.sketch(temp)
        self.dist(temp, outputFile)
        return self
