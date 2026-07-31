from datetime import datetime
from pathlib import Path
from subprocess import DEVNULL, run
from sys import stderr
from uuid import uuid4


class Skani:
    def __init__(self, skani, t, c, m, temp):
        self.skani = skani
        self.t = t
        self.c = c
        self.m = m
        self.temp = Path(temp)
        return None

    def parseInputFiles(self, files):
        temp = self.temp.joinpath(uuid4().hex)
        openFile = temp.open('w')
        for i in files:
            openFile.write(i + '\n')
        openFile.close()
        return temp

    def dist(self, inputFiles, outputFile):
        temp = self.parseInputFiles(inputFiles)
        completedProcess = run(
            [self.skani, 'dist', '--short-header', '-t', str(self.t), '--ql', str(temp), '--rl', str(temp), '-c', str(self.c), '-m', str(self.m), '-o', outputFile],
            stdout = DEVNULL, stderr = None
        )
        if completedProcess.returncode:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to run skani.', file = stderr, flush = True)
            exit(1)
        temp.unlink()
        return self
