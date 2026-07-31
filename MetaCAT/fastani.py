from datetime import datetime
from pathlib import Path
from subprocess import DEVNULL, run
from sys import exit, stderr
from uuid import uuid4


class Fastani:
    def __init__(self, fastani, t, k, fragLen, temp):
        self.fastani = fastani
        self.t = t
        self.k = k
        self.fragLen = fragLen
        self.temp = Path(temp)
        return None

    def parseInputFiles(self, files):
        temp = self.temp.joinpath(uuid4().hex)
        openFile = temp.open('w')
        for i in files:
            openFile.write(i + '\n')
        openFile.close()
        return temp

    def main(self, inputFiles, outputFile):
        temp = self.parseInputFiles(inputFiles)
        completedProcess = run(
            [self.fastani, '--rl', str(temp), '--ql', str(temp), '-k', str(self.k), '-t', str(self.t), '--fragLen', str(self.fragLen), '--minFraction', '0', '--maxRatioDiff', '100', '-o', outputFile],
            stdout = DEVNULL, stderr = None
        )
        if completedProcess.returncode:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: Failed to fastANI.', file = stderr, flush = True)
            exit(1)
        temp.unlink()
        return self
