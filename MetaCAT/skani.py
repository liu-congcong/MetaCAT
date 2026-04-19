import os
from subprocess import DEVNULL, run
import sys
from uuid import uuid4


def runSkani(skani, t, c, m, inputFiles):
    inputFile = uuid4().hex
    outputFile = uuid4().hex
    openFile = open(inputFile, 'w')
    for i in inputFiles:
        openFile.write(i + '\n')
    openFile.close()
    completedProcess = run(
        [skani, 'dist', '--short-header', '-t', str(t), '--ql', inputFile, '--rl', inputFile, '-c', str(c), '-m', str(m), '-o', outputFile],
        stdout = DEVNULL, stderr = None
    )
    os.remove(inputFile)
    if not completedProcess.returncode:
        print('An error has occurred while running skani.', flush = True)
        sys.exit()
    return outputFile
