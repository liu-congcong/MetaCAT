import os
from subprocess import DEVNULL, run
import sys
from uuid import uuid4


def runMash(mash, d, p, k, s, inputFiles):
    inputFile = uuid4().hex
    outputFile = uuid4().hex
    openFile = open(inputFile, 'w')
    for i in inputFiles:
        openFile.write(i + '\n')
    openFile.close()
    completedProcess = run(
        [mash, 'sketch', '-p', str(p), '-k', str(k), '-s', str(s), '-l', inputFile, '-o', outputFile],
        stdout = DEVNULL, stderr = None
    )
    os.remove(inputFile)
    if completedProcess.returncode:
        print('An error has occurred while running Mash.', flush = True)
        sys.exit()
    openFile = open(outputFile, 'wb')
    completedProcess = run(
        [mash, 'dist', '-p', str(p), '-d', str(d), outputFile + '.msh', outputFile + '.msh'],
        stdout = openFile, stderr = None
    )
    openFile.close()
    os.remove(outputFile + '.msh')
    if completedProcess.returncode:
        print('An error has occurred while running Mash.', flush = True)
        sys.exit()
    return outputFile
