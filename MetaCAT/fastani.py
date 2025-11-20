import os
from subprocess import DEVNULL, run
from uuid import uuid4


def runFastani(fastani, t, k, fragLen, inputFiles):
    inputFile = uuid4().hex
    outputFile = uuid4().hex
    openFile = open(inputFile, 'w')
    for i in inputFiles:
        openFile.write(i + '\n')
    openFile.close()
    completedProcess = run(
        [fastani, '--rl', inputFile, '--ql', inputFile, '-k', str(k), '-t', str(t), '--fragLen', str(fragLen), '--minFraction', '0', '--maxRatioDiff', '100', '-o', outputFile],
        stdout = DEVNULL, stderr = None
    )
    os.remove(inputFile)
    assert not completedProcess.returncode, 'An error has occurred while running fastANI.'
    return outputFile
