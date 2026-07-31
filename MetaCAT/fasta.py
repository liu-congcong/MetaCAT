import gzip
import os
from gzip import GzipFile
from math import ceil
from pathlib import Path
from uuid import uuid4


def readFastaFile(file):
    '''
    Parameters:
        file: the path to the (compressed) fasta file.
    Return:
        a generator (sequenceID, sequence)
    '''
    x = list()
    openFile = open(file, 'rb')
    magicCode = openFile.read(2)
    openFile.close()

    if magicCode == b'\x1f\x8b':
        openFile = gzip.open(file, mode = 'rt')
    else:
        openFile = open(file, mode = 'rt')
    for line in openFile:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if x:
                yield (i, ''.join(x))
            i = line.split(' ', maxsplit = 1)[0][1 : ]
            x.clear()
        else:
            x.append(line)
    if x:
        yield (i, ''.join(x))
    openFile.close()
    return None


def splitFastaFile(file, n, temp):
    '''
    Split a fasta into small files.
    Parameters:
        file: the path to the fasta file.
        n: the number of output files.
        temp: path to the temporary directory.
    Return:
        a generator of path of each output file.
    '''
    path = Path(file)
    temp = Path(temp)
    open4r = path.open('rb')
    magicCode = open4r.read(2)
    open4r.close()
    if magicCode == b'\x1f\x8b':
        _path = temp.joinpath(uuid4().hex)
        open4r = GzipFile(filename = path, mode = 'rb')
        open4w = _path.open('wb', buffering = 10485760)
        while open4w.write(open4r.read(10485760)):
            pass
        open4w.close()
        open4r.close()
        path = _path
        flag = 1
    else:
        flag = 0
    totalSize = path.stat().st_size
    blockSize = ceil(totalSize / n) # blockSize <= totalSize #
    filePosition = 0
    _filePosition = 0
    open4r = path.open('rb')
    while filePosition < totalSize:
        line = open4r.readline()
        filePosition += len(line)
        if line.startswith(b'>'):
            filePosition -= len(line)
            if filePosition > 0:
                open4r.seek(_filePosition, os.SEEK_SET)
                _path = temp.joinpath(uuid4().hex)
                open4w = _path.open('wb')
                while _filePosition < filePosition:
                    _filePosition += open4w.write(open4r.read(min(10485760, filePosition - _filePosition)))
                open4w.close()
                yield str(_path)
                # _filePosition will be equal to filePosition, open4r.tell() will be equal to filePosition #
            filePosition = open4r.seek(min(filePosition + blockSize, totalSize), os.SEEK_SET)
    open4r.seek(_filePosition, os.SEEK_SET)
    _path = temp.joinpath(uuid4().hex)
    open4w = _path.open('wb')
    while _filePosition < filePosition:
        _filePosition += open4w.write(open4r.read(min(10485760, filePosition - _filePosition)))
    open4w.close()
    yield str(_path)
    open4r.close()
    if flag:
        path.unlink()
    return None


def getNx(x, nx):
    x.sort(reverse = True)
    y = 0.01 * nx * sum(x)
    z = 0
    for i in x:
        z += i
        if z >= y:
            break
    return i
