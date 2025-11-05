import gzip
import os
from gzip import GzipFile
from math import ceil
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


def splitFastaFile(file, n):
    '''
    Split a fasta into small files.
    Parameters:
        file: the path to the fasta file.
        n: the number of output files.
    Return:
        a generator of path of each output file.
    '''

    open4r = open(file, 'rb')
    magicCode = open4r.read(2)
    open4r.close()
    if magicCode == b'\x1f\x8b':
        decompressedFile = uuid4().hex
        open4r = GzipFile(filename = file, mode = 'rb')
        open4w = open(decompressedFile, 'wb', buffering = 10485760)
        while open4w.write(open4r.read(10485760)):
            pass
        open4w.close()
        open4r.close()
        file = decompressedFile
        flag = 1
    else:
        flag = 0
    totalSize = os.path.getsize(file)
    blockSize = ceil(totalSize / n) # blockSize <= totalSize #
    filePosition = 0
    filePosition_ = 0
    open4r = open(file, mode = 'rb')
    while filePosition < totalSize:
        line = open4r.readline()
        filePosition += len(line)
        if line.startswith(b'>'):
            filePosition -= len(line)
            if filePosition > 0:
                open4r.seek(filePosition_, os.SEEK_SET)
                tempFile = uuid4().hex
                open4w = open(tempFile, 'wb')
                while filePosition_ < filePosition:
                    filePosition_ += open4w.write(open4r.read(min(10485760, filePosition - filePosition_)))
                open4w.close()
                yield tempFile
                # filePosition_ will be equal to filePosition, open4r.tell() will be equal to filePosition #
            filePosition = open4r.seek(min(filePosition + blockSize, totalSize), os.SEEK_SET)
    open4r.seek(filePosition_, os.SEEK_SET)
    tempFile = uuid4().hex
    open4w = open(tempFile, 'wb')
    while filePosition_ < filePosition:
        filePosition_ += open4w.write(open4r.read(min(10485760, filePosition - filePosition_)))
    open4w.close()
    yield tempFile
    open4r.close()
    if flag:
        os.remove(file)
    return None


def getNx(file, x):
    lengths = list()
    for _, sequence in readFastaFile(file):
        lengths.append(len(sequence))
    lengths.sort(reverse = True)
    lengthNx = 0.01 * x * sum(lengths)
    length = 0
    for i in lengths:
        length += i
        if length >= lengthNx:
            break
    return i
