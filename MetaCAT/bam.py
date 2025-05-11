import ctypes
import os
from gzip import GzipFile
from hashlib import md5
from math import ceil
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import RawArray, Value
from struct import Struct

import numpy

from .c import findLibrary
from .processbar import ProcessBar

bam = ctypes.cdll.LoadLibrary(findLibrary('bam'))
bam.indexBam.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
bam.indexBam.restype = ctypes.c_int

invalidFlag = 0x4 + 0x100 + 0x200 + 0x400
mapping = numpy.array([4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4], dtype = numpy.uint8)

struct1B = Struct('<B') # uint8 #
struct1c1i = Struct('<ci')
struct1H = Struct('<H') # uint16 #
struct1i = Struct('<i')
struct1I = Struct('<I') # uint32 #
struct2i = Struct('<2i')
struct1Q = Struct('<Q') # uint64 #
struct3Q = Struct('<QQQ')
struct1i3Q = Struct('<iQQQ')
struct2i2B2x2H1i = Struct('<2i2B2x2Hi') # 20 bytes #
struct1i2B2x2H1i = Struct('<i2B2x2Hi') # 16 bytes #
struct3c = Struct('<3c')
structHash = {
    b'A': (0, Struct('<c'), 1), b'B': (2, None, None)     ,
    b'c': (0, Struct('<b'), 1), b'C': (0, Struct('<B'), 1),
    b'f': (0, Struct('<f'), 4), b'H': (1, None, None)     ,
    b'i': (0, struct1i, 4)    , b'I': (0, struct1I, 4)    ,
    b's': (0, Struct('<h'), 2), b'S': (0, Struct('<H'), 2),
    b'Z': (1, None, None)     ,
}


def decodeValue(valueType, alignment, offset, skip = True):
    value = 0
    valueType, typeStruct, deltaOffset = structHash[valueType]
    if valueType == 0: # AcCfiIsS #
        if not skip:
            value = typeStruct.unpack_from(alignment, offset = offset)[0]
        offset += deltaOffset
    elif valueType == 1: # HZ #
        offset = alignment.index(b'\0', offset) + 1
    else: # B #
        subValueType, count = struct1c1i.unpack_from(alignment, offset = offset)
        offset += 5 + structHash[subValueType][2] * count
    return (offset, value)


def readAlignment(file, fileOffset, dataOffset, dataSize):
    openFile = open(file, 'rb')
    openFile.seek(fileOffset, 0)
    gzipFile = GzipFile(fileobj = openFile, mode = 'rb')
    gzipFile.read(dataOffset)
    while dataSize:
        n = struct1i.unpack(gzipFile.read(4))[0]
        alignment = gzipFile.read(n)
        yield (alignment, n)
        dataSize -= 4 + n
    gzipFile.close()
    openFile.close()
    return None


def indexProcess(queue, n, processBar):
    while True:
        file = queue.get()
        if file is None:
            break
        if not os.access(f'{file}.index', os.R_OK):
            assert bam.indexBam(file.encode(), f'{file}.index'.encode()) == 0, f'An error has occured while indexing "{file}".'
        n.acquire()
        n.value += 1
        processBar.plot(n.value)
        n.release()
    return None


def indexBam(files, threads):
    '''
    files: sorted bam files.
    threads: int
    '''
    n = Value(ctypes.c_int64, 0)
    queue = Queue(threads)
    processes = list()
    processBar = ProcessBar(len(files))
    for i in range(threads):
        processes.append(Process(target = indexProcess, args = (queue, n, processBar)))
        processes[-1].start()
    for i in files:
        queue.put(i)
    for process in processes:
        queue.put(None)
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def readIndex(file):
    sequences = list()
    lengths = list()
    fileOffsets = list()
    dataOffsets = list()
    dataSizes = list()
    marker = md5()
    gzipFile = GzipFile(filename = file, mode = 'rb', compresslevel = 9)
    for i in range(struct1I.unpack(gzipFile.read(4))[0]):
        sequence = gzipFile.read(struct1i.unpack(gzipFile.read(4))[0])
        sequences.append(sequence.decode('ascii'))
        length, fileOffset, dataOffset, dataSize = struct1i3Q.unpack(gzipFile.read(28))
        lengths.append(length)
        fileOffsets.append(fileOffset)
        dataOffsets.append(dataOffset)
        dataSizes.append(dataSize)
        marker.update(sequence)
        marker.update(length.to_bytes(length = 4, byteorder = 'little'))
    gzipFile.close()
    return (marker.digest(), sequences, lengths, fileOffsets, dataOffsets, dataSizes)


def readIndexProcess(queue, x, n, m, marker):
    '''
    n: #sequences
    m: #files
    '''
    x = numpy.ndarray(shape = (3, n, m), dtype = numpy.int64, buffer = x)
    while True:
        i, file = queue.get()
        if i is None:
            break
        markerI, sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndex(file)
        assert markerI == marker, 'All bam files should have the same header.'
        x[0, : , i] = fileOffsets[ : -1]
        x[1, : , i] = dataOffsets[ : -1]
        x[2, : , i] = dataSizes[ : -1]
    return None


def readIndices(files, threads):
    m = len(files)
    marker, sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndex(f'{files[0]}.index')
    n = len(sequences) - 1
    X = RawArray(ctypes.c_int64, 3 * n * m)
    x = numpy.ndarray(shape = (3, n, m), dtype = numpy.int64, buffer = X)
    x[0, : , 0] = fileOffsets[ : -1]
    x[1, : , 0] = dataOffsets[ : -1]
    x[2, : , 0] = dataSizes[ : -1]
    queue = Queue(threads)
    processes = list()
    for i in range(threads):
        processes.append(Process(target = readIndexProcess, args = (queue, X, n, m, marker)))
        processes[-1].start()
    for i, file in enumerate(files[1 : ], start = 1):
        queue.put((i, f'{file}.index'))
    for process in processes:
        queue.put((None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return (sequences[ : -1], lengths[ : -1], x[0], x[1], x[2])


def readNM(alignment, alignmentSize, offset):
    nm = 0
    while offset < alignmentSize: # optional fields #
        tagA, tagB, valueType = struct3c.unpack_from(alignment, offset = offset)
        offset += 3
        if tagA == b'N' and tagB == b'M':
            offset, nm = decodeValue(valueType, alignment, offset, skip = False)
            break
        else:
            offset = decodeValue(valueType, alignment, offset)[0]
    return nm


def decodeRead(file, fileOffset, dataOffset, dataSize, minMAPQ, minIdentity):
    ungappedRegions = list()
    for alignment, alignmentSize in readAlignment(file, fileOffset, dataOffset, dataSize):
        referencePosition, readIDLength, mapq, nCigars, flag, readLength = struct1i2B2x2H1i.unpack_from(alignment, offset = 4)
        if (flag & invalidFlag == 0) and (mapq >= minMAPQ) and nCigars:
            alignedLength = 0
            ungappedLength = 0
            nm = 0
            offset = 32 + readIDLength # 32 + readIDLength #
            cigarOffset = 4 * nCigars
            read = alignment[offset + cigarOffset : offset + cigarOffset + ceil(readLength * 0.5)]
            readPosition = 0
            for cigar in struct1I.iter_unpack(alignment[offset : offset + cigarOffset]):
                operations = cigar[0] >> 4
                operation = cigar[0] & 0x0F # MIDNSHP=X #
                if operation in {0, 7, 8}: # =MX #
                    ungappedRegions.append((referencePosition, referencePosition + operations, readPosition, readPosition + operations))
                    referencePosition += operations
                    readPosition += operations
                    ungappedLength += operations
                    alignedLength += operations
                elif operation == 2: # D #
                    referencePosition += operations
                    alignedLength += operations
                    nm -= operations
                elif operation == 1: # I #
                    readPosition += operations
                    alignedLength += operations
                    nm -= operations
                elif operation == 3: # N #
                    referencePosition += operations
                elif operation == 4: # S #
                    readPosition += operations
            offset += cigarOffset + ceil(readLength * 0.5) + readLength # cigar + seq + qual #
            nm += readNM(alignment, alignmentSize, offset)
            if ungappedLength - max(0, nm) >= alignedLength * minIdentity:
                read = numpy.asarray(list(read), dtype = numpy.uint8)[ : , None]
                yield (mapping[numpy.hstack((read >> 4, read & 0x0F), dtype = numpy.uint8).flatten()], ungappedRegions)
            ungappedRegions.clear()
    return None


def getUngappedRegions(file, fileOffset, dataOffset, dataSize, minMAPQ, minIdentity):
    ungappedRegions = list()
    for alignment, alignmentSize in readAlignment(file, fileOffset, dataOffset, dataSize):
        position, readIDLength, mapq, nCigars, flag, readLength = struct1i2B2x2H1i.unpack_from(alignment, offset = 4)
        if flag & 0x4 == 0:
            if (flag & invalidFlag == 0) and (mapq >= minMAPQ) and nCigars:
                offset = 32 + readIDLength # 32 + readIDLength #
                cigarOffset = 4 * nCigars
                alignedLength = 0
                ungappedLength = 0
                nm = 0
                for cigar in struct1I.iter_unpack(alignment[offset : offset + cigarOffset]):
                    operations = cigar[0] >> 4
                    operation = cigar[0] & 0x0F # MIDNSHP=X #
                    if operation in {0, 7, 8}: # =MX: consume reference #
                        ungappedRegions.append((position, position + operations))
                        ungappedLength += operations
                        position += operations
                        alignedLength += operations
                    elif operation == 2: # D: consume reference #
                        position += operations
                        alignedLength += operations
                        nm -= operations
                    elif operation == 1:
                        alignedLength += operations
                        nm -= operations
                    elif operation == 3: # N: consume reference #
                        position += operations
                offset += cigarOffset + ceil(readLength * 0.5) + readLength # cigar + seq + qual #
                nm += readNM(alignment, alignmentSize, offset)
                if ungappedLength - max(0, nm) < alignedLength * minIdentity:
                    ungappedRegions.clear()
            elif flag & 0x100:
                continue
            alignmentType = 1
        else:
            alignmentType = 0
        yield (alignmentType, ungappedRegions) # alignmentType: 0 - unmapped alignment, 1 - mapped alignment #
        ungappedRegions.clear()
    return None
