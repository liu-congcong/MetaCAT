import gzip
import os
from ctypes import c_int64
from hashlib import md5
from math import ceil
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value
from subprocess import DEVNULL, run
from uuid import uuid4

from .processbar import ProcessBar


def is_gzipped(input_file):
    open_file = open(input_file, 'rb')
    magic_code = open_file.read(2)
    open_file.close()
    return magic_code == b'\x1f\x8b'


def check_database(hmm: str, n: int = os.cpu_count()) -> tuple[bool, list[str]]:
    flag = True
    files = list()
    hmm = os.path.abspath(hmm)
    x = md5(os.stat(hmm).st_mtime_ns.to_bytes(16, byteorder = 'little'))
    x.update(n.to_bytes(16, byteorder = 'little'))
    database = os.path.join(os.path.dirname(hmm), x.hexdigest())
    if os.access(database, os.R_OK):
        files.append(database)
        for i in range(n):
            if (os.access(f'{database}.{i}.h3f', os.R_OK) and os.access(f'{database}.{i}.h3i', os.R_OK) and os.access(f'{database}.{i}.h3m', os.R_OK) and os.access(f'{database}.{i}.h3p', os.R_OK)):
                files.append(f'{database}.{i}')
            elif not os.access(f'{database}.{i}', os.R_OK):
                flag = False
                break
    else:
        flag = False
    return (flag, files)


def create_database(hmmpress: str, hmm: str, n: int = os.cpu_count()) -> None:
    hmm = os.path.abspath(hmm)
    x = md5(os.stat(hmm).st_mtime_ns.to_bytes(16, byteorder = 'little'))
    x.update(n.to_bytes(16, byteorder = 'little'))
    database = os.path.join(os.path.dirname(hmm), x.hexdigest())

    x = list()
    X = list()
    y = list()
    if is_gzipped(hmm): # gzipped file #
        openFile = gzip.open(hmm, mode = 'rt')
    else:
        openFile = open(hmm, 'r')
    for line in openFile:
        line = line.rstrip('\n')
        x.append(line + '\n')
        if line.startswith('NAME'):
            profile = line.split()[1]
            ga = ['0', '0']
            tc = ['0', '0']
            nc = ['0', '0']
        elif line.startswith('ACC'): # replace NAME with ACC #
            profile = line.split()[1]
        elif line.startswith('GA'):
            ga = line.split()[1 : ]
        elif line.startswith('TC'):
            tc = line.split()[1 : ]
        elif line.startswith('NC'):
            nc = line.split()[1 : ]
        elif line.startswith('//'):
            X.append(''.join(x))
            x.clear()
            y.append(f'{profile}\t{ga[0]}\t{ga[1]}\t{tc[0]}\t{tc[1]}\t{nc[0]}\t{nc[1]}\n')
    openFile.close()

    openFile = gzip.open(database, mode = 'wt', compresslevel = 9)
    openFile.write('profile\tga sequence\tga domain\ttc sequence\ttc domain\tnc_sequence\tnc_domain\n')
    openFile.writelines(y)
    openFile.close()

    step = ceil(len(X) / n)
    for i in range(n):
        openFile = open(f'{database}.{i}', 'w')
        x = X[i * step : (i + 1) * step]
        openFile.writelines(x)
        openFile.close()
        if x:
            run_hmmpress(hmmpress, f'{database}.{i}')
            os.remove(f'{database}.{i}')
    return None


def run_hmmpress(hmmpress, hmm):
    '''
    Parameters:
        hmmpress: hmmpress
        hmm: HMM file
    Return:
        None
    '''
    completedProcess = run([hmmpress, '-f', hmm], stdout = DEVNULL, stderr = None)
    assert not completedProcess.returncode, 'An error has occurred while running hmmpress.'
    return None


def read_score_file(score_file):
    '''
    Parameters:
        score_file: score file
    Return:
        profile2score: dict()
            profile: str
            scores: (sequence_score, domain_score)
    '''
    profile2score = dict()
    open_file = gzip.open(score_file, mode = 'rt')
    open_file.readline() # profile, ga sequence, ga domain, tc sequence, tc domain, nc_sequence, nc_domain #
    for line in open_file:
        lines = line.rstrip('\n').split('\t')
        if lines[0].startswith('TIGR') and lines[3] and lines[4]:
            profile2score[lines[0]] = (float(lines[3]), float(lines[4]))
        elif lines[1] and lines[2]: # for all profiles #
            profile2score[lines[0]] = (float(lines[1]), float(lines[2]))
        elif lines[3] and lines[4]: # for all profiles #
            profile2score[lines[0]] = (float(lines[3]), float(lines[4]))
        elif lines[5] and lines[6]: # for all profiles #
            profile2score[lines[0]] = (float(lines[5]), float(lines[6]))
        else: # for all profiles #
            profile2score[lines[0]] = None
    open_file.close()
    return profile2score


def dump_database(file, n):
    profile2score = dict()
    files = list()

    x = list()
    profile = ''
    ni = 0
    open4r = gzip.open(file, mode = 'rt')
    for line in open4r:
        line = line.rstrip('\n')
        x.append(line + '\n')
        if line.startswith('ACC'):
            profile = line.split()[1]
        elif profile.startswith('PF') and line.startswith('GA'):
            lines = line.split()
            profile2score[profile] = (float(lines[1]), float(lines[2]))
        elif profile.startswith('TIGR') and line.startswith('TC'):
            lines = line.split()
            profile2score[profile] = (float(lines[1]), float(lines[2]))
        elif line.startswith('//'):
            ni += 1
            if ni == n:
                files.append(uuid4().hex)
                open4w = open(files[-1], 'w')
                open4w.writelines(x)
                open4w.close()
                x.clear()
                ni = 0
    open4r.close()
    if ni:
        files.append(uuid4().hex)
        open4w = open(files[-1], 'w')
        open4w.writelines(x)
        open4w.close()
    return (profile2score, files)


def hmmsearchProcess(queue, hmmsearch, fastaFile, n, N):
    processBar = ProcessBar(N)
    while True:
        hmmFile, hitFile = queue.get()
        if hmmFile is None:
            break
        completedProcess = run(
            [hmmsearch, '--noali', '--notextw', '--cpu', '1', '--tblout', hitFile, hmmFile, fastaFile],
            stdout = DEVNULL, stderr = None
        )
        assert not completedProcess.returncode, 'An error has occurred while running hmmsearch.'
        n.acquire()
        n.value += 1
        processBar.plot(n.value)
        n.release()
    return None


def createProcesses(hmmsearch, fastaFile, n, threads):
    queue = Queue(threads)
    processes = list()
    N = Value(c_int64, 0)
    for i in range(threads):
        processes.append(Process(target = hmmsearchProcess, args = (queue, hmmsearch, fastaFile, N, n)))
        processes[-1].start()
    return (queue, processes, N)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def runHmmsearch(hmmsearch, hmmFiles, fastaFile, threads):
    n = len(hmmFiles)
    hitFiles = list()
    queue, processes, N = createProcesses(hmmsearch, fastaFile, n, threads)
    for hmmFile in hmmFiles:
        hitFiles.append(uuid4().hex)
        queue.put((hmmFile, hitFiles[-1]))
    freeProcesses(queue, processes)
    outputFile = uuid4().hex
    open4w = open(outputFile, 'wb', buffering = 10485760)
    for hitFile in hitFiles:
        open4r = open(hitFile, 'rb', buffering = 10485760)
        while open4w.write(open4r.read(10485760)):
            pass
        open4r.close()
        os.remove(hitFile)
    open4w.close()
    return outputFile


def get_valid_hits_(profile2score, profile_coverage, target_coverage, input_file):
    open_file = open(input_file, 'r')
    # target name, accession, tlen, query name, accession, qlen, E-value, score, bias, #, of, c-Evalue, i-Evalue, score, bias, from, to, from, to, from, to, acc, description of target #
    for line in open_file:
        if not line.startswith('#'):
            lines = line.rstrip('\n').split(maxsplit = 19)
            '''
            t name = lines[0]
            t len = int(lines[2])
            q name = lines[3]
            q accession = lines[4]
            q len: int(lines[5])
            sequence score: float(lines[7])
            domain score = float(lines[13])
            hmm from = int(lines[15])
            hmm to = int(lines[16])
            ali from = int(lines[17])
            ali to = int(lines[18])
            '''
            if (int(lines[16]) - int(lines[15]) + 1) / int(lines[5]) >= profile_coverage:
            # if (int(lines[18]) - int(lines[17]) + 1) / int(lines[2]) >= target_coverage:
                profile = lines[4] if lines[4] != '-' else lines[3]
                score_threshold = profile2score[profile]
                score = float(lines[7])
                if (score_threshold is None) or (score >= score_threshold[0]):
                    yield(lines[0], profile, score)
    open_file.close()
    return None


def get_valid_hits(profile2score, inputFile):
    openFile = open(inputFile, 'r')
    # target name, accession, query name, accession, E-value, score, bias, E-value, score, bias, exp, reg, clu, ov, env, dom, rep, inc, description of target #
    for line in openFile:
        if not line.startswith('#'):
            lines = line.rstrip('\n').split(maxsplit = 9)
            profile = lines[3] if lines[3] != '-' else lines[2]
            minScore = profile2score[profile]
            score = float(lines[5])
            if (minScore is None) or (score >= minScore[0]):
                yield(lines[0], profile, score)
    openFile.close()
    return None
