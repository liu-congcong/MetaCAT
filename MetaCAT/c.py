from datetime import datetime
from pathlib import Path
from platform import machine, system
from sys import exit, stderr

libraryHash = {
    ('darwin', 'arm64'): 'darwin-arm64.dylib',
    ('linux', 'x86_64'): 'linux-x86_64.so',
    ('windows', 'amd64'): 'windows-amd64.dll',
    ('windows', 'x86'): 'windows-x86.dll'
}


binaryHash = {
    ('darwin', 'arm64'): 'darwin-arm64',
    ('linux', 'x86_64'): 'linux-x86_64',
    ('windows', 'amd64'): 'windows-amd64',
    ('windows', 'x86'): 'windows-x86'
}


def findLibrary(name):
    x = libraryHash.get((system().lower(), machine().lower()), None)
    if x is None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: {system().lower()}-{machine().lower()} is not supported.', file = stderr, flush = True)
        exit(1)
    return str(Path(__file__).with_name(f'{name}-{x}'))


def findBinary(name):
    x = binaryHash.get((system().lower(), machine().lower()), None)
    if x is None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: {system().lower()}-{machine().lower()} is not supported.', file = stderr, flush = True)
        exit(1)
    return str(Path(__file__).with_name(f'{name}-{x}'))
