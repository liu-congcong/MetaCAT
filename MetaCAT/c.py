import os
import platform


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
    x = libraryHash.get((platform.system().lower(), platform.machine().lower()), None)
    assert x is not None, f'{platform.system().lower()}-{platform.machine().lower()} is not supported.'
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), f'{name}-{x}')


def findBinary(name):
    x = binaryHash.get((platform.system().lower(), platform.machine().lower()), None)
    assert x is not None, f'{platform.system().lower()}-{platform.machine().lower()} is not supported.'
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), f'{name}-{x}')
