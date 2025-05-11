import os
import platform


libraryHash = {
    ('Darwin', 'arm64'): 'Darwin-arm64.dylib',
    ('Linux', 'x86_64'): 'Linux-x86_64.so',
    ('Windows', 'AMD64'): 'Windows-AMD64.dll',
    ('Windows', 'x86'): 'Windows-x86.dll'
}


binaryHash = {
    ('Darwin', 'arm64'): 'Darwin-arm64',
    ('Linux', 'x86_64'): 'Linux-x86_64',
    ('Windows', 'AMD64'): 'Windows-AMD64',
    ('Windows', 'x86'): 'Windows-x86'
}


def findLibrary(name):
    x = libraryHash.get((platform.system(), platform.machine()), None)
    assert x is not None, f'{platform.system()}-{platform.machine()} is not supported.'
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), f'{name}-{x}')


def findBinary(name):
    x = binaryHash.get((platform.system(), platform.machine()), None)
    assert x is not None, f'{platform.system()}-{platform.machine()} is not supported.'
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), f'{name}-{x}')
