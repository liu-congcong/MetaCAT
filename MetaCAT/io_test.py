import os
from datetime import datetime
from pathlib import Path
from shutil import which
from sys import exit, stderr


def isReadable(x):
    '''
    Detect whether the path is readable.
    Parameters:
        x: the abspath to the file.
    Return:
        None
    '''
    if not Path(x).is_file():
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: \"{x}\" is not readable.', file = stderr, flush = True)
        exit(1)
    return None


def isWriteable(path, pathType = 'file'):
    '''
    Detect whether the path is writeable.
    Parameters:
        path: the abspath to the file / dirctory.
    Return:
        None
    '''
    path = Path(path)
    if pathType == 'file':
        try:
            path.touch()
            path.unlink()
        except Exception:
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: \"{path}\" is not writeable.', file = stderr, flush = True)
            exit(1)
    else:
        if not (path.is_dir() and os.access(path, os.W_OK)):
            print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: \"{path}\" is not writeable.', file = stderr, flush = True)
            exit(1)
    return None


def findExecutablePath(x):
    y = which(x)
    if y is None:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Error: \"{x}\" is not executable.', file = stderr, flush = True)
        exit(1)
    return str(Path(y).resolve())
