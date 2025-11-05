import os
from shutil import which


def isReadable(path):
    '''
    Detect whether the path is readable.
    Parameters:
        path: the abspath to the file.
    Return:
        None
    '''
    assert os.access(path, os.R_OK), f'\"{path}\" is not readable.'
    return None


def isWriteable(path, pathType = 'file'):
    '''
    Detect whether the path is writeable.
    Parameters:
        path: the abspath to the file / dirctory.
    Return:
        None
    '''
    try:
        if pathType == 'file':
            open4w = open(path, 'wb')
            open4w.close()
            os.remove(path)
        else:
            os.makedirs(path, exist_ok = True)
    except Exception:
        assert False, f'\"{path}\" is not writeable.'
    return None


def findExecutablePath(path, parameter):
    path = which(path, mode = os.X_OK)
    assert path is not None, f'"{parameter}" should be specified.'
    return os.path.realpath(path)
