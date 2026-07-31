from ctypes import c_float
from multiprocessing.sharedctypes import RawArray

import numpy

minSequenceLength = 200
maxSequenceLength = 102000
bias = numpy.array([[-1771.7327717505348], [452.3597982867218]], dtype = numpy.float32)
coef = numpy.array(
    [
        [[  5402.674858770911,  -4540.268014185855,     -1250.654519663401], [   -552.5237058125502,   91.72590989069444,       154.4644300086287]],
        [[-2392.3468804633594,   4607.097605942169,     103.30542639078023], [     160.781282896053, -250.97635101186478,     -13.035808035308687]],
        [[  883.2385844456049, -1768.1664944728725,      -5.76253805710207], [  -30.364569804827397,  108.59968171538006,       0.760294626101987]],
        [[-218.31321419172517,   351.8523009275184,     0.1979632815185405], [   3.3568860112960226,  -22.958200140330398,   -0.02732327534502954]],
        [[ 30.001572212577965,  -35.76672354478081, -0.0037983038476760594], [ -0.21517344714850234,   2.4492120522862595,  0.0005457023436405651]],
        [[-1.6970665454593554,   1.466417512938459,  3.152677349252168e-05], [0.0076066653553764395, -0.10541014090773845, -4.637715766033184e-06]]
    ],
    dtype = numpy.float32
) # (6, 2, 3) #


def parseW(w):
    sharedW = RawArray(c_float, w.size * 6)
    w_ = numpy.ndarray(shape = (6, w.size), dtype = numpy.float32, buffer = sharedW)
    numpy.clip(w, minSequenceLength, maxSequenceLength, out = w_[0])
    numpy.log10(w_[0], out = w_[0])
    numpy.power(w_[0], numpy.arange(2, 7)[ : , None], out = w_[1 : ])
    return sharedW


def model(z, i, x, y, w):
    '''
    z: temp array, (5, n)
    i: index, int
    x: kmer frequency, (n, d)
    y: sum of squared row of x, (n, )
    '''

    # calculate w #
    z[3 : 5] = bias
    for j in range(6):
        z[0, : i] = w[j, i]
        z[0, i : ] = w[j, i : ]
        z[1, : i] = w[j, : i]
        z[1, i : ] = w[j, i]
        numpy.multiply(w[j], w[j, i], out = z[2])
        numpy.dot(coef[j], z[ : 3], out = z[ : 2])
        z[3 : 5] += z[ : 2]

    # calculate distance #
    numpy.dot(x[i], x.T, out = z[0])
    z[0] *= -2
    z[0] += y[i]
    z[0] += y
    numpy.maximum(z[0], 0, out = z[0])
    z[0] **= 0.5

    # calculate probability #
    z[0] *= z[3]
    z[0] += z[4]
    with numpy.errstate(over = 'ignore', under = 'ignore'):
        numpy.exp(z[0], out = z[0])
    z[0] += 1

    # log probability #
    numpy.log(z[0], out = z[0])
    numpy.negative(z[0], out = z[0])

    return None