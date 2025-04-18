import numpy


def model(z, i, x, y):
    '''
    z: temp array, (4, n)
    i: index, int
    x: coverage mean, (n, d)
    y: coverage variance, (n, d)
    '''
    z[0] = 0
    for j in range(x.shape[1]):
        numpy.subtract(x[ : , j], x[i, j], out = z[1]) # μ1 - μ2 #
        z[1] **= 2 # (μ1 - μ2)^2 #
        z[1] /= -4 # -(μ1 - μ2)^2 / 4 #
        numpy.add(y[ : , j], y[i, j], out = z[2]) # σ1^2 + σ2^2 #
        z[1] /= z[2] # (μ1 - μ2)^2 / (-4(σ1^2 + σ2^2)) #
        z[0] += z[1]
        numpy.multiply(y[ : , j], y[i, j], out = z[1]) # σ1^2 * σ2^2 #
        with numpy.errstate(divide = 'ignore'):
            numpy.log(z[1 : 3], out = z[1 : 3]) # log(σ1^2 + σ2^2), 2log(σ1σ2) #
        z[1] *= 0.5 # log(σ1σ2) #
        z[1] += numpy.log(2) # log(2σ1σ2) #
        z[1] -= z[2] # log(2σ1σ2) - log(σ1^2 + σ2^2) #
        z[1] *= 0.5 # 0.5 * (log(2σ1σ2) - log(σ1^2 + σ2^2)) #
        z[0] += z[1]
    z[0] /= x.shape[1] # log probability #
    return None


def model_kl(z, i, x, y): # -kl #
    '''
    z: temp array, (4, n)
    i: index, int
    x: coverage mean, (n, d)
    y: coverage variance, (n, d)
    '''

    # kl(n1||n2) = log(σ2) - log(σ1) + (σ1 ** 2 + (μ1 - μ2) ** 2) / (2 * σ2 ** 2) - 0.5 #
    z[0] = 0
    for j in range(x.shape[1]):
        numpy.log(y[ : , j], out = z[1]) # 2 * log(σ2) #
        z[1] -= numpy.log(y[i, j]) # 2 * (log(σ2) - log(σ1)) #

        numpy.subtract(x[i, j], x[ : , j], out = z[2]) # μ1 - μ2 #
        z[2] **= 2 # (μ1 - μ2) ** 2 #
        z[2] += y[i, j] # (μ1 - μ2) ** 2 + σ1 ** 2 #
        z[2] /= y[ : , j] # ((μ1 - μ2) ** 2 + σ1 ** 2) / (σ2 ** 2) #

        z[0] -= z[1]
        z[0] -= z[2]
    z[0] /= 2 * x.shape[1]
    return None


def model_dist(z, i, x, y): # -dist #
    '''
    z: temp array, (4, n)
    i: index, int
    x: coverage mean, (n, d)
    y: coverage variance, (n, d)
    '''
    x = numpy.log(x + 1e-3)
    numpy.einsum('ij,ij->i', x, x, out = z[1])
    numpy.dot(x[i], x.T, out = z[0])
    z[0] *= -2
    z[0] += z[1, i]
    z[0] += z[1]
    numpy.maximum(z[0], 0, out = z[0])
    z[0] **= 0.5
    z[0] /= -(numpy.std(z[0]) + 10 * numpy.finfo(numpy.float32).eps)
    #print('0.90 0.95 0.98 0.99:', numpy.percentile(z[0], [90, 95, 98, 99]), 'max:', numpy.max(z[0]), flush = True)
    return None


def model_gaussian(z, i, x, y):
    '''
    z: temp array, (4, n)
    i: index, int
    x: coverage mean, (n, d)
    y: coverage variance, (n, d)
    '''
    z[0] = 0
    for j in range(x.shape[1]):
        numpy.subtract(x[ : , j], x[i, j], out = z[1])
        z[1] **= 2
        z[1] /= -0.5 * y[i, j]
        z[0] += z[1]
    z[0] /= x.shape[1]
    return None


def model_wasserstein(z, i, x, y):
    '''
    w(n1, n2) = |μ1 - μ2| + (σ1 ** 2 + σ2 ** 2 - 2 * σ1σ2) ** 0.5
              = |μ1 - μ2| + |σ1 - σ2|
    log(p(n1, n2)) = -w(n1, n2)

    z: temp array, (4, n)
    i: index, int
    x: coverage mean, (n, d)
    y: coverage variance, (n, d)
    '''

    z[0] = 0
    for j in range(x.shape[1]):
        numpy.subtract(x[ : , j], x[i, j], out = z[1]) # μ1 - μ2 #
        numpy.abs(z[1], out = z[1]) # |μ1 - μ2| #
        z[0] += z[1]
        numpy.sqrt(y[ : , j], out = z[1]) # σ1 #
        z[1] -= y[i , j] ** 0.5 # σ1 - σ2 #
        numpy.abs(z[1], out = z[1]) # |σ1 - σ2| #
        z[0] += z[1]
    z[0] /= x.shape[1]
    numpy.negative(z[0], z[0])
    return None
