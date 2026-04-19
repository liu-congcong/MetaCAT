import numpy
from numpy.typing import NDArray
from scipy.linalg import cholesky, solve_triangular
from scipy.sparse import csr_array
from scipy.special import digamma


def normalizeX(x):
    x -= numpy.max(x, axis = 1, keepdims = True)
    with numpy.errstate(divide = 'ignore', over = 'ignore', under = 'ignore'):
        sumExp = numpy.sum(numpy.exp(x), axis = 1, keepdims = True)
        logSumExp = numpy.log(sumExp)
        x -= logSumExp
        numpy.exp(x, out = x)
    # x[x < 1e-4] = 0.0 #
    return None


def choleskyDecomposition(x):
    diagAdd = -6
    while True:
        try:
            xCholesky = cholesky(x, lower = False)
            break
        except Exception: # rare conditions #
            x.flat[ : : x.shape[0] + 1] += 10 ** diagAdd
            diagAdd += 1
    return xCholesky # U #


class SWDPGMM:

    def __init__(
        self,
        minW: float = 0, neighbors: int = 5, batchSize: int = 100000,
        sigmaAdd: float = 1e-4, maxIterations = 50, tolerance: float = 1e-3
    ) -> None:
        self.minW = minW
        self.neighbors = neighbors
        self.batchSize = batchSize
        self.sigmaAdd = sigmaAdd
        self.maxIterations = maxIterations
        self.tolerance = tolerance
        self.k = 0
        return None


    def initializeHyperParameters(self, x, w):
        self.Alpha = 1.0
        self.Kappa = 1e-4
        self.Mu = numpy.average(x, weights = w, axis = 0)
        xMinusMu = x - self.Mu
        self.Sigma = xMinusMu.T * w @ xMinusMu
        self.Sigma /= self.n - 1
        self.M = self.d
        return self


    def initializeR(self, y):
        uniqueY = numpy.unique(y)
        data = numpy.zeros(shape = self.n, dtype = numpy.float32)
        indices = numpy.zeros(shape = self.n, dtype = numpy.int32)
        for i in uniqueY:
            if i:
                mask = y == i
                data[mask] = 1
                indices[mask] = self.k
                self.k += 1
        r = csr_array(
            (data, indices, numpy.arange(self.n + 1, dtype = numpy.int32)),
            shape = (self.n, self.k),
            dtype = numpy.float32
        )
        return r


    def initializeVariationalParameters(self):
        self.alphaBeta = numpy.empty(shape = (self.k, 2), dtype = numpy.float32)
        self.digammaAlphaBeta = numpy.empty(shape = (self.k, 2), dtype = numpy.float32)
        self.kappa = numpy.empty(shape = self.k, dtype = numpy.float32)
        self.m = numpy.empty(shape = self.k, dtype = numpy.float32)
        self.mu = numpy.empty(shape = (self.k, self.d), dtype = numpy.float32)
        self.precisionCholesky = numpy.empty(shape = (self.k, self.d, self.d), dtype = numpy.float32)
        self.logDeterminantPrecision = numpy.empty(shape = self.k, dtype = numpy.float32)
        return self


    def updateAlphaBeta(self):
        '''
        alpha_k = 1 + Σ{i:1->n}r_ik
        beta_k = Alpha + Σ{i:1->n}(Σ{k':k+1->K}(r_ik'))
        '''
        numpy.add(1.0, self.kr, out = self.alphaBeta[ : self.k, 0])
        # self.alphaBeta[ : , 1] = numpy.hstack((numpy.cumsum(self.kr[ : 0 : -1])[ : : -1], 0)) #
        if self.k > 1:
            numpy.cumsum(self.kr[-1 : 0 : -1], out = self.alphaBeta[self.k - 2 : : -1, 1])
        self.alphaBeta[self.k - 1, 1] = 0
        self.alphaBeta[ : self.k, 1] += self.Alpha
        digamma(self.alphaBeta[ : self.k], out = self.digammaAlphaBeta[ : self.k]) # digamma(α), digamma(β) #
        return self


    def updateKappa(self):
        '''
        kappa_k = Kappa + Σ{i:1->n}r_ik
        '''
        numpy.add(self.Kappa, self.krw, out = self.kappa[ : self.k])
        return self


    def updateM(self):
        numpy.add(self.M, self.kr, out = self.m[ : self.k])
        return self


    def updateMu(self, x, rw):
        '''
        mu_k = (Kappa * Mu + Σ{i:1->n}(r_ik * x_i)) / kappa_k
        '''
        self.kMu = rw.T.dot(x)
        numpy.add(self.kMu[ : self.k], self.Kappa * self.Mu, out = self.mu[ : self.k])
        self.kMu[ : self.k] /= self.krw[ : , None]
        self.mu[ : self.k] /= self.kappa[ : self.k, None]
        return self


    def updatePrecision(self, x, rw, Sigma):
        '''
        sigma#
        = Sigma + ∑(r_i * w_i0 * w_i# * (x_i# - kMu#)(x_i# - kMu#).T)
        + Kappa * ∑(r_i * w_i0 * w_i#) / (Kappa * ∑(r_i * w_i0 * w_i#)) * (kMu# - Mu#)(kMu# - Mu#).T
        '''
        for k in range(self.k):
            xMinusMu = x - self.kMu[k]
            sigma = rw[ : , [k]].multiply(xMinusMu).T @ xMinusMu
            sigma += Sigma
            muMinusMu = self.kMu[k] - self.Mu
            sigma += self.Kappa * self.krw[k] / self.kappa[k] * numpy.outer(muMinusMu, muMinusMu)
            sigma /= self.m[k]
            sigma.flat[ : : self.d + 1] += self.sigmaAdd
            sigmaCholesky = choleskyDecomposition(sigma) # U #
            precisionCholesky = solve_triangular(sigmaCholesky, numpy.identity(self.d), lower = False)
            '''
            sigma = LU
            inv(sigma) = inv(U)inv(U.T)
            '''
            self.precisionCholesky[k] = precisionCholesky
            self.logDeterminantPrecision[k] = 2.0 * numpy.sum(numpy.log(precisionCholesky.flat[ : : self.d + 1]))
        return self


    def updateR(self, x, w):
        n = x.shape[0]
        if self.neighbors > self.k:
            self.neighbors = self.k
            self.spIndptr[ : n + 1] = numpy.arange(0, n * self.neighbors + 1, self.neighbors, dtype = numpy.int32)

        # E[ln(v_t)] + E{i:1->t-1}[ln(1-v_i)]
        # = digamma(αt) - digamma(αt + βt) + Σ{i:1->t-1}(digamma(βi) - digamma(αi + βi))
        # = digamma(αt) + Σ{i:1->t-1}(digamma(βi)) - Σ{i:1->t}(digamma(αi + βi))
        if self.k > 1:
            self.digammaAlphaBeta[1 : self.k, 0] += numpy.cumsum(self.digammaAlphaBeta[ : self.k - 1, 1])
        self.digammaAlphaBeta[ : self.k, 0] -= numpy.cumsum(digamma(numpy.sum(self.alphaBeta[ : self.k], axis = 1)))

        for i in range(0, n, self.batchSize):
            j = min(n, i + self.batchSize)
            for k in range(self.k):
                self.batchArray[ : j - i, k] = numpy.sum(
                    numpy.square(x[i : j] @ self.precisionCholesky[k] - self.mu[k] @ self.precisionCholesky[k]),
                    axis = 1
                )
                self.batchArray[ : j - i, k] += self.d / self.kappa[k]
                self.batchArray[ : j - i, k] *= w[i : j]
                '''
                log(|sigma / m|) = log(|sigma|) - d * log(m)
                log(|sigma / m|) = - log("logDeterminantPrecision")
                    log(|sigma|) = - log("logDeterminantPrecision") + d * log(m)
                '''
                self.batchArray[ : j - i, k] += numpy.subtract(
                    self.d * numpy.log(self.m[k]) - self.logDeterminantPrecision[k],
                    numpy.sum(
                        digamma(0.5 * (self.m[k] - numpy.arange(self.d)[ : , None])),
                        axis = 0
                    )
                )
                self.batchArray[ : j - i, k] *= -0.5
            self.batchArray[ : j - i, : self.k] += self.digammaAlphaBeta[ : self.k, 0]
            self.spIndices[i : j, : self.neighbors] = numpy.argpartition(self.batchArray[ : j - i, : self.k], -self.neighbors, axis = 1)[ : , -self.neighbors : ]
            self.spData[i : j, : self.neighbors] = numpy.take_along_axis(self.batchArray[ : j - i, : self.k], self.spIndices[i : j, : self.neighbors], axis = 1)

        uniqueIndices = numpy.unique(self.spIndices[ : n, : self.neighbors])
        self.k = uniqueIndices.size
        for i, j in enumerate(uniqueIndices):
            if i != j:
                numpy.equal(self.spIndices[ : n, : self.neighbors], j, out = self.spMask[ : n, : self.neighbors])
                self.spIndices[ : n, : self.neighbors][self.spMask[ : n, : self.neighbors]] = i
        normalizeX(self.spData[ : n, : self.neighbors])

        r = csr_array(
            (
                self.spData[ : n, : self.neighbors].flatten(),
                self.spIndices[ : n, : self.neighbors].flatten(),
                self.spIndptr[ : n + 1]
            ),
            shape = (n, self.k),
            dtype = numpy.float32,
            copy = False
        )
        return r


    def isConverged(self, r0, r1):
        r0 -= r1
        return r0.data @ r0.data < self.tolerance


    def predictY(self, x, w):
        n = x.shape[0]
        y = numpy.empty(shape = n, dtype = numpy.int32)

        # E[ln(v_t)] + E{i:1->t-1}[ln(1-v_i)]
        # = digamma(αt) - digamma(αt + βt) + Σ{i:1->t-1}(digamma(βi) - digamma(αi + βi))
        # = digamma(αt) + Σ{i:1->t-1}(digamma(βi)) - Σ{i:1->t}(digamma(αi + βi))
        if self.k > 1:
            self.digammaAlphaBeta[1 : self.k, 0] += numpy.cumsum(self.digammaAlphaBeta[ : self.k - 1, 1])
        self.digammaAlphaBeta[ : self.k, 0] -= numpy.cumsum(digamma(numpy.sum(self.alphaBeta[ : self.k], axis = 1)))

        for i in range(0, n, self.batchSize):
            j = min(n, i + self.batchSize)
            for k in range(self.k):
                self.batchArray[ : j - i, k] = numpy.sum(
                    numpy.square(x[i : j] @ self.precisionCholesky[k] - self.mu[k] @ self.precisionCholesky[k]),
                    axis = 1
                )
                self.batchArray[ : j - i, k] += self.d / self.kappa[k]
                self.batchArray[ : j - i, k] *= w[i : j]
                '''
                log(|sigma / m|) = log(|sigma|) - d * log(m)
                log(|sigma / m|) = - log("logDeterminantPrecision")
                    log(|sigma|) = - log("logDeterminantPrecision") + d * log(m)
                '''
                self.batchArray[ : j - i, k] += numpy.subtract(
                    self.d * numpy.log(self.m[k]) - self.logDeterminantPrecision[k],
                    numpy.sum(
                        digamma(0.5 * (self.m[k] - numpy.arange(self.d)[ : , None])),
                        axis = 0
                    )
                )
                self.batchArray[ : j - i, k] *= -0.5
            self.batchArray[ : j - i, : self.k] += self.digammaAlphaBeta[ : self.k, 0]
            numpy.argmax(self.batchArray[ : j - i, : self.k], axis = 1, out = y[i : j])
        return y


    def main(self, x: NDArray, w: NDArray, y: NDArray) -> NDArray:
        '''
        x: (n, d)
        w: (n,  ), must be sorted from high to low
        y: (n,  ), 0 for unlabeled data
        '''
        self.n, self.d = x.shape
        n = numpy.count_nonzero(w >= self.minW)
        if not n:
            n = self.n
        self.initializeHyperParameters(x, w)

        # init e #
        r = self.initializeR(y) # csr #
        self.kr = r.sum(axis = 0) # array #
        self.kr += 10 * numpy.finfo(numpy.float32).eps # array #
        rw = r.multiply(w[ : , None]) # coo #
        rw = rw.tocsc() # csc #
        self.krw = rw.sum(axis = 0)
        self.krw += 10 * numpy.finfo(numpy.float32).eps

        # init m #
        self.initializeVariationalParameters()
        self.updateAlphaBeta()
        self.updateKappa()
        self.updateM()
        self.updateMu(x, rw)
        self.updatePrecision(x, rw, 0)

        self.neighbors = min(self.neighbors, self.k)
        self.spData = numpy.empty(shape = (n, self.neighbors), dtype = numpy.float32)
        self.spIndices = numpy.empty(shape = (n, self.neighbors), dtype = numpy.int32)
        self.spIndptr = numpy.arange(0, n * self.neighbors + 1, self.neighbors, dtype = numpy.int32)
        self.spMask = numpy.empty(shape = (n, self.neighbors), dtype = numpy.bool_)
        self.batchArray = numpy.empty(shape = (self.batchSize, self.k), dtype = numpy.float32)

        for i in range(self.maxIterations):
            # e step #
            r_ = r
            r = self.updateR(x[ : n], w[ : n]) # (n, k) #
            if r_.shape == r.shape and self.isConverged(r_, r):
                break
            self.kr = r.sum(axis = 0)
            self.kr += 10 * numpy.finfo(numpy.float32).eps
            rw = r.multiply(w[ : n, None]) # coo #
            rw = rw.tocsc() # csc #
            self.krw = rw.sum(axis = 0)
            self.krw += 10 * numpy.finfo(numpy.float32).eps

            # m step #
            self.updateAlphaBeta()
            self.updateKappa()
            self.updateM()
            self.updateMu(x[ : n], rw[ : n])
            self.updatePrecision(x[ : n], rw[ : n], self.Sigma)
        return self.predictY(x, w)


def generateRandomData(n, d, k, initN, random = 0):
    randomGenerator = numpy.random.default_rng(random)
    x = numpy.empty(shape = (n * k, d), dtype = numpy.float32)
    w = numpy.ones(shape = n * k, dtype = numpy.float32)
    y = numpy.zeros(shape = n * k, dtype = numpy.int32)

    means = randomGenerator.random(size = (k, d), dtype = numpy.float32)
    means += numpy.arange(k)[ : , None] * 10
    covariance = numpy.identity(d, dtype = numpy.float32) * 0.1
    for i in range(k):
        x[i * n : (i + 1) * n] = randomGenerator.multivariate_normal(mean = means[i], cov = covariance, size = n)
        y[i * n : i * n + initN] = i + 1
    return (x, w, y)
