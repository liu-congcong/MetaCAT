import cupy
import numpy
from cupy.linalg import cholesky
from cupyx.scipy.linalg import solve_triangular
from cupyx.scipy.sparse import csr_matrix
from cupyx.scipy.special import digamma
from numpy.typing import NDArray


def normalizeX(x):
    x -= cupy.max(x, axis = 1, keepdims = True)
    sumExp = cupy.sum(cupy.exp(x), axis = 1, keepdims = True)
    logSumExp = cupy.log(sumExp)
    x -= logSumExp
    cupy.exp(x, out = x)
    # x[x < 1e-4] = 0.0 # for fast convergence #
    return None


def choleskyDecomposition(x):
    xCholesky = cholesky(x)
    diagAdd = -6
    while cupy.any(cupy.isnan(xCholesky)): # rare conditions #
        x.flat[ : : x.shape[0] + 1] += 10 ** diagAdd
        xCholesky = cholesky(x)
        diagAdd += 1
    return xCholesky.T # U #


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
        self.Mu = cupy.average(x, weights = w, axis = 0)
        xMinusMu = x - self.Mu
        self.Sigma = xMinusMu.T * w @ xMinusMu
        self.Sigma /= self.n - 1
        self.M = self.d
        return self


    def initializeR(self, y):
        uniqueY = numpy.unique(y)
        data = cupy.zeros(shape = self.n, dtype = cupy.float32)
        indices = cupy.zeros(shape = self.n, dtype = cupy.int32)
        for i in uniqueY:
            if i:
                mask = y == i
                data[mask] = 1
                indices[mask] = self.k
                self.k += 1
        r = csr_matrix(
            (data, indices, cupy.arange(self.n + 1, dtype = cupy.int32)),
            shape = (self.n, self.k),
            dtype = cupy.float32
        )
        return r


    def initializeVariationalParameters(self):
        self.alphaBeta = cupy.empty(shape = (self.k, 2), dtype = cupy.float32)
        self.digammaAlphaBeta = cupy.empty(shape = (self.k, 2), dtype = cupy.float32)
        self.kappa = cupy.empty(shape = self.k, dtype = cupy.float32)
        self.m = cupy.empty(shape = self.k, dtype = cupy.float32)
        self.mu = cupy.empty(shape = (self.k, self.d), dtype = cupy.float32)
        self.precisionCholesky = cupy.empty(shape = (self.k, self.d, self.d), dtype = cupy.float32)
        self.logDeterminantPrecision = cupy.empty(shape = self.k, dtype = cupy.float32)
        return self


    def updateAlphaBeta(self):
        '''
        alpha_k = 1 + Σ{i:1->n}r_ik
        beta_k = Alpha + Σ{i:1->n}(Σ{k':k+1->K}(r_ik'))
        '''
        cupy.add(1.0, self.kr, out = self.alphaBeta[ : self.k, 0])
        # self.alphaBeta[ : , 1] = cupy.hstack((cupy.cumsum(self.kr[ : 0 : -1])[ : : -1], 0)) #
        if self.k > 1:
            cupy.cumsum(self.kr[-1 : 0 : -1], out = self.alphaBeta[self.k - 2 : : -1, 1])
        self.alphaBeta[self.k - 1, 1] = 0
        self.alphaBeta[ : self.k, 1] += self.Alpha
        digamma(self.alphaBeta[ : self.k], out = self.digammaAlphaBeta[ : self.k]) # digamma(α), digamma(β) #
        return self


    def updateKappa(self):
        '''
        kappa_k = Kappa + Σ{i:1->n}r_ik
        '''
        cupy.add(self.Kappa, self.krw, out = self.kappa[ : self.k])
        return self


    def updateM(self):
        cupy.add(self.M, self.kr, out = self.m[ : self.k])
        return self


    def updateMu(self, x, rw):
        '''
        mu_k = (Kappa * Mu + Σ{i:1->n}(r_ik * x_i)) / kappa_k
        '''
        self.kMu = rw.T.dot(x).astype(cupy.float32) # float32 -> float32 #
        cupy.add(self.kMu[ : self.k], self.Kappa * self.Mu, out = self.mu[ : self.k])
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
            sigma = rw[ : , [k]].multiply(xMinusMu).T.dot(xMinusMu).astype(cupy.float32) # float32 -> float32 #
            sigma += Sigma
            muMinusMu = self.kMu[k] - self.Mu
            sigma += self.Kappa * self.krw[k] / self.kappa[k] * cupy.outer(muMinusMu, muMinusMu)
            sigma /= self.m[k]
            sigma.flat[ : : self.d + 1] += self.sigmaAdd
            sigmaCholesky = choleskyDecomposition(sigma) # U #
            precisionCholesky = solve_triangular(sigmaCholesky, cupy.identity(self.d), lower = False)
            '''
            sigma = LU
            inv(sigma) = inv(U)inv(U.T)
            '''
            self.precisionCholesky[k] = precisionCholesky
            self.logDeterminantPrecision[k] = 2.0 * cupy.sum(cupy.log(precisionCholesky.flat[ : : self.d + 1]))
        return self


    def updateR(self, x, w):
        n = x.shape[0]
        if self.neighbors > self.k:
            self.neighbors = self.k
            self.spIndptr[ : n + 1] = cupy.arange(0, n * self.neighbors + 1, self.neighbors, dtype = cupy.int32)

        # E[ln(v_t)] + E{i:1->t-1}[ln(1-v_i)]
        # = digamma(αt) - digamma(αt + βt) + Σ{i:1->t-1}(digamma(βi) - digamma(αi + βi))
        # = digamma(αt) + Σ{i:1->t-1}(digamma(βi)) - Σ{i:1->t}(digamma(αi + βi))
        if self.k > 1:
            self.digammaAlphaBeta[1 : self.k, 0] += cupy.cumsum(self.digammaAlphaBeta[ : self.k - 1, 1])
        self.digammaAlphaBeta[ : self.k, 0] -= cupy.cumsum(digamma(cupy.sum(self.alphaBeta[ : self.k], axis = 1)))

        for i in range(0, n, self.batchSize):
            j = min(n, i + self.batchSize)
            for k in range(self.k):
                self.batchArray[ : j - i, k] = cupy.sum(
                    cupy.square(x[i : j] @ self.precisionCholesky[k] - self.mu[k] @ self.precisionCholesky[k]),
                    axis = 1
                )
                self.batchArray[ : j - i, k] += self.d / self.kappa[k]
                self.batchArray[ : j - i, k] *= w[i : j]
                '''
                log(|sigma / m|) = log(|sigma|) - d * log(m)
                log(|sigma / m|) = - log("logDeterminantPrecision")
                    log(|sigma|) = - log("logDeterminantPrecision") + d * log(m)
                '''
                self.batchArray[ : j - i, k] += cupy.subtract(
                    self.d * cupy.log(self.m[k]) - self.logDeterminantPrecision[k],
                    cupy.sum(
                        digamma(0.5 * (self.m[k] - cupy.arange(self.d)[ : , None])),
                        axis = 0
                    )
                )
                self.batchArray[ : j - i, k] *= -0.5
            self.batchArray[ : j - i, : self.k] += self.digammaAlphaBeta[ : self.k, 0]
            self.spIndices[i : j, : self.neighbors] = cupy.argpartition(self.batchArray[ : j - i, : self.k], -self.neighbors, axis = 1)[ : , -self.neighbors : ]
            self.spData[i : j, : self.neighbors] = cupy.take_along_axis(self.batchArray[ : j - i, : self.k], self.spIndices[i : j, : self.neighbors], axis = 1)
        uniqueIndices = cupy.unique(self.spIndices[ : n, : self.neighbors])
        self.k = uniqueIndices.size
        for i, j in enumerate(uniqueIndices):
            if i != j:
                cupy.equal(self.spIndices[ : n, : self.neighbors], j, out = self.spMask[ : n, : self.neighbors])
                self.spIndices[ : n, : self.neighbors][self.spMask[ : n, : self.neighbors]] = i
        normalizeX(self.spData[ : n, : self.neighbors])

        r = csr_matrix(
            (
                self.spData[ : n, : self.neighbors].flatten(),
                self.spIndices[ : n, : self.neighbors].flatten(),
                self.spIndptr[ : n + 1]
            ),
            shape = (n, self.k),
            dtype = cupy.float32,
            copy = False
        )
        return r


    def isConverged(self, r0, r1):
        r0 -= r1
        return r0.data @ r0.data < self.tolerance


    def predictY(self, x, w):
        n = x.shape[0]
        y = cupy.empty(shape = n, dtype = numpy.int32)
        # E[ln(v_t)] + E{i:1->t-1}[ln(1-v_i)]
        # = digamma(αt) - digamma(αt + βt) + Σ{i:1->t-1}(digamma(βi) - digamma(αi + βi))
        # = digamma(αt) + Σ{i:1->t-1}(digamma(βi)) - Σ{i:1->t}(digamma(αi + βi))
        if self.k > 1:
            self.digammaAlphaBeta[1 : self.k, 0] += cupy.cumsum(self.digammaAlphaBeta[ : self.k - 1, 1])
        self.digammaAlphaBeta[ : self.k, 0] -= cupy.cumsum(digamma(cupy.sum(self.alphaBeta[ : self.k], axis = 1)))

        for i in range(0, n, self.batchSize):
            j = min(n, i + self.batchSize)
            for k in range(self.k):
                self.batchArray[ : j - i, k] = cupy.sum(
                    cupy.square(x[i : j] @ self.precisionCholesky[k] - self.mu[k] @ self.precisionCholesky[k]),
                    axis = 1
                )
                self.batchArray[ : j - i, k] += self.d / self.kappa[k]
                self.batchArray[ : j - i, k] *= w[i : j]
                '''
                log(|sigma / m|) = log(|sigma|) - d * log(m)
                log(|sigma / m|) = - log("logDeterminantPrecision")
                    log(|sigma|) = - log("logDeterminantPrecision") + d * log(m)
                '''
                self.batchArray[ : j - i, k] += cupy.subtract(
                    self.d * cupy.log(self.m[k]) - self.logDeterminantPrecision[k],
                    cupy.sum(
                        digamma(0.5 * (self.m[k] - cupy.arange(self.d)[ : , None])),
                        axis = 0
                    )
                )
                self.batchArray[ : j - i, k] *= -0.5
            self.batchArray[ : j - i, : self.k] += self.digammaAlphaBeta[ : self.k, 0]
            cupy.argmax(self.batchArray[ : j - i, : self.k], axis = 1, out = y[i : j])
        return y


    def main(self, x: NDArray, w: NDArray, y: NDArray) -> NDArray:
        '''
        x: (n, d)
        w: (n,  ), must be sorted from high to low
        y: (n,  ), 0 for unlabeled data
        '''
        x = cupy.asarray(x, dtype = cupy.float32)
        w = cupy.asarray(w, dtype = cupy.float32)
        y = cupy.asarray(y, dtype = cupy.int32)

        self.n, self.d = x.shape
        n = int(cupy.count_nonzero(w >= self.minW))
        if not n:
            n = self.n
        self.initializeHyperParameters(x, w)

        # init e #
        r = self.initializeR(y) # csr #
        self.kr = r.sum(axis = 0, dtype = cupy.float32)[0]
        self.kr += 10 * cupy.finfo(cupy.float32).eps
        rw = r.multiply(w[ : , None])
        rw = rw.tocsc() # csc #
        self.krw = rw.sum(axis = 0, dtype = cupy.float32)[0]
        self.krw += 10 * cupy.finfo(cupy.float32).eps

        # init m #
        self.initializeVariationalParameters()
        self.updateAlphaBeta()
        self.updateKappa()
        self.updateM()
        self.updateMu(x, rw)
        self.updatePrecision(x, rw, 0)

        self.neighbors = min(self.neighbors, self.k)
        self.spData = cupy.empty(shape = (n, self.neighbors), dtype = cupy.float32)
        self.spIndices = cupy.empty(shape = (n, self.neighbors), dtype = cupy.int32)
        self.spIndptr = cupy.arange(0, n * self.neighbors + 1, self.neighbors, dtype = cupy.int32)
        self.spMask = cupy.empty(shape = (n, self.neighbors), dtype = cupy.bool_)
        self.batchArray = cupy.empty(shape = (self.batchSize, self.k), dtype = cupy.float32)

        for i in range(self.maxIterations):
            # e step #
            r_ = r
            r = self.updateR(x[ : n], w[ : n]) # (n, k) #
            if r_.shape == r.shape and self.isConverged(r_, r):
                break
            self.kr = r.sum(axis = 0, dtype = cupy.float32)[0]
            self.kr += 10 * cupy.finfo(cupy.float32).eps
            rw = r.multiply(w[ : n, None])
            rw = rw.tocsc() # csc #
            self.krw = rw.sum(axis = 0, dtype = cupy.float32)[0]
            self.krw += 10 * cupy.finfo(cupy.float32).eps

            # m step #
            self.updateAlphaBeta()
            self.updateKappa()
            self.updateM()
            self.updateMu(x[ : n], rw[ : n])
            self.updatePrecision(x[ : n], rw[ : n], self.Sigma)
        return cupy.asnumpy(self.predictY(x, w))


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
