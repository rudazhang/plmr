## Manifold Parzen windows (MParzen) algorithm [@Vincent2002]
## library(data.table)
## library(FNN)

#' Training MParzen model
#'
#' Modified to inlcude the point under consideration and center the local data matrix,
#' instead of center at the point.
#' @param X Matrix of training data, l-by-n.
#' @param d Dimension of manifold.
#' @param k Number of nearest neighbors to use in local covariance analysis.
#' @return A list of a (d, l) matrix `sdev` and an (n, d, l) array `tangent`.
#'         Top d standard deviations and directions of local covariance structure at each point.
#' @references @Vincent2002
MParzenTrain <- function(X, d, k) {
    ## Search nearest neighbors
    neighbor <- FNN::knn.index(X, k)
    LocalModel <- function(i) {
        ## Q-mode PCA: SVD on the data matrix
        pca <- prcomp(X[neighbor[i,],], rank = d, center = TRUE, scale = FALSE, retx = FALSE)
        ## Top singular values (standard deviations) and right singular vectors (tangent directions)
        rbind(sdev = pca$sdev[seq(d)], tangent = pca$rotation)
    }
    l <- nrow(X)
    n <- ncol(X)
    ret <- vapply(seq(l), LocalModel, FUN.VALUE = double(d * (1 + n)))
    dim(ret) <- c(1 + n, d, l)
    list(sdev = matrix(ret[1,,], d),
         tangent = array(ret[-1,,], c(n, d, l)))
}

#' Train MParzen model with fixed-radius nearest neighbors
MParzenRadiusTrain <- function(X, d, r) {
    ## Search fixed-radius nearest neighbors
    neighbor <- dbscan::frNN(X, r, sort = FALSE, approx = 0)$id
    nFewNeighbor <- sum(lengths(neighbor) <= 1L)
    if (nFewNeighbor > 0L) warning(nFewNeighbor, " points with only one or no neighbor.")
    LocalModel <- function(i) {
        ## Q-mode PCA: SVD on the data matrix
        pca <- prcomp(X[c(i, neighbor[[i]]),], rank = d, center = TRUE, scale = FALSE, retx = FALSE)
        ## Top singular values (standard deviations) and right singular vectors (tangent directions)
        rbind(sdev = pca$sdev[seq(d)], tangent = pca$rotation)
    }
    l <- nrow(X)
    n <- ncol(X)
    ret <- vapply(seq(l), LocalModel, FUN.VALUE = double(d * (1 + n)))
    dim(ret) <- c(1 + n, d, l)
    list(sdev = matrix(ret[1,,], d),
         tangent = array(ret[-1,,], c(n, d, l)))
}


#' Testing MParzen model
#'
#' Modified to add noise level $\min(\sigma, s_d)$ to the normal directions,
#' instead of adding $\sigma^2$ to all directions.
#' @param y density evaluation point.
#' @param X matrix of training data.
#' @param model Top d-th order local covariance structure; output of `MParzenTrain()`.
#' @param sigma Maximum noise level added to the normal space of each kernel,
#'        i.e. marginal standard deviation of an isotropic Gaussian. Should be greater than zero.
#' @return Manifold Parzen estimator at the evaluation point.
#' @references @Vincent2002
MParzenTest <- function(y, X, model, sigma) {
    n <- ncol(X)
    d <- nrow(model$sdev)
    ## Pointwise noise level not exceeding the d-th standard deviation.
    noise <- pmin(model$sdev[d,], sigma)
    coef <- n * log(2*pi) + 2 * colSums(log(model$sdev)) + 2 * (n - d) * log(noise)
    yX <- t(X) - y
    dt2 <- rowSums((t(model$sdev)^-2 - noise^-2) *
                   apply(model$tangent, 2, function(.) colSums(. * yX)^2))
    dn2 <- noise^-2 * colSums(yX^2)
    mean(exp(-0.5 * (coef + dt2 + dn2)))
}
## Unit tests
all.equal(MParzenTest(1, matrix(0), list(sdev = matrix(1), tangent = array(1, rep(1,3))), sigma = 1),
          dnorm(1))
all.equal(MParzenTest(1, matrix(c(0, 3), 2), sigma = 1,
                      list(sdev = matrix(c(1, 1), 1), tangent = array(c(1, 1), c(1, 1, 2)))),
          mean(c(dnorm(1), dnorm(2))))
all.equal(MParzenTest(c(1, 0), matrix(c(0, 0), 1), sigma = 1,
                      list(sdev = matrix(1), tangent = array(c(1, 0), c(2, 1, 1)))),
          dnorm(1) * dnorm(0))
all.equal(MParzenTest(c(1, 1), matrix(c(0, 0), 1), sigma = 0.5,
                      list(sdev = matrix(1), tangent = array(c(1, 0), c(2, 1, 1)))),
          dnorm(1) * dnorm(1, sd = 0.5))
all.equal(MParzenTest(c(1, 1), matrix(c(0, 3, 0, 0), 2), sigma = 0.5,
                      list(sdev = matrix(c(1, 1), 1), tangent = array(c(1, 0, 1, 0), c(2, 1, 2)))),
          mean(c(dnorm(1), dnorm(2))) * dnorm(1, sd = 0.5))
all.equal(MParzenTest(c(1, 1), matrix(c(0, 3, 0, 0), 2), sigma = 0.5,
                      list(sdev = matrix(c(1, 1), 1), tangent = array(c(1, 0, 0, 1), c(2, 1, 2)))),
          mean(c(dnorm(1) * dnorm(1, sd = 0.5), dnorm(1) * dnorm(2, sd = 0.5))))
all.equal(MParzenTest(c(1, 1), matrix(c(0, 3, 0, 0), 2), sigma = 0.1,
                      list(sdev = matrix(rep(c(1, 0.5), 2), 2),
                           tangent = array(c(1, 0, 0, 1, 0, 1, -1, 0), c(2, 2, 2)))),
          mean(c(dnorm(1) * dnorm(1, sd = 0.5), dnorm(1) * dnorm(2, sd = 0.5))))

#' Average `MParzenTest()` over all test points
#'
#' Useful for tuning hyperparameters: d, k, sigma.
#' @param Y Matrix of test data, n columns.
#' @param X Matrix of training data, n columns.
#' @param d Dimension of manifold.
#' @param k Number of nearest neighbors to use in local covariance analysis.
#' @param sigma Maximum noise level added to the normal space of each kernel,
#'        i.e. marginal standard deviation of an isotropic Gaussian. Should be greater than zero.
#' @return Average negative log likelihood of the test data given the MParzen model
#'         on the training data and the hyperparameters.
MParzenANLL <- function(Y, X, d, k, sigma) {
    model <- MParzenTrain(X, d, k)
    pMP <- apply(Y, 1, function(.) MParzenTest(., X, model, sigma))
    mean(-log(pMP))
}

MParzenRadiusANLL <- function(Y, X, d, r, sigma) {
    model <- MParzenRadiusTrain(X, d, r)
    pMP <- apply(Y, 1, function(.) MParzenTest(., X, model, sigma))
    mean(-log(pMP))
}

#' Parzen model
#'
#' @param y density evaluation point.
#' @param X matrix of training data.
#' @param sigma Marginal standard deviation of an isotropic Gaussian added to each point.
#'        Should be greater than zero.
#' @return Parzen estimator at the evaluation point.
ParzenTest <- function(y, X, sigma) {
    n <- ncol(X)
    coef <- n * log(2*pi) + 2 * n * log(sigma)
    d2 <- sigma^-2 * colSums((t(X) - y)^2)
    mean(exp(-0.5 * (coef + d2)))
}
## Unit tests
all.equal(ParzenTest(1, matrix(0), 1), dnorm(1))
all.equal(ParzenTest(c(1, 2), matrix(c(0, 0), 1), 1), dnorm(1) * dnorm(2))
all.equal(ParzenTest(c(1, 2), matrix(c(0, 3, 0, 0), 2), 1),
          mean(c(dnorm(1) * dnorm(2), dnorm(2) * dnorm(2))))

#' Average `ParzenTest()` over all test points
#'
#' Useful for tuning hyperparameters: d, k, sigma.
#' @param Y Matrix of test data, n columns.
#' @param X Matrix of training data, n columns.
#' @param sigma Marginal standard deviation of an isotropic Gaussian added to each point.
#'        Should be greater than zero.
#' @return Average negative log likelihood of the test data given the Parzen model
#'         on the training data and the bandwidth.
ParzenANLL <- function(Y, X, sigma) {
    pMP <- apply(Y, 1, function(.) ParzenTest(., X, sigma))
    mean(-log(pMP))
}

#' Leave-one-out cross validation.
#' @param X Data matrix, N-by-n.
#' @param sigma Marginal standard deviation of an isotropic Gaussian added to each point.
#'        Should be greater than zero.
#' @return Average negative log likelihood by Leave-one-out cross validation,
#'         given the Parzen model bandwidth.
LOOCVParzenANLL <- function(X, sigma) {
    pMP <- vapply(seq_len(nrow(X)), function(i) ParzenTest(X[i,], X[-i,], sigma), double(1))
    mean(-log(pMP))
}

#' Optimal Parzen bandwidth that maximizes log-likelihood on test set.
#' @param XTrain Training data.table.
#' @param XTest  Test data.table.
#' @param hlim   Range of bandwidth searching.
#' @return Data.table: minimum, optimal bandwidth; objective, average negative log-likelihood (ANLL).
MLBandwidth <- function(XTrain, XTest, hlim) {
    X <- as.matrix(XTrain)
    Y <- as.matrix(XTest)
    as.data.table(optimize(ParzenANLL, hlim, Y = Y, X = X))
}

#' Optimal Parzen bandwidth that maximizes log-likelihood on test set.
#' @param X    Data, data.table.
#' @param hlim Range of bandwidth searching.
#' @return Data.table: minimum, optimal bandwidth; objective, average negative log-likelihood (ANLL).
LOOCVMLBandwidth <- function(X, hlim) {
    X <- as.matrix(X)
    as.data.table(optimize(LOOCVParzenANLL, hlim, X = X))
}

#' Sample from MParzen model
#'
#' @param N Number of samples.
#' @param X matrix of training data.
#' @param model Top d-th order local covariance structure; output of `MParzenTrain()`.
#' @param sigma Maximum noise level added to the normal space of each kernel,
#'        i.e. marginal standard deviation of an isotropic Gaussian. Can be zero.
#' @return A sample of size N from the Manifold Parzen model.
MParzenSample <- function(N, X, model, sigma) {
    l <- nrow(X)
    n <- ncol(X)
    d <- nrow(model$sdev)
    noise <- pmin(model$sdev[d,], sigma)
    rsdev <- sqrt(model$sdev^2 - rep(noise^2, each = d))
    DT <- data.table(i = sample(seq(l), N, replace = TRUE))
    DT <- DT[, .(Ni = .N), keyby = .(i)]
    DT[, cumN0 := c(0, cumsum(Ni[-.N]))]
    newX <- matrix(double(n * N), n)
    if (d == 1L) {
        sampleIndex <- function(i, Ni, cumN0) {
            newX[,cumN0 + seq(Ni)] <<- X[i,] +
                matrix(model$tangent[,,i], n) %*% (rsdev[,i] * rnorm(Ni)) +
                noise[i] * matrix(rnorm(n * Ni), n)
        }
    } else {
        sampleIndex <- function(i, Ni, cumN0) {
            newX[,cumN0 + seq(Ni)] <<- X[i,] +
                model$tangent[,,i] %*% (rsdev[,i] * matrix(rnorm(d * Ni), d)) +
                noise[i] * matrix(rnorm(n * Ni), n)
        }
    }
    purrr::pwalk(DT, sampleIndex)
    return(newX)
}

#' Sample from MParzen model on tangent spaces
#'
#' 65% time than `MParzenSample(sigma = 0)`.
#' @param N Number of samples.
#' @param X matrix of training data.
#' @param model Top d-th order local covariance structure; output of `MParzenTrain()`.
#' @return A sample of size N from the Manifold Parzen model with no noise.
MParzenSampleNoNoise <- function(N, X, model) {
    l <- nrow(X)
    n <- ncol(X)
    d <- nrow(model$sdev)
    DT <- data.table(i = sample(seq(l), N, replace = TRUE))
    DT <- DT[, .(Ni = .N), keyby = .(i)]
    DT[, cumN0 := c(0, cumsum(Ni[-.N]))]
    newX <- matrix(double(n * N), n)
    if (d == 1L) {
        sampleIndex <- function(i, Ni, cumN0) {
            newX[,cumN0 + seq(Ni)] <<- X[i,] +
                matrix(model$tangent[,,i], n) %*% (model$sdev[,i] * rnorm(Ni))
        }
    } else {
        sampleIndex <- function(i, Ni, cumN0) {
            newX[,cumN0 + seq(Ni)] <<- X[i,] +
                model$tangent[,,i] %*% (model$sdev[,i] * matrix(rnorm(d * Ni), d))
        }
    }
    purrr::pwalk(DT, sampleIndex)
    return(newX)
}

## Performance benchmark
## system.time(mbm <- microbenchmark::microbenchmark(
##                 noise = MParzenSample(1e6, X, model, sigma = 0),
##                 NoNoise = MParzenSampleNoNoise(1e6, X, model),
##                 times = 30)) #18s
## ggplot2::autoplot(mbm)
