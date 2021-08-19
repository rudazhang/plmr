## Ridge Estimation
## library(data.table)
## library(purrr)
## library(irlba)
## library(RSpectra)
## library(dbscan)

#' Silverman's rule of thumb bandwidths
#'
#' The optimal bandwidths for Gaussian KDE on a Gaussian random vector
#' in terms of the mean integrated squared error (MISE).
#' @param N Sample size.
#' @param n Dimension of observation.
#' @param sd Data standard deviations.
#' @references @Silverman1986
silverman <- function(N, n, sd) (N * (n+2) / 4)^(-1/(n+4)) * sd

#' Distance to other points in a set
#' @param X Data, an N-by-n matrix.
mindist <- function(X) as.vector(FNN::knn.dist(X, 1))

#' Detect Outers
#'
#' @param X Data, an N-by-n matrix.
#' @param sds Multiple of standard deviations away from the median distance to point cloud
#' for a data point to be outlier.
#' @return Logical vector indicating outliers in data for a smoother KDE.
#' @references @Genovese2014
outlier <- function(X, sds = 6) {
    d <- mindist(X)
    d > median(d) + sds * sd(d)
    ## n <- nrow(X)
    ## N <- ncol(X)
    ## if (is.null(t)) t <- 1.5 / N
    ## if (is.null(h)) h <- silverman(N, n, apply(X, 1, sd))
    ## Kernel density
    ## Dh <- as.matrix(dist(t(X / h)))
    ## K <- exp(- Dh^2 / 2)
    ## p <- rowMeans(K)
    ## which(p < t)
}


#' (KDE-based) Mean shift
#'
#' @param X       Data, an N-by-n matrix.
#' @param h       Bandwidth for isotropic Gaussian kernel density estimation, a scalar.
#' @param Y0      Initial points, an M-by-n matrix, defaults to X.
#' @param minD    Convergence criterion, length of mean shift vector.
#' @param maxIter Maximum number of iterations.
#' @return A list: `Y`, the final points; `updatedPoints`, the number of updated points per update.
#' @references @Comaniciu2002
MeanShift <- function(X, h, Y0 = NULL, minD = 1e-12, maxIter = 100L, silent = FALSE) {
    stopifnot(is.matrix(X))
    n <- ncol(X)
    N <- nrow(X)
    if(is.null(Y0)) Y0 <- X
    else stopifnot(ncol(Y0) == n)
    M <- nrow(Y0)
    ## Turn to data-contiguous storage for ease of computation.
    tX <- t(X)
    Yd <- rbind(t(Y0), delta = rep(NA_real_, M))
    updatePoint <- function(i) {
        y <- Yd[-nrow(Yd), i]
        Z <- (tX - y) / h
        ## Kernel density
        c <- exp(-0.5 * colSums(Z^2))
        avgc <- mean(c)
        avgcZ <- drop(Z %*% c / length(c))
        ## Mean-shift (update) vector, based on log of Gaussian KDE
        m <- h * avgcZ / avgc
        c(y + m, delta = sqrt(sum(m^2)))
    }
    if (!silent) message("Updating points: ", M)
    iterCounter <- 1L
    updatedPoints <- integer(maxIter)
    updatedPoints[[iterCounter]] <- M
    Yd <- vapply(seq_len(M), updatePoint, double(n + 1))
    cols <- which(Yd["delta",] > minD)
    while(length(cols) > 0L & iterCounter < maxIter) {
        if (!silent) message("...", length(cols))
        iterCounter <- iterCounter + 1L
        updatedPoints[[iterCounter]] <- length(cols)
        Yd[,cols] <- vapply(cols, updatePoint, double(n + 1))
        cols <- which(Yd["delta",] > minD)
    }
    list(Y = Yd[-nrow(Yd),], updatedPoints = updatedPoints)
}

#' KDE-based Subspace Constrained Mean Shift
#'
#' Compuational complexity per step: O(M * N * n^3), where M is for every query point,
#' N for every data point, n^3 for eigen-decomposition.
#' If partial/truncated eigen-decomposition is used, the n^3 term would become d*n^2
#' (e.g. the Lanczos algorithm, which further reduces by the sparsity of the matrix).
#'
#' In the return, `lambdaNormal` and `eigengap` are both useful for checking
#' if the final point is a ridge point.
#' @param Y0      Initial points, an n-by-M matrix.
#' @param X       Data, an n-by-N matrix.
#' @param d       Dimension of the density ridge, an integer.
#' @param h       Bandwidth for isotropic Gaussian kernel density estimation, a scalar.
#' @param log     Whether to use the log Gaussian KDE. Default to TRUE.
#' @param r       Radius of nearest neighbors to use in evaluating density and derivatives.
#'                Default to use all data.
#' @param omega   Stepsize relaxation factor.
#' @param maxStep Maximum stepsize. You probably don't need to use both `omega` and `maxStep`.
#' @param minCos  Convergence criterion, cosine of the angle between gradient and
#'                its component in the "local normal space".
#' @param maxIter Maximum number of iterations.
#' @param returnV Whether to return the eigenvector at Y; an (n,n,M) array.
#' @return A list: `Y`, the final points; `cosNormal`, convergence score;
#'                 `lambdaNormal`, largest Hessian eigenvalue in the normal space;
#'                 `eigengap`, d-th eigengap of the Hessian;
#'                 `updatedPoints`, the number of updated points per update;
#' @references [@Ozertem2011]
#' @export
scms <- function(Y0, X, d, h, log = TRUE, r = NULL, omega = 1, maxStep = NULL,
                 minCos = 0.01, maxIter = 100L, returnV = FALSE, silent = FALSE) {
    ## N <- ncol(X)
    n <- nrow(X)
    stopifnot(nrow(Y0) == n)
    M <- ncol(Y0)
    if (is.null(r)) neighbor <- NULL
    else {
        ## Fixed-radius nearest neighbors to initial points: transpose to "data matrices".
        neighbor <- dbscan::frNN(t(X), r, query = t(Y0), sort = FALSE, approx = 0)$id
        nFewNeighbor <- sum(lengths(neighbor) <= 1L)
        if (nFewNeighbor > 0L) warning(nFewNeighbor, " points with only one or no neighbor.")
    }
    extraVars <- 3L
    ## Prepare coordinate-cosine matrix for updates later.
    Yc <- rbind(Y0, matrix(rep(NA_real_, M * extraVars), extraVars))
    if (returnV) V <- array(rep(NA_real_, n * n * M), c(n, n, M))
    ## Update a point by index.
    ## Return new coordinates and cosine between gradient and local normal space.
    updatePoint <- function(i) {
        y <- Yc[seq_len(n), i]
        idNN <- neighbor[[i]]
        if (is.null(idNN)) Z <- (X - y) / h
        else if (length(idNN) > 1L) Z <- (X[, idNN] - y) / h
        else if (length(idNN) == 1L) Z <- (matrix(X[, idNN], n) - y) / h
        else return(c(y, cosNormal = 0))
        ## Kernel density.
        c <- exp(-0.5 * colSums(Z^2))
        avgc <- mean(c)
        avgcZ <- drop(Z %*% c / length(c))
        if (log) {
            ## Log of Gaussian KDE, top d eigenvectors of the Hessian.
            rootpcZ <- rep(sqrt(c / sum(c)), each = n) * Z
            sumpcZ <- avgcZ / avgc
            spectr <- eigen(tcrossprod(rootpcZ) - tcrossprod(sumpcZ), symmetric = TRUE)
            Vd <- spectr$vectors[, seq(d)]
            ## Hessian eigenvalues should minus 1 to account for the omitted identity matrix.
            lambdaNormal <- spectr$values[[d + 1L]] - 1
            eigengap <- spectr$values[[d]] - spectr$values[[d + 1L]]
            ## Partial eigen-decomposition
            ## Input to `RSpectra::eigs_sym()` must be a square matrix of size at least 3.
            ## spectr <- RSpectra::eigs_sym(tcrossprod(rootpcZ) - tcrossprod(sumpcZ), d)
            ## Vd <- spectr$vectors
        } else {
            ## Gaussian KDE, top d eigenvectors of the Hessian.
            rootcZ <- rep(sqrt(c), each = n) * Z
            sVd <- svd(rootcZ, nu = d, nv = 0)
            ## Partial SVD
            ## If d >= n/2: "You're computing too large a percentage of total singular values..."
            ## sVd <- irlba::irlba(rootcZ, nu = d, nv = 0)
            Vd <- sVd$u
        }
        ## NOTE: special case for d == 1L
        dim(Vd) <- c(n, d)
        ## Mean-shift (update) vector
        m <- h * avgcZ / avgc
        ## SCMS (update) vector
        dy <- m - Vd %*% colSums(Vd * m)
        ## Stepsize relaxation
        dy <- dy * omega
        ## Stepsize regularization
        if(!is.null(maxStep)) {
            relStep <- sqrt(sum(dy^2)) / maxStep
            if(relStep > 1) dy <- dy / relStep
        }
        yNew <- y + drop(dy)
        ## Cosine of the angle between gradient and the "local normal space".
        cosNormal <- sum(m * dy) / (sqrt(sum(m^2) * sum(dy^2)) + .Machine$double.eps)
        ## Return vd1 at if state is final.
        if (returnV & (cosNormal <= minCos | iterCounter == maxIter))
            V[,,i] <<- spectr$vectors
        c(yNew, cosNormal = cosNormal, lambdaNormal = lambdaNormal, eigengap = eigengap)
    }
    if (!silent) message("Updating points:\n...", M)
    iterCounter <- 1L
    updatedPoints <- integer(maxIter)
    updatedPoints[[iterCounter]] <- M
    Yc <- vapply(seq_len(M), updatePoint, double(n + extraVars))
    cols <- which(Yc["cosNormal",] > minCos)
    while(length(cols) > 0L & iterCounter < maxIter) {
        if (!silent) message("...", length(cols))
        iterCounter <- iterCounter + 1L
        updatedPoints[[iterCounter]] <- length(cols)
        Yc[,cols] <- vapply(cols, updatePoint, double(n + extraVars))
        cols <- which(Yc["cosNormal",] > minCos)
    }
    ret <- list(Y = Yc[seq_len(n),], cosNormal = Yc["cosNormal",],
                lambdaNormal = Yc["lambdaNormal",], eigengap = Yc["eigengap",],
                updatedPoints = updatedPoints)
    if (returnV) ret$V <- V
    ret
}

## `tcrossprod()` takes ~55% time than `. %*% t(.)`.
## A <- matrix(rnorm(40000*200), 200)
## A2 <- Matrix::Matrix(A)
## system.time(mbm <- microbenchmark::microbenchmark(
##                 tcrossprod = tcrossprod(A),
##                 ## tprod = A %*% t(A),
##                 MatrixTCP = Matrix::tcrossprod(A2),
##                 times = 30)) #14.7s
## ggplot2::autoplot(mbm)


#' MParzen KDE-based Subspace Constrained Mean Shift
#'
#' Using the log Manifold Parzen windows (MParzen) Gaussian KDE.
#' @param Y0 Initial points, an n-by-M matrix.
#' @param X Data, an n-by-N matrix.
#' @param model Top d-th order local covariance structure; output of `MParzenTrain()`.
#' @param sigma Maximum noise level added to the normal space of each kernel,
#'        i.e. marginal standard deviation of an isotropic Gaussian. Can be zero.
#' @param h Bandwidth for isotropic Gaussian kernel density estimation, a scalar.
#'          Determines step size.
#' @param r Radius of nearest neighbors to use in evaluating density and derivatives.
#'          Default to use all data.
#' @param minCos  Convergence criterion, cosine of the angle between gradient and
#'                its component in the "local normal space".
#' @param maxIter Maximum number of iterations.
#' @return A list: `Y`, the final points; `cosNormal`, convergence score;
#'                 `updatedPoints`, the number of updated points per update.
scmsMParzen <- function(Y0, X, model, sigma, h, r = NULL, minCos = 0.01, maxIter = 100L) {
    N <- ncol(X)
    M <- ncol(Y0)
    n <- nrow(X)
    stopifnot(nrow(Y0) == n)
    d <- nrow(model$sdev)
    noise <- pmin(model$sdev[d,], sigma)
    invrsd <- sqrt(rep(noise ^ -2, each = d) - model$sdev ^ -2)
    if (is.null(r)) neighbor <- NULL
    else {
        ## Fixed-radius nearest neighbors to initial points: transpose to "data matrices".
        neighbor <- dbscan::frNN(t(X), r, query = t(Y0), sort = FALSE, approx = 0)$id
        nFewNeighbor <- sum(lengths(neighbor) <= 1L)
        if (nFewNeighbor > 0L) warning(nFewNeighbor, " points with only one or no neighbor.")
    }
    ## Prepare coordinate-cosine matrix for updates later.
    Yc <- rbind(Y0, rep(NA_real_, M))
    ## Update a point by index.
    ## Return new coordinates and cosine between gradient and local normal space.
    updatePoint <- function(i) {
        y <- Yc[-nrow(Yc), i]
        idNN <- neighbor[[i]]
        if (is.null(idNN)) yX <- X - y
        else if (length(idNN) > 1L) yX <- X[, idNN] - y
        else if (length(idNN) == 1L) yX <- matrix(X[, idNN] - y, n)
        else return(c(y, cosNormal = 0))
        if (d > 1L) {
            residual <- function(i) {
                invrvVY <- invrsd[,i]^2 * colSums(model$tangent[,,i] * yX[,i])
                model$tangent[,,i] %*% invrvVY
            }
        } else {
            residual <- function(i) {
                model$tangent[,,i] * invrsd[,i]^2 * sum(model$tangent[,,i] * yX[,i])
            }
        }
        if (is.null(idNN)) idNN <- seq_len(N)
        noise <- noise[idNN]
        Z <- yX * rep(noise ^ -2, each = n) - vapply(idNN, residual, double(n))
        ## Kernel density.
        sdev <- model$sdev[,idNN]
        if (d > 1L)
            fc <- exp(-0.5 * colSums(yX * Z) - (n - d) * log(noise) - colSums(log(sdev)))
        else
            fc <- exp(-0.5 * colSums(yX * Z) - (n - d) * log(noise) - log(sdev))
        gc <- drop(Z %*% fc / length(idNN))
        ## Log MParzen KDE, top d eigenvectors of the Hessian.
        avgfc <- mean(fc)
        rootpf <- sqrt(fc / sum(fc))
        rootpfZ <- rep(rootpf, each = n) * Z
        sumpfZ <- gc / avgfc
        rootpfInvrsdV <- rep(rootpf, each = n * d) * rep(invrsd[,idNN], each = n) *
            model$tangent[,,idNN]
        dim(rootpfInvrsdV) <- c(n, d * length(idNN))
        spectr <- eigen(tcrossprod(rootpfZ) + tcrossprod(rootpfInvrsdV) - tcrossprod(sumpfZ),
                        symmetric = TRUE)
        Vd <- spectr$vectors[, seq(d)]
        ## Partial eigen-decomposition
        ## Input to `RSpectra::eigs_sym()` must be a square matrix of size at least 3.
        ## spectr <- RSpectra::eigs_sym(tcrossprod(rootpcZ) - tcrossprod(sumpcZ), d)
        ## Vd <- spectr$vectors
        ## NOTE: special case for d == 1L
        dim(Vd) <- c(n, d)
        ## Mean-shift (update) vector
        m <- gc / avgfc * h^2
        ## SCMS (update) vector
        dy <- m - drop(Vd %*% colSums(Vd * m))
        yNew <- y + dy
        ## Cosine of the angle between gradient and the "local normal space".
        cosNormal <- sum(m * dy) / sqrt(sum(m^2) * sum(dy^2))
        c(yNew, cosNormal = cosNormal)
    }
    message("Updating points: ", M)
    Yc <- vapply(seq_len(M), updatePoint, double(n + 1))
    updatedPoints <- M
    cols <- which(Yc["cosNormal",] > minCos)
    while(length(cols) > 0L & length(updatedPoints) <= maxIter) {
        message("...", updatedPoints[length(updatedPoints)])
        Yc[,cols] <- vapply(cols, updatePoint, double(n + 1))
        updatedPoints[length(updatedPoints) + 1L] <- length(cols)
        cols <- which(Yc["cosNormal",] > minCos)
    }
    list(Y = Yc[-nrow(Yc),], cosNormal = Yc["cosNormal",], updatedPoints = updatedPoints)
}

#' Ridge Estimation
#' @param X0    Initial data.table, N-by-n.
#' @param d     Dimension of density ridge.
#' @param sigma KDE bandwidth for SCMS.
#' @param omega Stepsize relaxation factor, default to 1.
#' @param returnV Whether to return the eigenvector at ridge; an (n,n,N) array.
#' @return      Data.table of estimated density ridge.
#' @export
estimateRidge <- function(X0, d, sigma, minCos = 0.05, omega = 1, maxStep = NULL,
                          maxIter = 100L, returnV = FALSE, silent = FALSE) {
    X0 <- copy(X0)
    N <- nrow(X0)
    n <- ncol(X0)
    R0 <- t(as.matrix(X0))
    X0[, outX := outlier(t(R0))]
    Xs <- t(as.matrix(X0[!(outX), -c("outX")]))
    llog <- scms(R0, Xs, d = d, h = sigma, minCos = minCos,
                 omega = omega, maxStep = maxStep, maxIter = maxIter,
                 returnV = returnV, silent = silent)
    X0[, paste0("r", seq_len(n)) := as.data.table(t(llog$Y))]
    cols <- c("lambdaNormal", "eigengap", "cosNormal")
    X0[, (cols) := llog[cols]]
    if (returnV) return(list(X0 = X0, V = llog$V))
    X0[]
}
