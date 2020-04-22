## Manifold Sampling: Gaussian KDE projected on truncated diffusion maps basis.

#' Expectations of some random variables.
#'
#' E|z|, absolute value of standard Gaussian.
#' E(χ_2), Rayleigh distribution, or chi with 2 degrees of freedom, χ(2).
#' @keywords internal
E <- list(absZ = sqrt(2/pi), X2 = sqrt(2) * gamma(3/2))

#' Preprocessing
#' @param X0 The original sample matrix, m-by-n; can be a data.frame.
#' @param scale A logical for whether the variables should be scaled.
#' @return A list: pc, affine transformation to normalized principal component (PC) coordinates;
#'                 H0, normalized sample matrix (zero-mean, unit-variance), m-by-n;
#' @details All PC coordinates are currently assumed to be non-degenerate.
#'          Further generalization is needed.
preprocess <- function(X0, scale = FALSE) {
    pc <- prcomp(X0, scale. = scale)
    H0 <- pc$x / rep(pc$sdev, rep(nrow(X0), ncol(X0)))
    dimnames(H0) <- NULL
    pc$x <- NULL
    list(pc = pc, H0 = H0)
}

#' Manifold learning via diffusion map
#' @param H0 The normalized sample matrix, m-by-n.
#' @param mult Multiples of the mean square distance within the set, determines the scale parameter.
#' @param kappa Analysis scale, as "relative time" for the diffusion to reach the whole set.
#' @param delta Convergence criterion
#' @return A list: lambda, eigenvalues of the diffusion map;
#'                 g, diffusion map basis;
#'                 k, minimal truncation order satisfying the convergence criterion.
manifoldLearn <- function(H0, mult = 3, kappa = 1, delta = 1e-3) {
    ## Sample size
    stopifnot(nrow(H0) > ncol(H0))
    m <- nrow(H0)
    ## Distance matrix (the same for X0 and H0 if no scaling in preprocessing)
    D <- as.matrix(dist(H0))
    dimnames(D) <- NULL
    ## Determine scale parameter
    mindist <- apply(D, 1, function(x) min(x[x > 0]))
    eps <- mean((mindist * mult)^2)
    ## Gaussian kernel: K_epsilon
    K <- exp(- D^2 / eps)
    ## Approximate density: p_epsilon
    p <- rowSums(K)
    ## Density-neutralizing operator: \tilde{K}_epsilon
    K2 <- K / (p %o% p)
    ## Diffusion matrix: \bar{A}_epsilon
    v <- sqrt(rowSums(K2))
    A <- K2 / (v %o% v)
    ## Eigen-decomposition: decreasing eigenvalues (can be slow)
    AEigen <- eigen(A, symmetric = TRUE)
    stopifnot(abs(AEigen$values[[1]] - 1) < 5 * .Machine$double.eps,
              AEigen$values[[m]] > -.Machine$double.eps) ## Eigenvalues in (0,1].
    ## Transform to eigenvectors of the transition matrix, and normalize them
    PEigenvectors <- AEigen$vectors / v
    PEigenvectors <- PEigenvectors / rep(sqrt(colSums(PEigenvectors^2)), rep(m, m))
    ## Diffusion map basis
    g <- PEigenvectors * rep(AEigen$values^((kappa * m) %/% mult), rep(m, m))
    ## Convergence criterion for truncation order
    c0 <- H0 %*% t(H0) / (m - 1)
    k <- 1L
    repeat {
        Hm <- manifoldSample(H0, g, k)
        cm <- Hm %*% t(Hm) / (m - 1)
        if (norm(cm - c0, 'f') / norm(c0, 'f') < eps0) break
        k <- k + 1L
    }
    list(lambda = AEigen$values, g = g, k = k, vectors = PEigenvectors)
}

#' Manifold sampling: Gaussian KDE projected on truncated diffusion maps basis.
#' @param H0 The normalized sample matrix, n-by-m.
#' @param g Diffusion map basis
#' @param m Truncation order
#' @param mult Multiple in sample size of the new sample; if 0, project the original sample.
#' @param avgScatter Average relative scattering distance from the original sample.
#' @param seed Seed for random number generation.
#' @return A list: kde, KDE sample matrix, n-by-m x mult;
#'                 kdem, manifold-reduced sample matrix, n-by-m x mult.
manifoldSample <- function(H0, g, m, mult = 0L, avgScatter = 0.1, seed = 42L) {
    ## Sample size
    n <- nrow(H0)
    m <- ncol(H0)
    stopifnot(n < m)
    ## Truncation
    gm <- g[, seq(m)]
    ## Projection
    Pgm <- gm %*% chol2inv(chol(t(gm) %*% gm)) %*% t(gm)
    if (mult == 0L) return(H0 %*% Pgm)
    set.seed(seed)
    ## Retaining order
    H00 <- H0[, rep(seq.int(m), mult)]
    kde <- H00 + rnorm(length(H00)) * avgScatter / E$X2
    dim(kde) <- c(n, m, mult)
    kdem <- apply(kde, 3, `%*%`, Pgm)
    dim(kdem) <- c(n, m * mult)
    dim(kde) <- c(n, m * mult)
    list(kde = kde, kdem = kdem)
}

#' Reverse prediction: principal components to original coordinates
#' @param H Principal component coordinates, m-by-n.
#' @param pc Affine transformation to principal component coordinates.
#' @param scale Whether restore original scale.
#' @return Original coordinates, m-by-n.
revpredict <- function(H, pc) {
    ## Use `rep()` for better performance and readability.
    ## https://stackoverflow.com/a/32364355/5407633
    reps <- rep(nrow(H), ncol(H))
    if(is.null(pc$scale))
        return((H * rep(pc$sdev, reps)) %*% t(pc$rotation) + rep(pc$center, reps))
    (H * rep(pc$sdev, reps)) %*% t(pc$rotation) * rep(pc$scale, reps) + rep(pc$center, reps)
}
