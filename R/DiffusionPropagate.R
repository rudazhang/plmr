## Geometric Diffusion on Manifolds

#' Estimate probability measure on mesh by fine Gaussian kernels
#' @param M Mesh points, an n-by-m' matrix.
#' @param X Sample data, an n-by-m matrix.
#' @param mult Multiples of the mean square distance within the set, determines the scale parameter.
#' @return A list: p, probability mass on mesh; f, probability density on mesh.
pmesh <- function(M, X, mult = 1/2) {
    ## Sample size
    n <- nrow(X)
    m <- ncol(X)
    mp <- ncol(M)
    stopifnot(m > n, mp > nrow(M), nrow(M) == n)
    ## Square distance matrix
    D2 <- as.matrix(Rfast::dista(t(M), t(X), square = TRUE))
    ## Determine scale parameter
    mindist2 <- apply(D2, 1, function(x) min(x[x > 0]))
    eps <- mean(mindist2) * mult^2
    ## Gaussian kernel: K_epsilon(M, X)
    K <- exp(- D2 / eps)
    ## Normalization
    v <- colSums(K)
    ## Approximate probability mass: p_epsilon(M)
    p <- rowSums(K / rep(v, rep(mp, m)))
    ## Approximate probability density: f_epsilon(M)
    f <- rowSums(K)
    list(p = p / sum(p), f = f)
}

#' Heat diffusion of the probability measure/mass on mesh
#'
#' @param M Mesh points, an n-by-m' matrix.
#' @param p Initial probability vector.
#' @param t Diffusion time, as "relative time" for the diffusion to reach the whole set.
#' @param mult Multiples of the mean square distance within the set, determines the scale parameter.
#' @param smooth Multiple of the scale parameter for smoothing density estimate.
#' #' @return Probability masses on the mesh points.
gd <- function(M, p0, t, mult = 3, smooth = 1, scale = NULL) {
    ## Sample size
    n <- nrow(M)
    mp <- ncol(M)
    stopifnot(mp > n)
    ## Distance matrix
    D <- as.matrix(dist(t(M)))
    dimnames(D) <- NULL
    ## Determine scale parameter
    mindist <- apply(D, 1, function(x) min(x[x > 0]))
    eps <- mean((mindist * mult)^2)
    ## Gaussian kernel: K_epsilon
    K <- exp(- D^2 / eps)
    Ks <- exp(- D^2 / (eps * smooth^2))
    ## Approximate density: p_epsilon
    p <- rowSums(K)
    ## Density-neutralizing operator: \tilde{K}_epsilon
    K2 <- K / (p %o% p)
    ## Diffusion matrix: \bar{A}_epsilon
    v <- sqrt(rowSums(K2))
    A <- K2 / (v %o% v)
    ## Eigen-decomposition: decreasing eigenvalues (can be slow)
    AEigen <- eigen(A, symmetric = TRUE) # Slow!
    stopifnot(abs(AEigen$values[[1]] - 1) < 30 * .Machine$double.eps,
              AEigen$values[[mp]] > -20 * .Machine$double.eps) ## Eigenvalues in (0,1].
    ## Transform to eigenvectors of the transition matrix, and normalize them
    PEigenvectors <- AEigen$vectors / v
    PEigenvectorsNorm <- sqrt(colSums(PEigenvectors^2))
    PEigenvectors <- PEigenvectors / rep(PEigenvectorsNorm, rep(mp, mp))
    PEigenvectorsInv <- PEigenvectorsNorm * t(v * AEigen$vectors)
    if (is.null(scale)) scale <- (t * mp) %/% mult
    pt <- as.vector(p0 %*% PEigenvectors %*% (AEigen$values^scale * PEigenvectorsInv))
    ## Peron vector: the stationary distribution (scale -> Inf)
    peron <- PEigenvectorsInv[1,] / sum(PEigenvectorsInv[1,])
    pt <- pt / sum(pt)
    list(pt = pt, ft = as.vector(Ks %*% pt), peron = peron, fPeron = as.vector(Ks %*% peron))
}
