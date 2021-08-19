#' Construct a smooth frame from an initial orthonormal frame
#'
#' A general function that constructs a smooth frame of a vector bundle over a smooth manifold.
#' It can be applied to the normal or tangent bundle of an estimated density ridge.
#' @param V Initial orthonormal frame, an (n, d, N) array for a rank-d bundle.
#' @param Mridge Ridge points, an N-by-n matrix.
#' @param j Index of the reference frame.
#' @references See Supplementary Materials of [@ZhangRD2021nbb].
#' For the `align` sub-procedure, see [@Rheinboldt1988] cited therein.
#' @export
SmoothFrame <- function(V, Mridge, j = NULL) {
    N <- dim(V)[[3]]
    E <- array(rep(NA_real_, length(V)), dim = dim(V))
    q <- rep(NA_real_, N)
    b <- 1L
    if (is.null(j)) j <- round(N / 2)
    E[, , j] <- V[, , j]
    ret <- FNN::get.knn(Mridge, N - 1)
    K <- ret$nn.index
    D <- ret$nn.dist
    k <- rep(1L, N)
    while(b < N) {
        q[[b]] <- j
        K[K == j] <- NA_integer_
        a <- q[seq(b)]
        for (i in a) {
            while(is.na(K[i ,k[[i]]])) {
                k[[i]] <- k[[i]] + 1L
            }
        }
        i <- a[[which.min(D[cbind(a, k[a])])]]
        j <- K[i, k[[i]]]
        ## Align(j, i)
        Theta <- t(V[, , j]) %*% E[, , i]
        ret <- svd(Theta, ncol(Theta))
        Q <- ret$u %*% t(ret$v)
        E[, , j] <- V[, , j] %*% Q
        ## End Align(j, i)
        b <- b + 1L
    }
    return(E)
}

