#' Normal-bundle bootstrap.
#'
#' A direct implementation of NBB without using smooth frame.
#' @param X0    Initial data/sample, a data.table: x, y.
#' @param sigma KDE bandwidth for SCMS dynamical system.
#' @param d     Dimension of density ridge.
#' @param knn   Number of nearest neighbors to use, default to all.
#' @return New data: I, index of ridge point; k, order of projection vector by distance on ridge.
#' @references Algorithm 3.1 in [@ZhangRD2021nbb].
#' @seealso \code{\link{scms}} for ridge estimation;
#' \code{\link{SmoothFrame}} for smooth frame construction.
#' @export
NBB <- function(X0, sigma = NULL, d = 1L, knn = NULL) {
    N <- nrow(X0)
    if (is.null(sigma)) {
        ## KDE bandwidth by maximizing likelihood with a sample of 300.
        sigmaN300 <- 4e-2 / h
        oversmoothing <- 2.5 #2-2.5
        sigma <- sigmaN300 * sqrt(300L / N) * oversmoothing
    }
    ## SCMS
    R0 <- t(as.matrix(X0))
    X0$outX <- outlier(t(R0))
    Xs <- t(as.matrix(X0[!(outX), .(x, y)]))
    llog <- scms(R0, Xs, d = d, h = sigma, minCos = 0.05)
    X0[, c("xr", "yr") := list(llog$Y[1,], llog$Y[2,])]
    X0[, c("dx", "dy") := list(x - xr, y - yr)]
    if (is.null(knn)) knn <- N - 1L
    idNeighbors <- FNN::knn.index(as.matrix(X0[, .(xr, yr)]), knn)
    X0[, I := .I]
    X <- cbind(X0[rep(seq(.N), knn), .(I, xr, yr)],
               X0[as.vector(idNeighbors), .(dx, dy, k = rep(seq(knn), each = N))])
    X[, c("x", "y") := .(xr + dx, yr + dy)]
    X[, .(I, k, x, y)]
}
