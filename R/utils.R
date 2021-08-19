## Internal helper functions

#' Euclidean norm of 2-d vectors
norm2 <- function(x, y) sqrt(x^2 + y^2)

#' Mode estimation in 1d
#'
#' Assumes a unimodal distribution.
#' @param x    Data in 1d, a vector.
#' @param minD Convergence criterion for mean shift.
#' @return     Estimated mode.
mode1d <- function(x, minD, silent = FALSE) {
    ## Use Silverman's rule-of-thumb bandwidth,
    ## assuming distributions in the normal spaces are approximately Gaussian.
    sigma <- silverman(length(x), 1, sd(x))
    X <- matrix(x, ncol = 1L)
    ## Because there is only one mode, to reduce computation,
    ## pick one starting point near both ends of the data.
    Y0 <- matrix(quantile(x, c(0.1, 0.9)), ncol = 1L)
    ## If you care less about accuracy: (not recommended, better reduce minD)
    ## Y0 <- matrix(median(x), ncol = 1L)
    ret <- MeanShift(X, h = sigma, Y0 = Y0, minD = minD, silent = silent)
    if (!silent) cat(".")
    mean(ret$Y)
}

#' Bootstrap statistics of sample distribution of the mean
#' @param x  Sample of a random variable.
#' @param B  Bootstrap size for estimating sampling distribution of the mean, default to sample size.
#' @param seed Random seed for bootstrapping.
#' @return   Data.table of sample mean, standard error, bootstrap range and 90% confidence intervals.
#'           If B = 100, bootstrap range is ~98% CI.
bootStatsMean <- function(x, B = NULL, seed = 42L) {
    N <- length(x)
    if (is.null(B)) B <- N
    set.seed(seed)
    avgs <- replicate(B, mean(sample(x, N, replace = TRUE)))
    data.table(avg = mean(x), stderr = sd(avgs),
               bmin = min(avgs), bmax = max(avgs),
               bp05 = quantile(avgs, 0.05), bp95 = quantile(avgs, 0.95))
}

#' Transform data space
#'
#' From a 2D histogram of grid size 0.02 / 0.03 (smoothing error < ~7.7e-3 / 1.15e-2)
#' @param x     original x coordinate
#' @param y     original y coordinate
#' @param theta rotation angle, clockwise, in radius
#' @param right offset rightward
#' @param up    offset upward
#' @return A list: (xt, yt), transformed x and y coordinates.
affineTransform <- function(x, y, theta, right = 0, up = 0) {
    xr <- x * cos(theta) + y * sin(theta)
    yr <- x * -sin(theta) + y * cos(theta)
    list(xt = xr + right, yt = yr + up)
}

