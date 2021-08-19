## Probabilistic models used in the papers:
## data generating process, helper functions.

#' Data Generating Process for the Parabola Example
#'
#' y = x^2
#' x ~ U[0,1]
#' x' = x + e * (1 + 2 x)
#' y' = y + e * (1 + 2 y)
#' @param N sample size
dgpParabola <- function(N) {
    DT <- data.table::data.table(x0 = runif(N, min = 0, max = 1),
                                 ex = rnorm(N, sd = 0.03),
                                 ey = rnorm(N, sd = 0.03))
    DT[, y0 := x0^2]
    DT[, c("x", "y") := list(x0 + ex * (1 + 2 * x0), y0 + ey * (1 + 2 * y0))]
    DT[, .(x0, y0, ex, ey, x, y)]
}

#' Radial coordinate of the parabola example.
#'
#' The radial partitions the data set rather evenly.
#' @return radial angle of input points, between 0 and pi, mostly around pi/2.
radial <- function(x, y) {
    topleft <- c(0, 1)
    offset <- 0.2
    p0 <- topleft + offset * c(-1, 1)
    p1 <- p0 + c(1, 1)
    d01sq <- sum((p0 - p1)^2)
    d0sq <- (x - p0[1])^2 + (y - p0[2])^2
    d1sq <- (x - p1[1])^2 + (y - p1[2])^2
    cosine <- (d0sq + d01sq - d1sq) / (2 * sqrt(d0sq * d01sq))
    acos(cosine)
}
