## Parametrized submanifolds (d = 1, 2) and sampling on them.

#' Spherical to Cartesian coordinates
spherical2cartesian <- function(theta, phi) {
    list(x = sin(pi * phi) * cos(2 * pi * theta),
         y = sin(pi * phi) * sin(2 * pi * theta),
         z = cos(pi * phi))
}

#' Cartesian to Spherical coordinates
cartesian2spherical <- function(x, y, z) {
    r12 <- sqrt(x^2 + y^2)
    r <- sqrt(x^2 + y^2 + z^2)
    list(r = r,
         theta = ifelse(y > 0, acos(x / r12), 2 * pi - acos(x / r12)),
         phi = acos(z / r))
}

#' Helix 1: cylindrical coordinates
#' @param N Sample size.
#' @param invCDF Inverse CDF, for probability integral transformation.
#' @param a Major radius.
#' @param b Minor radius.
#' @param n Number of periods.
rhelix1 <- function(N, invCDF, a = 2, b = 1, n = 7) {
    DT <- data.table::data.table(s = seq(N)/N)
    DT[, t := vapply(s, invCDF, double(1))]
    DT[, theta := n * 2 * pi * t]
    DT[, r := a + b * cos(2 * pi * t)]
    DT[, z := 3 * b * sin(2 * pi * t)]
    DT[, x := r * cos(theta)]
    DT[, y := r * sin(theta)]
    DT[]
}

#' Helix 2: cylindrical coordinates
rhelix2 <- function(N, invCDF, a = 2, b = 1, n = 8) {
    DT <- data.table::data.table(s = seq(N)/N)
    DT[, t := vapply(s, invCDF, double(1))]
    DT[, theta := 2 * pi * t]
    DT[, r := a + b * cos(n * 2 * pi * t)]
    DT[, z := b * sin(n * 2 * pi * t)]
    DT[, x := r * cos(theta)]
    DT[, y := r * sin(theta)]
    DT[]
}

#' Trefoil curve
rtrefoil <- function(N, invCDF, a = 1, b = 2, n = 3) {
    DT <- data.table::data.table(s = seq(N)/N)
    DT[, t := vapply(s, invCDF, double(1))]
    DT[, x := sin(2 * pi * t) + b * sin(b * 2 * pi * t)]
    DT[, y := cos(2 * pi * t) - b * cos(b * 2 * pi * t)]
    DT[, z := sin(n * 2 * pi * t)]
    DT[]
}

## Density along parameter t
## pdft <- function(t) 5/3 * cos(pi * t)^2 + 1/6
## Distribution along parameter t
## cdft <- function(t) sin(2*pi*t) * 5 / (12*pi) + t
## Probability integral transformation
## invCDF <- function(u) uniroot(function(t) cdft(t) - u, c(0, 1))$root
## Visualization
## DT[, rgl::plot3d(x, y, z, aspect = FALSE, col = hsv(h = t))]

#' Sampling uniformly on the unit circle.
#' @param N Sample size.
#' @return  Data.table: theta, azimuth angle; x, y, Cartesian coordinates.
runifcircle <- function(N) {
    DT <- data.table::data.table(theta = runif(N) * 2 * pi)
    DT[, c("x", "y") := .(cos(theta), sin(theta))]
    DT[]
}

#' Sampling on a circular Gaussian distribution.
#' @param N     Sample size.
#' @param sigma Standard deviation of Gaussian noise in the radial direction.
#' @return Data.table: theta, azimuth angle; r, radial coordinate; x, y, Cartesian coordinates.
rCircularGauss <- function(N, sigma) {
    X <- data.table(theta = runif(N) * 2 * pi, r = 1 + rnorm(N) * sigma)
    X[, c("x", "y") := .(r * cos(theta), r * sin(theta))]
    X[]
}

#' Sampling on a circular chi distribution in a flow model.
#'
#' Circular chi: uniform distribution in azimuth angle,
#' mode-normalized chi distribution in the radial coordinate.
#' Density ridge is r=1.
#' @param N  Sample size.
#' @param df Degree of freedom of the chi distribution; df >= 2, integer.
#' @return   Data.table: theta, azimuth angle; r, radial coordinate; x, y, Cartesian coordinates.
rCircularChi <- function(N, df) {
    DT <- data.table::data.table(theta = runif(N) * 2 * pi,
                                 r = sqrt(rchisq(N, df = df) / (df - 1)))
    DT[, c("x", "y") := .(r * cos(theta), r * sin(theta))]
    DT[]
}

#' Log density of the unit circle with Gaussian noise, disregarding the constant term.
#' @param r radial coordinate
#' @param h bandwidth of Gaussian noise
logPNoisyCircle <- function(r, h = 0.2) {
    scaledBesselI0 <- besselI(r / h^2, 0, expon.scaled = TRUE)
    logBesselI0 <- log(scaledBesselI0) + r / h^2
    logBesselI0 - r^2 / (2 * h^2)
}

#' Partial derivative of the log density in the radial direction.
drLogPNoisyCircle <- function(r, h = 0.2) {
    scaledBesselI1 <- besselI(r / h^2, 1, expon.scaled = TRUE)
    scaledBesselI0 <- besselI(r / h^2, 0, expon.scaled = TRUE)
    (scaledBesselI1 / scaledBesselI0 - r) / h^2
}


#' Sampling uniformly on the unit sphere, as a flow model.
#' @param N Sample size.
#' @return  Data.table: theta, azimuth angle; phi, polar angle.
runifsphere <- function(N) {
    data.table::data.table(theta = runif(N) * 2 * pi, phi = acos(1 - 2 * runif(N)))
}

#' Sampling from a distribution on the unit sphere
#' @param N Sample size.
#' @param f Probability density as a function of spherical coordinates (theta, phi).
rsphere <- function(N, f) {
    maxf <- -optim(c(1, 1), function(x) -f(x[[1]], x[[2]]))$value
    DT <- data.table::data.table(theta = NA, phi = NA, accept = FALSE)
    while (sum(DT$accept) < N) {
        DT2 <- data.table::data.table(theta = runif(N), phi = runif(N))
        DT2[, f := f(theta, phi)]
        DT2[, accept := runif(.N) < sin(pi * phi) * f / maxf]
        DT <- rbind(DT, DT2[, -"f"])
    }
    DT[(accept)][-sample(.N, .N - N, replace = FALSE), .(theta, phi)]
}

#' Sampling from a distribution on a torus
#' @param N Sample size.
#' @param f Probability density as a function of spherical coordinates (theta, phi).
#' @param a Radius of big circle.
#' @param b Radius of small circle.
rtorus <- function(N, f, a, b) {
    maxf <- -optim(c(1, 1), function(x) -f(x[[1]], x[[2]]))$value
    DT <- data.table::data.table(theta = NA, phi = NA, accept = FALSE)
    while (sum(DT$accept) < N) {
        DT2 <- data.table::data.table(theta = runif(N), phi = runif(N))
        DT2[, f := f(theta, phi)]
        DT2[, accept := runif(.N) < (a + b * sin(2 * pi * phi)) / (a + b) * f / maxf]
        DT <- rbind(DT, DT2[, -"f"])
    }
    DT[(accept)][-sample(.N, .N - N, replace = FALSE), .(theta, phi)]
}

## N <- 4000
## a <- 1
## b <- 0.9
## ## Sphere: Uniform on parameter space
## DT <- data.table::data.table(theta = runif(N), phi = runif(N))
## ## Sphere: Uniform on manifold
## DT <- rsphere(N, function(theta, phi) 1)
## ## Sphere: sin(phi) on manifold
## DT <- rsphere(N, function(theta, phi) sin(pi * phi))
## ## Plots: parameter space vs. embedding space.
## DT[, plot(theta, phi, pch = 19, col = gray(0, .3))]
## rect(0,0,1,1, col = gray(0, .2), border = FALSE)
## DT[, rgl::plot3d(sin(pi * phi) * cos(2 * pi * theta),
##                  sin(pi * phi) * sin(2 * pi * theta),
##                  cos(pi * phi), alpha = .3, aspect = FALSE)]
## ## Torus: Uniform on parameter space
## DT <- data.table::data.table(theta = runif(N), phi = runif(N))
## ## Torus: Uniform on manifold
## DT <- rtorus(N, function(theta, phi) 1, a, b)
## ## Plots: embedding space.
## DT[, rgl::plot3d((a + b * sin(2 * pi * phi)) * cos(2 * pi * theta),
##                  (a + b * sin(2 * pi * phi)) * sin(2 * pi * theta),
##                  b * cos(2 * pi * phi), alpha = .3, aspect = FALSE)]
