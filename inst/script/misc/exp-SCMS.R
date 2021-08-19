## Examples for the graphical interpretation of SCMS.
library(data.table)
library(purrr)
library(rgl)

## ------------------------------------------------------------------------------------------
## Ridge estimation on small data
## When data size is small and oversmothing is ~1, ridge estimation approximates (piecewise) PCA.
## Circle in the plane:
##     2, always line; 3, pinched triangle to three rays; 4, always square;
##     5 slowly smoothing; 6+, polygon to circle;
nSparse <- 2L
nDense <- 180L
sigmaN300 <- 4e-2
oversmoothing <- 1 #1; 1.25; 2; 2.5
rad <- seq(0, nDense - 1) / nDense * 2 * pi
DTcircle <- data.table(x = cos(rad), y = sin(rad))
DTcircle[, sub := FALSE]
DTcircle[seq(1, .N, by = .N / nSparse), sub := TRUE]
R0 <- t(as.matrix(DTcircle[, .(x, y)]))
Xs <- t(as.matrix(DTcircle[(sub), .(x, y)]))
sigma <- sigmaN300 * sqrt(300L / nSparse) * oversmoothing
llog <- scms(R0, Xs, d = 1, h = sigma, minCos = 0.05)
DTcircle[, c("xr", "yr") := list(llog$Y[1,], llog$Y[2,])]
DTcircle[, plot(x, y, asp = 1)]
DTcircle[(sub), points(x, y, pch = 19)]
DTcircle[, points(xr, yr, pch = 19, col = "gold")]
## Try other oversmoothing values
oversmoothing <- 2 #1.25;2.5
sigma <- sigmaN300 * sqrt(300L / nSparse) * oversmoothing
llog <- scms(R0, Xs, d = 1, h = sigma, minCos = 0.05)
DTcircle[, c("xr", "yr") := list(llog$Y[1,], llog$Y[2,])]
DTcircle[, points(xr, yr, pch = 19, col = "cyan")]

## Sphere in the 3-space
##     6, increasingly pinched smoothed octahedron; 8, always cube;
##     12, icosahedron to sphere; 20, dodecahedron (rather stable) to sphere;
library(rgl)
nDenseLine <- 30L
sigmaN300 <- 4e-2
oversmoothing <- 1 #1; 1.25; 2; 2.5
coord <- seq(0, nDenseLine) / nDenseLine * 2 - 1
DTsquare <- as.data.table(expand.grid(x = coord, y = coord))
DTcube <- rbind(DTsquare[, .(x, y, z = -1)], DTsquare[, .(x, y, z = 1)],
                DTsquare[, .(x, y = -1, z = y)], DTsquare[, .(x, y = 1, z = y)],
                DTsquare[, .(x = -1, y, z = x)], DTsquare[, .(x = 1, y, z = x)])
DTcube <- unique(DTcube)
nDense <- DTcube[, .N]
DTcube[, sub := FALSE]
## DTcube[x^2 + y^2 + z^2 == 1, sub := TRUE]
## DTcube[abs(x) + abs(y) + abs(z) == 3, sub := TRUE]
## nSparse <- DTcube[(sub), .N]
goldenRatio <- (sqrt(5) - 1)/2
DTicosahedron <- as.data.table(expand.grid(x = c(-1, 1) * goldenRatio, y = c(-1, 1)))
## DTicosahedron[, len := sqrt(x^2 + y^2)]
## DTicosahedron <- DTicosahedron[, .(x = x / len, y = y / len)]
DTicosahedron <- rbind(DTicosahedron[, .(x, y, z = 0)],
                       DTicosahedron[, .(x = 0, y = x, z = y)],
                       DTicosahedron[, .(x = y, y = 0, z = x)])
## nSparse <- DTicosahedron[, .N]
DTdodecahedron <- as.data.table(expand.grid(x = c(-1, 1) / goldenRatio, y = c(-1, 1) * goldenRatio))
DTdodecahedron <- rbind(DTdodecahedron[, .(x, y, z = 0)],
                        DTdodecahedron[, .(x = 0, y = x, z = y)],
                        DTdodecahedron[, .(x = y, y = 0, z = x)],
                        expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)))
nSparse <- DTdodecahedron[, .N]
R0 <- t(as.matrix(DTcube[, .(x, y, z)]))
## Xs <- t(as.matrix(DTcube[(sub), .(x, y, z)]))
## Xs <- t(as.matrix(DTicosahedron[, .(x, y, z)]))
Xs <- t(as.matrix(DTdodecahedron[, .(x, y, z)]))
sigma <- sigmaN300 * sqrt(300L / nSparse) * oversmoothing
system.time(llog <- scms(R0, Xs, d = 2, h = sigma, minCos = 0.05)) #4s; 1.7; 3.6s
DTcube[, c("xr", "yr", "zr") := list(llog$Y[1,], llog$Y[2,], llog$Y[3,])]
DTcube[, rgl::plot3d(x, y, z, asp = 1, type = 'n', box = FALSE, axes = FALSE,
                     xlab = '', ylab = '', zlab = '')]
## DTcube[(sub), rgl::plot3d(x, y, z, size = 10, add = TRUE)]
## DTicosahedron[, rgl::plot3d(x, y, z, size = 10, add = TRUE)]
DTdodecahedron[, rgl::plot3d(x, y, z, size = 10, add = TRUE)]
DTcube[, rgl::plot3d(xr, yr, zr, col = "red", add = TRUE)]
## Try other oversmoothing values
oversmoothing <- 6
sigma <- sigmaN300 * sqrt(300L / nSparse) * oversmoothing
system.time(llog <- scms(R0, Xs, d = 2, h = sigma, minCos = 0.05)) #43s; 4s
DTcube[, c("xr", "yr", "zr") := list(llog$Y[1,], llog$Y[2,], llog$Y[3,])]
DTcube[, rgl::plot3d(xr, yr, zr, col = "orange", add = TRUE)]

## Check stability of critical points
deviation <- 0.06
set.seed(42L)
DTnoise <- data.table(dx = rnorm(nDense), dy = rnorm(nDense), dz = rnorm(nDense))
cols <- c("dx", "dy", "dz")
DTnoise[, len := sqrt(dx^2 + dy^2 + dz^2)]
DTnoise[, (cols) := lapply(.SD, `*`, deviation / len), .SDcols = cols]
DTnoisyridge <- cbind(DTcube, DTnoise)[, .(x = xr + dx, y = yr + dy, z = zr + dz)]
## DTnoisyridge[, rgl::plot3d(x, y, z, add = TRUE)]
R0 <- t(as.matrix(DTnoisyridge[, .(x, y, z)]))
system.time(llog <- scms(R0, Xs, d = 1, h = sigma, minCos = 0.05)) #3s
DTnoisyridge[, c("xr", "yr", "zr") := list(llog$Y[1,], llog$Y[2,], llog$Y[3,])]
DTnoisyridge[, rgl::plot3d(xr, yr, zr, col = "darkred", add = TRUE)]


## von Mises distribution with small Gaussian noise----------------------------------------
## DTbessel <- data.table(x = seq(0, 1000))
## DTbessel[, K0scale := besselI(x, 0, expon.scaled=TRUE)]
## DTbessel[, K1scale := besselI(x, 1, expon.scaled=TRUE)]

#' log density
LogP <- function(r, theta, eps = 0.1, kappa = 1, flip = FALSE) {
    a <- sqrt(kappa^2 + 2 * kappa * r / eps^2 * cos(theta) + r^2 / eps^4)
    logP <- - (r^2 + 1) / (2 * eps^2) + log(besselI(a, 0, expon.scaled=TRUE)) + a
    if (flip) return(-logP)
    logP
}

#' Determine radial derivative of density (NOT RIDGE!)
f <- function(r, theta, eps = 0.04) {
    a <- sqrt(1 + 2 * r / eps^2 * cos(theta) + r^2 / eps^4)
    (cos(theta) + r / eps^2) * besselI(a, 1, expon.scaled=TRUE) -
        r * a * besselI(a, 0, expon.scaled=TRUE)
}

#' Determine radial derivative of log density (NOT RIDGE!)
g <- function(r, theta, eps = 0.1, kappa = 1) {
    a <- sqrt(1 + 2 * r / eps^2 * cos(theta) + r^2 / eps^4)
    besselI(a, 1, expon.scaled=TRUE) / besselI(a, 0, expon.scaled=TRUE) -
        r * sqrt(1 + (kappa * sin(theta) / (r / eps^2 + kappa * cos(theta)))^2)
}

## Density plot
coord <- seq(-1.5, 1.5, by = 0.01)
DTplot <- data.table::CJ(x = coord, y = coord)
DTplot[, r := sqrt(x^2 + y^2)]
DTplot[, theta := Arg(x + y * 1i)]
DTplot[, logP := LogP(r, theta, eps = 0.4, kappa = 2)]
rgl::open3d()
rgl::material3d(col = "black")
DTplot[, rgl::persp3d(coord, coord, (logP), col = "lightblue")]

## Ridge by root finding
DT <- data.table::data.table(theta = seq(0, 2*pi, length.out=120))

cols <- c("root", "f.root", "iter", "init.it", "estim.prec")
DT[, (cols) := purrr::map_dfr(theta, ~uniroot(purrr::partial(g, theta = ., eps = 0.4, kappa = 2),
                                              c(0.5, 1.3)))]
DT[, plot(theta, 1 - root, type = 'l')]
DT[, x := root * cos(theta)]
DT[, y := root * sin(theta)]
DT[, plot(cos(theta), sin(theta), type = 'l', asp = 1)]
DT[, lines(x, y, col = "red")]

## Ridge by optimization
cols <- c("minimum", "objective")
DT[, (cols) := purrr::map_dfr(theta, ~optimize(purrr::partial(LogP, theta = ., eps = 0.4, kappa = 2,
                                                              flip = TRUE), c(0.5, 1)))]
DT[, plot(theta, 1 - minimum, type = 'l')]
DT[, x := minimum * cos(theta)]
DT[, y := minimum * sin(theta)]
DT[, plot(cos(theta), sin(theta), type = 'l', bty = 'n', asp = 1)]
DT[, lines(x, y, col = "red")]
