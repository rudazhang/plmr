## Examples of Probability Estimation on Manifold via Diffusion
source("parameterized-manifold.R")
source("diffusion.R")
source("ridge.R")
source("data.R")

library(data.table)
library(rgl)

## Circle with doubled density on the upper half -----------------------------------------------
## Data generation
N <- c(train = 2000L, validate = 2000L)
DT <- data.table(type = rep(factor(names(N)), times = N))
avgRelErr <- 0.03
set.seed(420L)
DT[, s := runif(.N)]
DT[, relErr := rnorm(.N) / E$absZ * avgRelErr]
DT[, x := cos(3 * pi * s) * (1 + relErr)]
DT[, y := sin(3 * pi * s) * (1 + relErr)]
DTcircle <- data.table::data.table(s = seq(0, 300) / 300)
DTcircle[, x := cos(2 * pi * s)]
DTcircle[, y := sin(2 * pi * s)]

## Optimal Parzen window
X <- as.matrix(DT[type == "train", .(x, y)])
Y <- as.matrix(DT[type == "validate", .(x, y)])
DTparzen <- data.table(sigma = avgRelErr * seq(0.2, 10, by = 0.2))
DTparzen[, anll := purrr::map_dbl(sigma, ParzenANLL, Y = Y, X = X)]
## DTparzen[anll < 0.3, plot(sigma, anll, type = 'b')]
## DTparzen[, abline(h = min(anll), col = "red")]
optimize(ParzenANLL, avgRelErr * c(0.2, 10), Y = Y, X = X) # sigma = 0.025; ANLL = -0.003340218
## Compare with Silverman's Rule of Thumb for 1d Gaussian distributions.
h <- silverman(nrow(X), ncol(X), apply(X, 2, sd))
cat("Oversmoothing factors:", h / DTparzen[which.min(anll), sigma], "\n") # 8.3, 8.0

## Ridge finding
X <- t(X)
outX <- outlier(X)
Xs <- X[, !outX]
Xs <- X[, -outlier(X)]
system.time(llog <- scms(Xs, Xs, d = 1, h = 0.025, minCos = 0.05)) #(2k, 6s)

M <- l$M[-nrow(l$M), ]
## M <- M[, -outlier(M)]
system.time(l1 <- pmesh(M, X, mult = 1)) #(4k, 0.7s)
time <- 0.1
system.time(l2 <- gd(M, l1$p, t = time, mult = 60, smooth = 4)) #(4k, 31s)
DTmesh <- data.table::data.table(x = M[1, ], y = M[2, ], pEps = l1$p, fEps = l1$f,
                                 pt = l2$pt, ft = l2$ft, pInf = l2$peron, fInf = l2$fPeron)
DTmesh[, s := ifelse(y > 0, acos(x / sqrt(x^2 + y^2)), 2 * pi - acos(x / sqrt(x^2 + y^2))) / (2*pi)]

## Visualizations
barplot(l$updatedPoints, names.arg = seq_along(l$updatedPoints),
        xlab = "Iteration", ylab = "Updated pionts", main = "Ridge Estimation")
## Mesh plot
svglite::svglite("circle-scms.svg", width = 10, height = 10)
DT[type == "train", plot(x, y, asp = 1, pch = 19, col = adjustcolor("red", .2),
                         main = "Circle with doubled density on the upper half")]
DT[type == "validate", points(x, y, pch = 19, col = adjustcolor("blue", .2))]
## DT[, rgl::plot3d(x, y, z, asp = FALSE)]
## DTwire[, rgl::lines3d(x, y, z, col = "red", lwd = 3)]
points(M[1,], M[2,], pch = 19, col = adjustcolor("blue", .3))
segments(DT$x, DT$y, l$M[1,], l$M[2,], col = gray(0, .5))
DTcircle[, lines(x, y, type = 'l', lwd = 2, col = "white")]
abline(h = 0, lty = 2)
legend("bottomright", legend = c("data", "mesh", "manifold"),
       col = c("red", "blue", gray(.8)), pch = 19)
dev.off()
## Viz: CDF
svglite::svglite("circle-CDF.svg", width = 10, height = 10)
DTmesh[order(s)][, plot(s, cumsum(pInf), type = 's', lwd = 3, col = "red",
                        xlim = c(0,1), ylim = c(0,1), xaxs = 'i', yaxs = 'i',
                        main = "CDF along arc length")]
segments(0, 0, 1/2, 2/3); segments(1, 1, 1/2, 2/3)
DTmesh[order(s)][, segments(s[1], pInf[1], s[.N], 1)]
DTmesh[order(s)][, lines(s, .I/.N, type = 's', lwd = 1)]
DTmesh[order(s)][, lines(s, cumsum(pEps), type = 's', lwd = 1, col = "magenta")]
DTmesh[order(s)][, lines(s, cumsum(pt), type = 's', lwd = 2, col = "blue")]
legend("bottomright", legend = paste0("t = ", c("epsilon", time, "Inf")),
       col = c("magenta", "blue", "red"), lty = 1, lwd = 2)
dev.off()
## Viz: PDF
svglite::svglite("circle-PDF.svg", width = 10, height = 10)
DTmesh[, plot(s[order(s)], fEps[order(s)] / max(fEps), type = 'l', lwd = 1, col = "magenta",
              bty = 'n', main = "PDF along arc length (smoothed, normalized)")]
DTmesh[, lines(s[order(s)], fInf[order(s)] / max(fInf[order(s)]), lwd = 3, col = "red")]
DTmesh[, lines(s[order(s)], ft[order(s)] / max(ft[order(s)]), lwd = 3, col = "blue")]
abline(v = c(0, 1/2, 1), h = c(0, 1/2, 1), lty = 2)
legend("bottomright", legend = paste0("t = ", c("epsilon", time, "Inf")),
       col = c("magenta", "blue", "red"), lty = 1, lwd = 2)
dev.off()


## Sphere with density sin(phi) -----------------------------------------------------------------
## Data generation
N <- 4000L
set.seed(420L)
## Sphere: sin(phi) on manifold
DT <- rsphere(N, function(theta, phi) sin(pi * phi))
DT[, c("x", "y", "z") := spherical2cartesian(theta, phi)]
## Plots: parameter space vs. embedding space.
DT[, plot(theta, phi, pch = 19, col = gray(0, .3))]
rect(0,0,1,1, col = gray(0, .2), border = FALSE)
DT[, rgl::plot3d(x, y, z, aspect = FALSE,
                 col = colorRampPalette(c("red", "blue"))(100)[100 * phi + 1], alpha = .5)]
## Ambient noise
avgRelErr <- 0.03
c <- avgRelErr / E$absZ / sqrt(3 - 1)
set.seed(420L)
DT[, c("x", "y", "z") := lapply(.SD, function(x) x + c * rnorm(.N)), .SDcols = c("x", "y", "z")]
## Data plot
DT[, rgl::plot3d(x, y, z, col = ifelse(x^2+y^2+z^2 > 1, "red", "blue"), alpha = .5, aspect = FALSE)]
rgl::view3d(theta = 0, phi = 0, zoom = 0.6)
rgl::rgl.postscript("sphere.svg", "svg")
## 3D data plot
## htmlwidgets::saveWidget(rgl::rglwidget(width = 1000, height = 1000), "sphere-data.html")

## Computation: manifold estimation
X <- t(as.matrix(DT[, .(x, y, z)]))
## (1) Data as initial mesh
## if (length(ol <- outlier(X, 6))) Xs <- X[, -ol]
system.time(l <- scms(X, X, d = 2, epsilon = 0.2, maxIter = 30)) #(4k, 17s)
M <- l$M[-nrow(l$M), ]
(iqr <- diff(quantile(sqrt(colSums(M[seq(3),]^2)), c(1/4, 3/4))))
(mbias <- median(sqrt(colSums(M[seq(3),]^2))) - 1)
## (IQR, median bias): (5.77e-3, 2.73%)
system.time(l <- scms(X, X, d = 2, h = rep(4 * c, 3), epsilon = 0.2, maxIter = 30)) #(4k, 20s)
## (mult, IQR, median bias): (5, 5.77e-3, 1.65%); (4, 6.76e-3, 1.0%); (3, 9.2e-3, 5.0e-3)
r <- sqrt(colSums(M^2))
d <- mindist(M)
plot(d, r, xlab = "minimum distance to set", ylab = "radius")
abline(h = mean(r) + sd(r) * c(-3, 3, 4), v = mean(d) + sd(d) * c(-3, 3))
## Remove outliers
M <- M[, -which(abs((r - mean(r)) / sd(r)) > 3)]
## (2) Mesh with uniform distribution on sphere
DTmesh <- rsphere(N, function(theta, phi) 1)
DTmesh[, c("x", "y", "z") := spherical2cartesian(theta, phi)]
M <- t(as.matrix(DTmesh[, .(x, y, z)]))
## (3) Regular mesh on sphere
## TODO...

## Viz: 3D mesh plot
rgl::plot3d(M[1,], M[2,], M[3,], col = "green", alpha = .3)
## htmlwidgets::saveWidget(rgl::rglwidget(width = 1000, height = 1000), "sphere-mesh.html")
## Viz: mesh radius
svglite::svglite("sphere-radius-ECDF.svg", width = 8, height = 8)
DT[, .(r = sqrt(x^2+y^2+z^2))
   ][order(r)][, plot.ecdf(r, col = ifelse(r < 1, "blue", "red"), lwd = 2, bty = 'n',
                 xlab = "radius", main = "ECDF of radius")]
plot.ecdf(sqrt(colSums(M[seq(3),]^2)), add = TRUE, col = "green", lwd = 2)
abline(v = 1, lty = 2)
legend("bottomright", legend = c("Data (outside sphere)", "Data (inside sphere)", "Mesh"),
       col = c("red", "blue", "green"), lty = 1, lwd = 2)
dev.off()
## Relative standard deviation
(sd(sqrt(colSums(M[seq(3),]^2))) / DT[, sd(sqrt(x^2+y^2+z^2))]) #0.19

## Computation: intial probability estimation
system.time(l1 <- pmesh(M, X, mult = 1/2)) #(4k, 1s)
## Computation: probability diffusion
## (mult = 3, scale = 1); (mult = 2, scale = 10); smooth = 4
mult <- 2L
scale <- 10L
smooth <- 3L
mp <- ncol(M)
system.time(l2 <- gd(M, l1$p, t = NULL, mult = mult, smooth = smooth, scale = scale)) #(4k, 30s)
DTmesh <- data.table::data.table(x = M[1, ], y = M[2, ], z = M[3, ], p0 = l1$p, f0 = l1$f,
                                 pt = l2$pt, ft = l2$ft, pInf = l2$peron, fInf = l2$fPeron)
DTmesh[, c("theta", "phi") := cartesian2spherical(x, y, z)[c("theta", "phi")]]

## Viz: PDF
rgl::par3d(userMatrix = diag(rep(1,4))[c(2,3,1,4),] + 1e-4)
DTmesh[, rgl::plot3d(theta, phi, fInf, alpha = .3)]
rgl::open3d()
DTmesh[, rgl::plot3d(theta, phi, ft, alpha = .3)]
## Viz: Comparison
svglite::svglite("sphere-PDF.svg", width = 10, height = 10)
DTmesh[theta > 4, plot(phi, ft, type = 'n',
                       ylim = c(0, max(ft)*1.02), yaxs = 'i', bty = 'n',
                       main = "PDF along zenith angle (smoothed)")]
DTmesh[theta > 4, points(phi, f0 / max(f0) * max(ft), pch = 19, col = gray(0, 0.2))]
DTmesh[theta > 4, points(phi, ft, pch = 19, col = adjustcolor("blue", 0.2))]
DTmesh[theta > 4, points(phi, fInf, pch = 19, col = adjustcolor("red", 0.2))]
curve(DTmesh[theta > 4 & phi > 0.4 * pi & phi < 0.6 * pi, mean(ft)] * sin(x),
      0, pi, lwd = 2, col = "blue", add = TRUE)
abline(h = DTmesh[theta > 4, median(fInf)], lwd = 2, col = "red")
legend("topright", legend = paste0("t = ", c(0, round(scale * mult / mp, digits = 3), "Inf")),
       col = c("gray", "blue", "red"), pch = 19)
dev.off()


## Ghanem's two-circle data ------------------------------------------------------------
DT <- data.table::fread("two-circles.csv")
DT[, plot(x, y, asp = 1, pch = 19, col = gray(0, 0.3))]

## Density Ridge
X <- t(as.matrix(DT[, .(x, y)]))
M0 <- X
outliers <- outlier(X)
Xs <- if(length(outliers) == 0) X else X[, -outliers]
## Convergence criterion
system.time(l <- scms(M0, Xs, d = 1)) #(231, 2.5s)
system.time(l <- scms(M0, Xs, d = 1, epsilon = 0.055)) #(231, 1.8s)
system.time(l <- scms(M0, Xs, d = 1, epsilon = 0.45)) #(231, 0.083s)
M <- l$M[-nrow(l$M),]
points(t(M), pch = 19, col = "blue")

## Density kernel bandwidth
## Reducing on Silverman's rule of thumb: h is too large; h * 0.1 is too small.
h <- silverman(ncol(X), nrow(X), apply(X, 1, sd))
system.time(l <- scms(M0, Xs, d = 1, h = h * 0.2)) #(231, 2.5s)
M <- l$M[-nrow(l$M), ]
points(t(M), pch = 19, col = "blue")
symbols(rep(0, 2), rep(0, 2), circles = h, inches = FALSE, add = TRUE)
symbols(rep(0, 2), rep(0, 2), circles = h * 0.2, inches = FALSE, add = TRUE)

## Sampling
set.seed(42L)
mErr <- matrix(rnorm(length(M0)), nrow = 2) * h * 0.1
M1 <- M0 + mErr
system.time(l <- scms(M1, Xs, d = 1, h = h * 0.2)) #(231, 2.5s)
M <- l$M[-nrow(l$M), ]

## Plot
svglite::svglite("twocircles-mssm.svg", width = 24, height = 12)
DT[, plot(x, y, asp = 1, pch = 19, col = gray(0, 0.3))]
points(t(M1), pch = 19, col = "blue")
points(t(M), pch = 19, col = "red")
segments(M1[1,], M1[2,], M[1,], M[2,])
dev.off()

## Single point perturbation
DTp <- DT[c(20, 115, 142, 230), .(x0 = x, y0 = y)]
set.seed(42L)
DTp <- DTp[rep(seq(.N), each = 100)
           ][, c("xErr", "yErr") := replicate(2, rnorm(.N) * min(h*0.166), simplify = FALSE)][]
DTp[, x := x0 + xErr][, y := y0 + yErr][]
M1 <- t(as.matrix(DTp[, .(x, y)]))
system.time(l <- scms(M1, Xs, d = 1, h = h * 0.4)) #(, s)
M <- l$M[-nrow(l$M), ]

svglite::svglite("twocircles-mssm-noise-contraction.svg", width = 24, height = 12)
rbind(DTp[, .(x, y)], DT, list(x = M[1,], y = M[2,]))[, plot(x, y, asp = 1, type = 'n')]
DT[, points(x, y, asp = 1, pch = 19, col = gray(0, 0.3))]
segments(M1[1,], M1[2,], M[1,], M[2,], col = "gray75")
DTp[, points(x, y, pch = 19, col = adjustcolor("red", alpha.f = 0.3))]
points(t(M), pch = 19, col = adjustcolor("blue", alpha.f = 0.3))
DTp[, points(x0, y0, pch = 19, col = "black")]
dev.off()
