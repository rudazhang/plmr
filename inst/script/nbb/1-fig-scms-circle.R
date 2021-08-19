setwd("~/repo/R/plmr")
devtools::load_all()
library(data.table)

## The unit circle with radial Gaussian noise
## Figure: SCMS dynamical system.
sigmaNoise <- 0.2
f <- function(x, sigma = sigmaNoise) {
    e1 <- exp(-0.5 * ((x - 1) / sigma)^2)
    e2 <- exp(-0.5 * ((x + 1) / sigma)^2)
    (e1 - e2) / (e1 + e2) - x
}
## Ridge numerically remains at 1 until noise is really large.
rRidge <- uniroot(f, c(1e-3, 1))$root
sigmaN256 <- 0.144
oversmoothing <- 2
## Unit circle
DTcircle <- data.table(theta = seq(180) / 180 * 2 * pi)
DTcircle[, c("x", "y") := .(cos(theta), sin(theta))]
## Sample
N <- 300L
set.seed(42L)
DTsample <- rCircularGauss(N, sigma = sigmaNoise)
sigmaSCMS <- sigmaN256 * oversmoothing * (N / 256L)^(-0.24)
DTridge <- estimateRidge(DTsample[, .(x, y)], d = 1, sigma = sigmaSCMS)

## SCMS trajectory: outside
n <- 140L
xrange <- c(-1.5, 1.5)
a <- seq(xrange[[1]], xrange[[2]], length.out = n)
DTboundary <- data.table(x = c(rep(a, 2), rep(xrange, each = n)),
                         y = c(rep(xrange, each = n), rep(a, 2)))
DTboundary <- unique(DTboundary)
maxStepSize <- 0.05
maxSteps <- 100L
X <-t(as.matrix(DTsample[, .(x, y)]))
Y0 <- t(as.matrix(DTboundary[, .(x, y)]))
Y <- array(double(1L), c(dim(Y0), maxSteps))
t <- 1L
Y[,,t] <- Y0
maxCos <- 1
while(maxCos > 0.05 & t < maxSteps) {
    t <- t + 1L
    llog <- scms(Y0, X, d = 1, h = sigmaSCMS, maxStep = maxStepSize, minCos = 0.05, maxIter = 1L)
    maxCos <- max(llog$cosNormal)
    Y0 <- llog$Y
    Y[,,t] <- Y0
}

#' Select trajectories whose endpoints are no less than minDist.
select <- function(Y, minDist = 0.04) {
    M <- Y[,,t]
    s <- rep(NA, dim(Y)[[2]])
    s[[1]] <- TRUE
    while(is.na(s[[length(s)]])) {
        i <- max(which(s))
        j <- i + 1
        while(sqrt(sum((M[,i] - M[,j])^2)) < minDist) {
            s[[j]] <- FALSE
            j <- j + 1
            if (j > length(s)) break
        }
        if (j > length(s)) break
        s[[j]] <- TRUE
    }
    s
}

## Subset
YOuter <- Y[,,seq_len(t)]
cosRidgeOuter <- YOuter[1,,t]
sinRidgeOuter <- YOuter[2,,t]
acosRidgeOuter <- acos(cosRidgeOuter)
thetaRidgeOuter <- acosRidgeOuter * ifelse(sinRidgeOuter >= 0, 1, -1)
YOuterSorted <- YOuter[, order(thetaRidgeOuter),]
YplotOuter <- YOuterSorted[, select(YOuterSorted, 0.07),]

## SCMS trajectory: inside
radiusInner <- 0.1
n <- 1000
DTinner <- data.table(theta = seq(n) / n * 2 * pi)
DTinner <- DTinner[, .(x = radiusInner * cos(theta),
                       y = radiusInner * sin(theta))]
Y0 <- t(as.matrix(DTinner[, .(x, y)]))
Y <- array(double(1L), c(dim(Y0), maxSteps))
t <- 1L
Y[,,t] <- Y0
maxCos <- 1
maxStepSize <- 0.01
maxSteps <- 100L
system.time(while(maxCos > 0.05 & t < maxSteps) {
    t <- t + 1L
    llog <- scms(Y0, X, d = 1, h = sigmaSCMS, maxStep = maxStepSize, minCos = 0.05, maxIter = 1L)
    maxCos <- max(llog$cosNormal)
    Y0 <- llog$Y
    Y[,,t] <- Y0
}) #10s
YInner <- Y[,,seq_len(t)]
cosRidgeInner <- YInner[1,,t]
sinRidgeInner <- YInner[2,,t]
acosRidgeInner <- acos(cosRidgeInner)
thetaRidgeInner <- acosRidgeInner * ifelse(sinRidgeInner >= 0, 1, -1)
YInnerSorted <- YInner[, order(thetaRidgeInner),]
YplotInner <- YInnerSorted[, select(YInnerSorted, 0.07),]

## SCMS trajectory: insider reverse
Y0 <- YplotInner[,,1]
Y <- array(double(1L), c(dim(Y0), maxSteps))
t <- 1L
Y[,,t] <- Y0
maxDelta <- maxStepSize
minStepSize <- 1e-4
maxSteps <- 100L
while(maxDelta > minStepSize & t < maxSteps) {
    t <- t + 1L
    llog <- scms(Y0, X, d = 1, h = sigmaSCMS, maxStep = maxStepSize, minCos = 0.05, maxIter = 1L)
    dy <- llog$Y - Y0
    maxDelta <- max(norm2(dy[1,], dy[2,]))
    Y0 <- Y0 - dy #Reverse steps
    Y[,,t] <- Y0
}
YplotRev <- Y

## Plot parameters
cexPoint <- 1.5
lwdRegular <- 3
lwdBold <- 6

## Plot
## png("dynamical-system-circle.png", width = 800, height = 800)
op  <- par(mar = c(2,2.5,0.1,0.1))
DTcircle[, plot(x, y, type = 'n', asp = 1, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '',
                xlim = xrange, ylim = xrange, xaxs = 'i', yaxs = 'i')]
axis(1, at = seq(-1, 1), cex.axis = 2)
axis(2, at = seq(-1, 1), cex.axis = 2)
DTcircle[c(seq(.N), 1), lines(x, y, lwd = lwdBold, col = gray(0, 0.2))]
DTsample[, points(x, y, pch = 19, col = adjustcolor("blue", 0.2), cex = cexPoint)]
## DTridge[, points(r1, r2, pch = 19, col = "red")]
purrr::walk(seq_len(dim(YplotOuter)[[2]]), ~lines(YplotOuter[1,.,], YplotOuter[2,.,], col = "orange"))
purrr::walk(seq_len(dim(YplotInner)[[2]]), ~lines(YplotInner[1,.,], YplotInner[2,.,], col = "orange"))
purrr::walk(seq_len(dim(YplotRev)[[2]]), ~lines(YplotRev[1,.,], YplotRev[2,.,], col = "orange"))
nPoints <- dim(YplotOuter)[[2]]
t <- dim(YplotOuter)[[3]]
Yend <- YplotOuter[,c(seq(nPoints), 1),t]
lines(Yend[1,], Yend[2,], col = "red", lwd = lwdRegular)
par(op)
## dev.off()
