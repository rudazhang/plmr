## Quadratic function with heteroscedastic error in both variables.
## Figures: data generating process; SCMS dynamical system.
setwd("~/repo/R/plmr")
devtools::load_all()
library(data.table)

## Data
N <- c(train = 300L, validate = 300L)
DT <- data.table(type = rep(factor(names(N)), times = N))
set.seed(42L)
DT <- cbind(DT, dgpParabola(DT[, .N]))

## Plot parameters
cexPoint <- 1.5
lwdRegular <- 3
lwdBold <- 6

## Subplot 1: Data generating process ----------------------------------------
## png("data-parabola.png", width = 800, height = 800)
op <- par(mar = c(2,2.5,0.1,0.1))
DT[type == "train", plot(x, y, asp = 1, cex = cexPoint, pch = 19,
                         col = adjustcolor("blue", alpha.f = 0.3), xaxt = 'n', yaxt = 'n')]
axis(1, cex.axis = 2)
axis(2, cex.axis = 2)
DT[type == "train"][order(x0)][seq(1, .N, length.out = 20), lines(x0, y0, lwd = lwdBold)]
DT[type == "train", segments(x0, y0, x, y, lwd = lwdRegular, col = gray(0, .3))]
par(op)
## dev.off()

## Parzen model: optimal bandwidth by log-likelihood on validation set.
X <- as.matrix(DT[type == "train", .(x, y)])
Y <- as.matrix(DT[type == "validate", .(x, y)])
DTparzen <- data.table(sigma = seq(1.5e-2, 5e-2, by = 1e-3))
DTparzen[, anll := purrr::map_dbl(sigma, ParzenANLL, Y = Y, X = X)]
DTparzen[anll < -1, plot(sigma, anll, type = 'b')]
DTparzen[, abline(h = min(anll), col = "red")]
optimize(ParzenANLL, c(3e-2, 5e-2), Y=Y, X=X) # sigma = 0.04; ANLL = -1.015
## Compare with Silverman's Rule of Thumb for 1d Gaussian distributions.
h <- silverman(nrow(X), ncol(X), apply(X, 2, sd))
cat("Oversmoothing factors:", h / DTparzen[which.min(anll), sigma], "\n")

## SCMS
Y0 <- t(as.matrix(DT[type == "train", .(x, y)]))
DT[type == "train", outX := outlier(X)]
Xs <- t(as.matrix(DT[!(outX), .(x, y)]))
## llog <- scms(Y0, Xs, d = 1, h = 4e-2 * 2.25, minCos = 0.05)
llog <- scms(Y0, Xs, d = 1, h = 4e-2 * 2.5, minCos = 0.05)

## Subplot: Moving Data to Density Ridge ----------------------------------------
## png("project-parabola.png", width = 800, height = 800)
op <- par(mar = c(2,2.5,0.1,0.1))
DT[!(outX), plot(x, y, asp = 1, cex = cexPoint, pch = 19, col = adjustcolor("blue", .3),
                 xaxt = 'n', yaxt = 'n')]
axis(1, cex.axis = 2)
axis(2, cex.axis = 2)
## DT[!(outX)][order(x0), lines(x0, y0, lwd = 3, lty = 3, col = gray(0, 0.5))]
DT[!(outX)][order(x0), lines(x0, y0, lwd = lwdBold)]
## text(0, 1, 2.5, cex = 4)
points(llog$Y[1,], llog$Y[2,], cex = cexPoint, pch = 19, col = adjustcolor("red", .3))
segments(Y0[1,], Y0[2,], llog$Y[1,], llog$Y[2,], lwd = lwdRegular, col = gray(0, .3))
par(op)
## dev.off()

## SCMS trajectory
n <- 140L
xrange <- c(-0.4, 1.25)
a <- seq(xrange[[1]], xrange[[2]], length.out = n)
DTboundary <- data.table(x = c(rep(a, 2), rep(xrange, each = n)),
                         y = c(rep(xrange, each = n), rep(a, 2)))
DTboundary <- unique(DTboundary)
maxStepSize <- 0.03
maxSteps <- 500L
Y0 <- t(as.matrix(DTboundary[, .(x, y)]))
Y <- array(double(1L), c(dim(Y0), maxSteps))
t <- 1L
Y[,,t] <- Y0
maxCos <- 1
system.time(while(maxCos > 0.05 & t < maxSteps) {
    t <- t + 1L
    llog <- scms(Y0, Xs, d = 1, h = 4e-2 * 3, maxStep = maxStepSize, minCos = 0.1, maxIter = 1L)
    maxCos <- max(llog$cosNormal)
    ## dY <- llog$Y - Y0
    ## Y0 <- Y0 + dY * min(maxStepSize / sqrt(colSums(dY^2)), 1)
    Y0 <- llog$Y
    Y[,,t] <- Y0
}) #27s

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
Y <- Y[,,seq_len(t)]
tooShort <- sqrt(colSums((Y[,,t] - Y[,,1])^2)) < 0.02
ridgeLeft <- which(Y[1,,t] %between% (1e-2 * c(-1, 1)))
ridgeTop <- which(Y[2,,t] %between% (xrange[[2]] + 1e-3 * c(-1, 1)))
ymin <- mean(Y[2, ridgeLeft, t])
xmax <- mean(Y[1, ridgeTop, t])
fromAbove <- ((Y[1,,1] == xrange[[1]]) & (Y[2,,1] > ymin) |
              (Y[2,,1] == xrange[[2]]) & (Y[1,,1] < xmax))
YUpper <- Y[, !(tooShort) & fromAbove, ]
YLower <- Y[, !(tooShort) & !(fromAbove), ]
YUpper <- YUpper[, order(YUpper[1,,t]),]
YLower <- YLower[, order(YLower[1,,t]),]
stopifnot(all(diff(YUpper[1,,t]) > 0) & all(diff(YLower[1,,t]) > 0))

## Subplot: SCMS Dynamical System ----------------------------------------
## png("dynamical-system-parabola.png", width = 800, height = 800)
op <- par(mar = c(2,2.5,0.1,0.1))
## DTgrid[, plot(x, y, asp = 1, type = 'n', xaxs = 'i', yaxs = 'i')]
DT[type == "train", plot(x, y, asp = 1, type = 'n', xaxt = 'n', yaxt = 'n',
                         xlim = xrange, ylim = xrange, xaxs = 'i', yaxs = 'i')]
axis(1, cex.axis = 2)
axis(2, cex.axis = 2)
DT[!(outX)][order(x0), lines(x0, y0, lwd = lwdBold, col = gray(0, 0.3))]
DT[type == "train", points(x, y, pch = 19, col = adjustcolor("blue", 0.2), cex = cexPoint)]
## points(llog$Y[1,], llog$Y[2,], pch = 19, col = adjustcolor("red", .3), cex = cexPoint)
## Plot all
## purrr::walk(seq_len(dim(Y)[[2]]), ~lines(Y[1,.,], Y[2,.,], col = gray(0, .3)))
## Plot subset
Yplot <- YUpper[,select(YUpper, 0.04),]
purrr::walk(seq_len(dim(Yplot)[[2]]), ~lines(Yplot[1,.,], Yplot[2,.,], col = "orange"))
Yend <- Yplot[,,t]
Yplot <- YLower[,select(YLower, 0.04),]
purrr::walk(seq_len(dim(Yplot)[[2]]), ~lines(Yplot[1,.,], Yplot[2,.,], col = "orange"))
Yend <- cbind(Yend, Yplot[,,t])
## DTgrid[order(yr)[c(seq(1, .N, by = .N %/% 150), .N)], lines(xr, yr, lwd = 3, col = "red")]
Yend <- Yend[,order(Yend[2,])]
lines(Yend[1,], Yend[2,], col = "red", lwd = lwdRegular)
par(op)
## dev.off()
