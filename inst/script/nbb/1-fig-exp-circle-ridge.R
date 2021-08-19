library(plmr)
library(FNN)
options("digits" = 5L)

## The unit circle with radial Gaussian noise
N <- 2^7
B <- N
sigmaNoise <- 0.2
set.seed(42L)
DTsample <- plmr:::rCircularGauss(N, sigma = sigmaNoise)
DTsample[, I := .I]

#' Ridge estimation
#' @param DT sample data, a data.table.
#' @param oversmooth oversmoothing factor.
#' @param silent wheather warning messages shall be printed.
getRidge <- function(DT, oversmooth = 2, silent = FALSE) {
    stopifnot(all(c("x", "y") %in% names(DT)))
    X <- DT[, .(x, y)]
    hLoocv <- plmr::LOOCVMLBandwidth(X, hlim = c(0.05, 0.5))
    h <- hLoocv$minimum * oversmooth
    DTridge <- plmr::estimateRidge(X, d = 1, sigma = h, maxStep = .01, maxIter = 500, silent = TRUE)
    ## Check if all points are on the ridge.
    DTtail <- DTridge[order(eigengap)[seq(.N * .2)], .(I = .I, eigengap)]
    lmfit <- DTtail[, lm(eigengap ~ I)]
    slope <- coefficients(lmfit)[["I"]]
    cutoff <- DTtail[, which.max(eigengap - 4 * slope * I)]
    if (!silent) message("Non-ridge points: ", cutoff - 1L)
    ## Drop non-ridge points.
    DTridge[order(eigengap), isRidge := ifelse(.I < cutoff, FALSE, TRUE)]
    DT <- cbind(DT, DTridge[, .(rx = r1, ry = r2, isRidge)])
    DT <- DT[(isRidge), -"isRidge"]
    DT[]
}

## NBB CI of estimated ridge
alpha <- 3
DT <- getRidge(DTsample, oversmooth = alpha)
## Rectify the projection vectors of neighbors onto the normal space, 1d case.
## In case the normal space is one-dimensional, rectification simply means
## taking the length while keeping relative direction with the normal vector.
DT[, rr := plmr:::norm2(rx, ry)]
DT[, dr := plmr:::norm2(x - rx, y - ry)]
DT[, proj := ifelse(r > rr, dr, - dr)]
## Bias of estimated density ridge
offsetNBBMode <- plmr:::mode1d(DT$proj, minD = 1e-9)
## Bootstrap confidence interval of bias
set.seed(42L)
offsetNBBModeB <- replicate(B, plmr:::mode1d(sample(DT$proj, replace = TRUE), minD = 1e-6)) #1.3s
c90 <- quantile(offsetNBBModeB, c(0.05, 0.95))
DT[, c("rxp05", "ryp05") := .(rx + (x - rx) * c90[[1]] / proj, ry + (y - ry) * c90[[1]] / proj)]
DT[, c("rxp95", "ryp95") := .(rx + (x - rx) * c90[[2]] / proj, ry + (y - ry) * c90[[2]] / proj)]

## Bootstrap CI of estimated ridge
set.seed(42L)
system.time(DTb <- purrr::map_dfr(seq(B), function(b) {
    DT <- getRidge(DTsample[sample(.N, replace = TRUE),], oversmooth = alpha)
    DT[, B := b]
})) #26s
## Distance from bootstrap ridge points to estimated ridge
Mr <- as.matrix(DT[, .(rx, ry)])
dist <- function(b) {
    Mrb <- as.matrix(unique(DTb[B == b, .(rx, ry)]))
    FNN::knnx.dist(Mr, Mrb, 1)
}
distB <- purrr::map_dfr(seq(B), ~data.table::data.table(B = ., d = as.vector(dist(.))))
data.table::setorder(distB, B, d)
distB[, ecdf := .SD[, .I / .N], keyby = .(B)]
distB[, score := (ecdf - 1.0 * d / max(d)), keyby = .(B)]
distB[, outlier := .SD[, .I > which.max(score)], keyby = .(B), .SDcols = c("score")]
c90b <- distB[!(outlier), quantile(d, .9)]
DT[, c("rxp05b", "ryp05b") := .(rx - (x - rx) * c90b / proj, ry - (y - ry) * c90b / proj)]
DT[, c("rxp95b", "ryp95b") := .(rx + (x - rx) * c90b / proj, ry + (y - ry) * c90b / proj)]

## Save the results
## data.table::fwrite(DT, "1-exp-circle-CI.csv")

## Plots --------------------------------------------------------------------------------
colPoint <- adjustcolor("blue", 0.5) #gray(0, .5)
colRidge <- "red"
lwdRidge <- 2
colCI <- gray(0, 0.4)
major <- c(-1, 0, 1)
cexText <- 2
Ncircle <- 120
DTcircle <- data.table::data.table(theta = seq(Ncircle) / Ncircle * 2 * pi)
DTcircle[, c("x", "y") := .(cos(theta), sin(theta))]

## Plot: NBB confidence band
## svglite::svglite("1-exp-circle-CI-NBB.svg", width = 8, height = 8)
op <- par(mar = c(2, 2, .2, .2), las = 1)
DT[, plot(x, y, type = 'n', asp = 1, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')]
axis(1, labels = FALSE)
axis(1, at = major, tick = FALSE, line = -0.1, cex.axis = cexText)
axis(2, labels = FALSE)
axis(2, at = major, tick = FALSE, line = -0.5, cex.axis = cexText)
DTcircle[c(seq(.N), 1), lines(x, y, lwd = lwdRidge)]
DT[, points(x, y, pch = 19, col = colPoint)]
DT[, points(rx, ry, pch = 19, col = colRidge)]
DT[order(theta)[c(seq(.N), 1, 1, seq(.N, 1))],
       polygon(c(rxp95[seq(1, .N/2)], rxp05[-seq(1, .N/2)]),
               c(ryp95[seq(1, .N/2)], ryp05[-seq(1, .N/2)]), col = colCI, border = NA)]
## Ridge + mode (in normal space)
## DT[, points(rx + (x - rx) * offsetNBBMode / proj,
##             ry + (y - ry) * offsetNBBMode / proj, pch = 19, col = "black")]
par(op)
## dev.off()

## Plot: Bootstrap confidence band
## svglite::svglite("1-exp-circle-CI-bootstrap.svg", width = 8, height = 8)
op <- par(mar = c(2, 2, .2, .2), las = 1)
DT[, plot(x, y, type = 'n', asp = 1, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')]
axis(1, labels = FALSE)
axis(1, at = major, tick = FALSE, line = -0.1, cex.axis = cexText)
axis(2, labels = FALSE)
axis(2, at = major, tick = FALSE, line = -0.5, cex.axis = cexText)
DTcircle[c(seq(.N), 1), lines(x, y, lwd = lwdRidge)]
DT[, points(x, y, pch = 19, col = colPoint)]
DT[, points(rx, ry, pch = 19, col = colRidge)]
DT[order(theta)[c(seq(.N), 1, 1, seq(.N, 1))],
       polygon(c(rxp95b[seq(1, .N/2)], rxp05b[-seq(1, .N/2)]),
               c(ryp95b[seq(1, .N/2)], ryp05b[-seq(1, .N/2)]), col = colCI, border = NA)]
par(op)
## dev.off()
