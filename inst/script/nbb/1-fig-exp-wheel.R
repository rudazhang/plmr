## Figure: DNN regression with NBB vs. smooth bootstrap
library(plmr)

## Plot: training and validation error --------------------------------------------------
libpath <- find.package("plmr")
getdata <- function(filename) data.table::fread(file.path(libpath, "extdata/nbb", filename))
DThistoryRaw <- getdata("history-raw.csv")
DThistoryNBB <- getdata("history-nbb.csv")
DThistorySB <- getdata("history-sb.csv")

DThistoryRaw[, valMaeM20 := data.table::frollapply(val_mae, 20, min, align = "right")]
rollmin <- function(x) {
    m3 <- data.table::frollapply(x, 3, min, align = "right")
    m10 <- data.table::frollapply(x, 10, min, align = "right")
    m20 <- data.table::frollapply(x, 20, min, align = "right")
    m100 <- data.table::frollapply(x, 100, min, align = "right")
    ## pmin(m3, m20, m100, na.rm = TRUE)
    pmin(m3, m10, m20, m100, na.rm = TRUE)
}
DThistoryRaw[, maeRM := rollmin(mae)]
DThistoryNBB[, valMaeRM := rollmin(val_mae)]
DThistoryNBB[, maeRM := rollmin(mae)]
DThistorySB[, valMaeRM := rollmin(val_mae)]
DThistorySB[, maeRM := rollmin(mae)]

cexText <- 2 #1.5
lwdLine <- 2
ltyTrain <- 1
ltyVal <- 3
colRaw <- "red"
colNBB <- "blue"
colSB <- "orange"
xmin <- 3

## svglite::svglite("1-exp-wheel-error.svg", width = 8, height = 8)
op <- par(las = 1, mar = c(3, 2.5, 1.5, 1.5))
DThistoryRaw[, plot(epoch, mae, type = 'n', xaxs = 'i', yaxs = 'i',
                    xlim = c(xmin, max(epoch)), log = 'x', ylim = c(0.1, 0.35),
                    xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")]
xmajor <- 10^c(1:3)
xminor <- rep(xmajor, each = 9) * rep(seq(2, 10) / 10, 3)
axis(1, at = xminor, labels = FALSE)
axis(1, at = xmajor, tick = FALSE, line = -0.2, cex.axis = cexText)
ymajor <- 0.1 * c(1:3)
yminor <- axTicks(2)
axis(2, at = yminor, labels = FALSE)
axis(2, at = ymajor, tick = FALSE, line = -0.5, cex.axis = cexText)
abline(h = ymajor[-1], v = xmajor, col = gray(0, .2))
mtext("Mean Absolute Error", side = 3, at = 3, adj = 0, cex = cexText)
mtext("Epoch", side = 1, line = 2 - 0.2, at = max(xmajor), adj = 1, cex = cexText)
DThistoryRaw[, lines(epoch, valMaeM20, col = colRaw, lwd = lwdLine, lty = ltyVal)]
DThistoryRaw[, lines(epoch, maeRM, col = colRaw, lwd = lwdLine, lty = ltyTrain)]
DThistoryNBB[, lines(epoch, valMaeRM, col = colNBB, lwd = lwdLine, lty = ltyVal)]
DThistoryNBB[, lines(epoch, maeRM, col = colNBB, lwd = lwdLine, lty = ltyTrain)]
DThistorySB[, lines(epoch, valMaeRM, col = colSB, lwd = lwdLine, lty = ltyVal)]
DThistorySB[, lines(epoch, maeRM, col = colSB, lwd = lwdLine, lty = ltyTrain)]
legend("topright", c("train", "validate", "NBB train", "NBB validate", "SB train", "SB validate"),
       col = rep(c(colRaw, colNBB, colSB), each = 2),
       lty = rep(c(ltyTrain, ltyVal), 3),
       lwd = lwdLine, cex = cexText, bty = 'n')
par(op)
## dev.off()

## Plot: augmented points --------------------------------------------------
Ntrain <- 32L
knn <- Ntrain / 2
k <- knn
id <- seq(Ntrain)
l <- 8L
n <- 2 * l + 1
ctx <- c(n, 1)

Mplot <- Mnbbsf[Mnbbsf[,"i"] %in% id & Mnbbsf[,"k"] <= k, ctx]
xrange <- range(Mplot[,1])
colData <- "red"
colNew <- adjustcolor("blue", 0.2)
colSBNew <- gray(0.9)
cexText <- 2 #1.5

## svglite::svglite("1-exp-wheel-augmented.svg", width = 8, height = 8)
op <- par(las = 1, mar = c(3, 3.5, 1.5, 0.5))
## plot(Mplot, type = 'n', xlim = c(0, 2.1), asp = 1, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
plot(Mplot, type = 'n', xlim = xrange + c(-1, 1) * diff(xrange) / 100, xaxs = 'i', asp = 1,
     xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
axis(1, labels = FALSE)
axis(1, tick = FALSE, line = -0.2, cex.axis = cexText)
axis(2, labels = FALSE)
axis(2, tick = FALSE, line = -0.5, cex.axis = cexText)
mtext("theta", side = 1, line = 2 - 0.2, at = xrange[2], adj = 1, cex = cexText)
mtext("x1", side = 3, at = xrange[1], adj = 1, cex = cexText)
points(MtrainSB[-seq(nrow(Mtrain)), c("thetaObs", "xObs_1")], pch = 19, col = colSBNew)
curve(cos(pi * (x + gamma[1])), from = 0, to = 2, add = TRUE)
points(Mplot, pch = 19, col = colNew)
points(Mtrain[, c("thetaObs", "xObs_1")], pch = 19, col = colData)
legend("bottomright", c("original", "NBB augmented", "SB augmented"),
       col = c(colData, colNew, colSBNew),
       pch = 19, cex = cexText, bty = 'n')
par(op)
## dev.off()
