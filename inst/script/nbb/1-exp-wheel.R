## Data augmentation for regression by deep neural network:
## A rotating wheel example: d = 1, n = 17 (discretized Hilbertian data), N = 25.
## Requires TensorFlow and `keras` for deep neural net training.
## TensorFlow tutorials:
##   https://www.tensorflow.org/tutorials/keras/regression
##   https://www.tensorflow.org/tutorials/keras/overfit_and_underfit
##   https://www.tensorflow.org/tutorials/images/data_augmentation
library(plmr)
options("digits" = 5L)

## Data sample ------------------------------------------------------------
l <- 8L
n <- 1 + 2 * l
N <- 1000
## Add noise to discretization of theta to avoid singularity.
noiseGamma <- 1 / N

## No need to normalize range: theta, (0, 2); x, y, (-1, 1).
## Ambient Gaussian noise of equal size.
rangeXYT <- 2
noiseRel <- 0.05
noiseCoord <- rangeXYT * noiseRel
noiseTheta <- rangeXYT * noiseRel

set.seed(42L)
gamma <- (seq(l) - noiseGamma * runif(l)) / l * rangeXYT
theta <- seq(N) / N * rangeXYT
## Although angle is periodic, modulation breaks the data manifold.
## If desirable, extend data by half period on both ends for regression and ridge estimation.
thetaObs <- (theta + noiseTheta * rnorm(N))

## Long table
DT <- data.table::CJ(i = seq(N), j = seq(l))
DT[, thgam := theta[i] + gamma[j]]
DT[, x := cos(pi * thgam)]
DT[, y := sin(pi * thgam)]
DT[, xObs := x + noiseCoord * rnorm(.N)]
DT[, yObs := y + noiseCoord * rnorm(.N)]
DT[, thgamObs := thetaObs[i] + gamma[j]]
DT[, xTrue := cos(pi * thgamObs)]
DT[, yTrue := sin(pi * thgamObs)]
DT[, thgamObs := NULL]

## Wide table
DTwide <- data.table::dcast(DT, i ~ j, value.var = c("xObs", "yObs"))
DTwide[, thetaObs := thetaObs[i]]
DTwide[, i := NULL]
data.table::setorder(DTwide, thetaObs)
Ntrain <- 32L
thetaGrid <- seq(0, 2, length.out = Ntrain + 1)
idTrain <- DTwide[, .(.I, thetaObs)][list(t = thetaGrid), on = .(thetaObs = t), roll = TRUE, I[-1]]
idTrain <- idTrain + 1L #Shift a little
## Shuffle
set.seed(42L)
DTtrain <- DTwide[idTrain,][sample(.N),]
DTval <- DTwide[-idTrain,][sample(.N),]
Mtrain <- as.matrix(DTtrain)
Mval <- as.matrix(DTval)

## Regression views: x-t, y-t
alpha <- 0.6
plot(Mval[,"thetaObs"], Mval[,"xObs_1"], pch = 19, col = gray(0, 0.1), asp = 1)
points(Mtrain[,"thetaObs"], Mtrain[,"xObs_1"], pch = 19, col = adjustcolor("red", alpha))
curve(cos(pi * (x + gamma[1])), from = 0, to = 2, add = TRUE)
DT[i %in% idTrain & j == 1, points(thetaObs[i], xTrue, pch = 19, col = adjustcolor("blue", alpha/2))]
dev.new()
plot(Mval[,"thetaObs"], Mval[,"yObs_1"], pch = 19, col = gray(0, 0.1), asp = 1)
points(Mtrain[,"thetaObs"], Mtrain[,"yObs_1"], pch = 19, col = adjustcolor("red", alpha))
curve(sin(pi * (x + gamma[1])), from = 0, to = 2, add = TRUE)
DT[i %in% idTrain & j == 1, points(thetaObs[i], yTrue, pch = 19, col = adjustcolor("blue", alpha/2))]

## Phase space views: x-y (one, all, noiseless)
plot(Mtrain[,"xObs_1"], Mtrain[,"yObs_1"], asp = 1)
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
dev.new()
DT[i %in% idTrain, plot(xObs, yObs, asp = 1)]
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
dev.new()
DT[i %in% idTrain, plot(x, y, asp = 1)]
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
DT[i %in% idTrain & j == 1, points(x, y, pch = 19)]

## 3D views: (x, y, t)
library(rgl)
DT[i %in% idTrain, rgl::plot3d(xObs, yObs, thetaObs[i])]
rgl::plot3d(Mtrain[, c("xObs_1", "yObs_1", "thetaObs")])

## Data augmentation by NBB ----------------------------------------------------------------------
## Extend by half period on both ends to avoid boundary effect.
DTtrainpre <- DTtrain[thetaObs > 1]
DTtrainpre[, thetaObs := thetaObs - 2]
DTtrainpost <- DTtrain[thetaObs < 1]
DTtrainpost[, thetaObs := thetaObs + 2]
DTtrainext <- rbind(DTtrain, DTtrainpre, DTtrainpost)
Mtrainext <- as.matrix(DTtrainext)

## Bandwidth selection
hVal <- plmr::MLBandwidth(DTtrain, DTval, hlim = c(0.05, 0.5)) # preferred
hLoocv <- plmr::LOOCVMLBandwidth(DTtrainext, hlim = c(0.05, 0.5)) # larger
oversmooth <- 4 #3.5,4
sigma <- hVal$minimum * oversmooth

## Ridge estimation
ret <- plmr::estimateRidge(DTtrainext, 1, sigma, minCos = 0.01,
                     maxStep = 0.01, maxIter=500, returnV = TRUE)
DTret <- ret$X0
colRidge <- paste0("r", seq(n))
DTridge <- DTret[, .SD, .SDcols = colRidge]
Mridge <- as.matrix(DTridge)
V <- ret$V

colRidgeTheta <- paste0("r", 2*l + 1)
colRidgeXY <- paste0("r", seq(2*l))
colTheta <- c(colRidgeTheta, "thetaObs")
Mdiff <- Mtrainext - Mridge

## Plot data--ridge pairs
cxyt <- c(1, 9, 17)
rgl::plot3d(Mtrainext[, cxyt], asp = "iso")
rgl::plot3d(Mridge[, cxyt], col = "red", add = TRUE)
DT[j == 1, rgl::lines3d(x, y, theta[i])]
Mseg <- matrix(t(cbind(Mtrainext[, cxyt], Mridge[, cxyt])), byrow = TRUE, ncol = 3)
rgl::segments3d(Mseg, col = "red")

## Plot densified ridge
Nridge <- 100
R0 <- t(as.matrix(DTval[order(thetaObs)[seq(1, .N, length.out = Nridge)],]))
X0 <- t(as.matrix(DTtrainext))
llog <- scms(R0, X0, d = 1, h = sigma, minCos = 0.01,
             omega = 1, maxStep = NULL, maxIter = 100L,
             returnV = FALSE, silent = FALSE)
Y <- t(llog$Y)
rgl::plot3d(Y[, c("xObs_1", "yObs_1", "thetaObs")], col = "red", add = TRUE)

## Normal-bundle bootstrap: without smooth frame
library(FNN)
knn <- ceiling(Ntrain * 0.5)
ret <- FNN::get.knn(Mridge, knn) # Find neighbors
idNeighbors <- ret$nn.index
Dist <- ret$nn.dist
Mnbb <- cbind(Mridge[rep(seq(nrow(Mridge)), knn), ] + Mdiff[as.vector(idNeighbors), ],
              k = rep(seq(knn), each = nrow(Mridge)),
              i = rep(seq(nrow(Mridge)), knn))

## Plot NBB augmented points
ct <- 17L
ck <- ct + 1L
id <- seq(nrow(Mtrain))
rgl::plot3d(Mnbb[Mnbb[,"i"] %in% id, cxyt], col = "blue", add = TRUE)

## Test out smooth frame construction on tangent frame ##
## Plot tangent vectors of the ridge
M <- Mridge[id, cxyt]
V1 <- V[cxyt, 1, id]
Mtan <- matrix(rbind(t(M), t(M) + V1), byrow = TRUE, ncol = 3)
rgl::segments3d(Mtan, col = "magenta")
## Construct a smooth frame of the tangent bundle
V1 <- V[,1,]
dim(V1) <- c(n, 1, dim(V)[[3]])
E1 <- SmoothFrame(V1, Mridge)
## Plot the smooth frame of the tangent bundle: (x, y, z)
M <- Mridge[id, cxyt]
E1 <- E1[cxyt, 1, id]
MtanSmooth <- matrix(rbind(t(M), t(M) + E1), byrow = TRUE, ncol = 3)
rgl::segments3d(MtanSmooth, col = "green")

## Normal-bundle bootstrap with smooth frame
E <- SmoothFrame(V[,-1,], Mridge) # Smooth frame of the normal bundle
Mcoord <- vapply(seq(nrow(Mdiff)), function(i) t(E[,,i]) %*% Mdiff[i,], double(dim(E)[2]))
Mnormal <- vapply(seq(nrow(Mridge)),
                  function(i) E[,,i] %*% Mcoord[, idNeighbors[i,]],
                  double(ncol(Mridge) * knn))
dim(Mnormal) <- c(ncol(Mridge), knn * nrow(Mridge))
Mnormal <- t(Mnormal)
Mnbbsf <- cbind(Mridge[rep(seq(nrow(Mridge)), each = knn), ] + Mnormal,
                k = rep(seq(knn), nrow(Mridge)),
                i = rep(seq(nrow(Mridge)), each = knn))

## Plot augmented points
rgl::plot3d(Mnbbsf[Mnbbsf[,"i"] %in% id, cxyt], col = "orange", add = TRUE)

## Regression views: x1
## k <- knn / 2
k <- knn
ctx <- c(17, 1)
## plot(Mval[,"thetaObs"], Mval[,"xObs_1"], pch = 19, col = gray(0, 0.1))
plot(Mval[,"thetaObs"], Mval[,"xObs_1"], type = 'n', xlim = c(0, 2.1))
points(Mnbbsf[Mnbbsf[,"i"] %in% id, ctx], pch = 19, col = adjustcolor("blue", 0.2))
points(Mnbb[Mnbb[,"i"] %in% id, ctx], col = "magenta")
points(Mtrain[,"thetaObs"], Mtrain[,"xObs_1"], pch = 19, col = "red")
curve(cos(pi * (x + gamma[1])), from = 0, to = 2, add = TRUE)
## Mnbbsf, Mtrain, gamma

## Regression views: y1
cty <- c(17, 9)
plot(Mval[,"thetaObs"], Mval[,"yObs_1"], pch = 19, col = gray(0, 0.1))
points(Mtrain[,"thetaObs"], Mtrain[,"yObs_1"], pch = 19, col = adjustcolor("red", alpha))
curve(sin(pi * (x + gamma[1])), from = 0, to = 2, add = TRUE)
points(Mnbb[Mnbb[,"i"] %in% id, cty], col = "magenta")
points(Mnbbsf[Mnbbsf[,"i"] %in% id, cty], col = "blue")
## NOTE: For this problem, NBB results do not appear to differ much with or without smooth frame.

## DNN training ----------------------------------------------------------------------
## Trained models are often not globally optimal,
## so larger nets does not guarantee better trained models.
library(keras)
library(tibble)
library(ggplot2)

k <- knn   #16-fold
MtrainNBB <- rbind(Mtrain, Mnbbsf[Mnbbsf[,"i"] %in% id & Mnbbsf[,"k"] <= k, colRidge])

## Comparison: smooth bootstrap (k-fold)
set.seed(42L)
MtrainRep <- Mtrain[rep(seq(nrow(Mtrain)), each = k), ]
## Use h: leave-one-out cross validation maximum likelihood bandwidth
Mnoise <- hLoocv$minimum * rnorm(length(MtrainRep))
## Use validation maximum likelihood bandwidth (~14% smaller than h)
## Mnoise <- hVal$minimum * rnorm(length(MtrainRep))
MtrainSB <- rbind(Mtrain, MtrainRep + Mnoise)
## Much larger scatter than NBB.
points(MtrainSB[,"thetaObs"], MtrainSB[,"xObs_1"], pch = 19, col = "orange")

## Use a sequential model with two densely connected hidden layers.
NewModel <- function() {
    model <- keras_model_sequential() %>%
        layer_dense(units = 256, activation = "relu", input_shape = 1) %>%
        layer_dense(units = 128, activation = "relu") %>%
        layer_dense(units = 64, activation = "relu") %>%
        layer_dense(units = 32, activation = "relu") %>%
        layer_dense(units = 2*l)
    model <- compile(model, loss = "mse", #mean_absolute_error
                     optimizer = optimizer_adam(lr = 0.015, beta_1 = 0.9, beta_2 = 0.9),
                     metrics = list("mse", "mean_absolute_error"))
}
modelRaw <- NewModel()
summary(modelRaw)
modelNBB <- NewModel()
## summary(modelNBB)
modelSB <- NewModel()

## Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
    on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
})

## Train the model and store training stats
## `validation_split = 0.2`: it takes the latter 20%, not randomized.
eps <- 1000
historyRaw <- modelRaw %>% fit(
  Mtrain[,n], Mtrain[,-n],
  epochs = eps,
  validation_data = list(Mval[,n], Mval[,-n]),
  verbose = 0,
  callbacks = list(print_dot_callback)
)
fitRaw <- modelRaw

historyNBB <- modelNBB %>% fit(
  MtrainNBB[,n], MtrainNBB[,-n],
  epochs = eps,
  validation_data = list(Mval[,n], Mval[,-n]),
  verbose = 0,
  callbacks = list(print_dot_callback)
)
fitNBB <- modelNBB
## DThistoryRaw <- data.table::as.data.table(historyRaw$metrics)
## cols <- sub("mean_absolute_error", "mae", names(DThistoryRaw))
## data.table::setnames(DThistoryRaw, cols)
## DThistoryRaw[, epoch := .I]
## DThistoryNBB <- data.table::as.data.table(historyNBB$metrics)
## data.table::setnames(DThistoryNBB, cols)
## DThistoryNBB[, epoch := .I]

historySB <- modelSB %>% fit(
  MtrainSB[,n], MtrainSB[,-n],
  epochs = eps,
  validation_data = list(Mval[,n], Mval[,-n]),
  verbose = 0,
  callbacks = list(print_dot_callback)
)
fitSB <- modelSB

dev.new()
with(historyRaw$metrics, {
    plot(mean_absolute_error, col = adjustcolor("red", alpha), pch = 19,
         ## xlim = c(0, eps), #eps, 500, 100
         xlim = c(5, eps), log = 'x', #eps, 500, 100
         ylim = c(0.1, .35), xaxs = 'i', yaxs = 'i');
    points(val_mean_absolute_error, col = gray(0, alpha), pch = 19)
    abline(h = min(val_mean_absolute_error), col = "red")
})
with(historyNBB$metrics, {
    points(mean_absolute_error, col = adjustcolor("red", alpha), pch = 15,
         xlim = c(0, length(loss)), ylim = c(0, .6), xaxs = 'i', yaxs = 'i');
    points(val_mean_absolute_error, col = adjustcolor("blue", alpha), pch = 15)
    abline(h = min(val_mean_absolute_error), col = "blue")
})

with(historySB$metrics, {
    points(mean_absolute_error, col = adjustcolor("red", alpha), pch = 17,
         xlim = c(0, length(loss)), ylim = c(0, .6), xaxs = 'i', yaxs = 'i');
    points(val_mean_absolute_error, col = adjustcolor("blue", alpha), pch = 17)
    abline(h = min(val_mean_absolute_error), col = "blue")
})
DThistorySB <- as.data.table(historySB$metrics)
fwrite(DThistorySB, "history-sb.csv")

with(historyRaw$metrics, {
    plot(mse, col = adjustcolor("red", alpha), pch = 19,
         xlim = c(0, eps), #eps, 500, 100
         ylim = c(0, .6), xaxs = 'i', yaxs = 'i');
    points(val_mse, col = gray(0, alpha), pch = 19)
    abline(h = min(val_mse), col = "red")
})
with(historyNBB$metrics, {
    points(mse, col = adjustcolor("red", alpha), pch = 15,
         xlim = c(0, length(loss)), ylim = c(0, .6), xaxs = 'i', yaxs = 'i');
    points(val_mse, col = adjustcolor("blue", alpha), pch = 15)
    abline(h = min(val_mse), col = "blue")
})

j <- 2
t <- seq(-0.5, 2.5, by = 0.01)
x <- fitRaw %>% predict(t)
plot(Mtrain[,n], Mtrain[,j], pch = 19, col = "red")
curve(cos(pi * (x + gamma[j])), from = 0, to = 2, add = TRUE) #j=1...8
## curve(sin(pi * (x + gamma[j-l])), from = 0, to = 2, add = TRUE) #j=9...16
## lines(t, x[,j], col = "red")
points(Mtrain[,n], Mtrain[,j])
points(Mval[,n], Mval[,j], pch = 19, col = gray(0, 0.1))
x <- fitNBB %>% predict(t)
lines(t, x[,j], col = "blue")
## Noise smoothed out the peaks.
x <- fitSB %>% predict(t)
lines(t, x[,j], col = "green")
