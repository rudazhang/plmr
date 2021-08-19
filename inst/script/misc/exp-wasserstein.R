## Optimal transport and Wasserstein distances.
## ------------------------------------------------------------------------------------------
## Empirical complexities of wasserstein():
## "wpp", networkflow, ~2,  *p = 1: N = ~530, 0.09s, IQR (0.086, 0.093); (1024, ~530), 0.183;
##                          *p = 2: N = ~530, 0.08s, IQR (0.076, 0.083); (1024, ~530), 0.158;
## "pp", networkflow        *p = 1: N = 1024, 0.09s;
##                          *p = 2: N = 1024, 0.08s;
## "pp", auction, 2.3        p = 1: N = 1024, 1.29s, IQR (1.022, 1.347);
##                           p = 2: N = 1024, 1.76s, IQR (1.355, 2.021);
## "pgrid", networkflow,    *p = 1: 1024, 0.47 (0.40, 0.55);
##                           p = 2: 1024, 0.52 (0.44, 0.60); 4096, 11 (9, 14);
## "pgrid", revised simplex, p = 1: 1024, 0.416s; 4096, 6.95s; 128^2, 710s; 25120, 22min;
##                           2-3.33, 2.25-2.5; relative error (-4.7e-3, 3.6e-3).
## "pgrid", shielding, 2.2, ~p = 2: 1024, 0.25 (0.24, 0.29); 4096, 7 (6.3, 8.0); 128^2, 106s;
## "pgrid", aha,             p = 2: 1024, 8.5s;
## semi-discrete, p = 1: 2.4-2.5;
##                timing: (1024, 30 * 2^3), 16.036s, 90% CI (14.150, 18.605); (1024, 1024), ~734s;
## Integer mass is faster.
## On subwasserstein():
## Estimating Wasserstein distance between measures on a large graph by their empirical distributions
## tends to larger, and the relative error is not very small. (w1 = 2: min -11%, median +20%)
## On partitioning the data space for faster computation of Wasserstein metrics:
## Larger than the exact result, and the error increases with the number of partitions.
## W1 optimal transport plan often involves large movements, while W2 does not.
## ------------------------------------------------------------------------------------------
library(transport)

## Getting started
set.seed(27)
x <- transport::pp(matrix(runif(400), 200, 2))
y <- transport::pp(matrix(runif(400), 200, 2))
result <- transport::transport(x, y, method = "networkflow", threads = 1, fullreturn = TRUE)
par(mai = rep(0.02, 4))
plot(x, y, result)

## wasserstein() vs subwasserstein()
system.time(w <- wasserstein(random64a, random64b)) #4.9s, 4.4s
system.time(sw <- replicate(21, subwasserstein(random64a, random64b, S = 1000))) #4.7s
stats <- bootStatsMean(sw)
stats[, plot(0, avg, ylim = c(0, round(bmax, 1)), yaxs = 'i')]
stats[, segments(0, bmin, 0, bmax)]
points(0, w)

#' Random pgrid generation
#' @param n pixels in each dimension
rpgrid <- function(n = 32L) {
    N <- n^2
    g <- seq.int(n) - 1L
    g <- as.data.table(expand.grid(x = g, y = g))
    alpha <- sample(10, 4, replace = TRUE)
    X <- data.table(x = rbeta(N * 2^4, alpha[1], alpha[2]),
                    y = rbeta(N * 2^4, alpha[2], alpha[4]))
    X[, x := quantize(x, 1/n)]
    X[, y := quantize(y, 1/n)]
    X <- X[, .N, keyby = .(x, y)]
    Xg <- X[g, on = .(x, y)][is.na(N), N := 0L][]
    pgrid(matrix(Xg$N, n), boundary = c(0, n-1, 0, n-1))
}
