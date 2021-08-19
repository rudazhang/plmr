### Bent circle with doubled density on the upper half ----------------------------------------
## Mean relative error
avgRelErr <- 0.05
c <- avgRelErr / E$absZ / sqrt(3 - 1)
DT <- data.table::data.table(s = runif(N))
DT[, x := cos(3 * pi * s)]
DT[, y := sin(3 * pi * s)]
DT[, z := ifelse(y > 0, y, 0)]
DT[y > 0, y := 0]
DT[, x := x + c * rnorm(N)]
DT[, y := y + c * rnorm(N)]
DT[, z := z + c * rnorm(N)]
DTwire <- data.table::data.table(s = seq(0, 300) / 300)
DTwire[, x := cos(2 * pi * s)]
DTwire[, y := sin(2 * pi * s)]
DTwire[, z := ifelse(y > 0, y, 0)]
DTwire[y > 0, y := 0]
