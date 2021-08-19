## Earthquake Location Distribution along Plates/Faults.
library(data.table)
library(rgl)
library(sf)
library(httr)
library(purrr)

#' Compose an API request URL for the ANSS ComCat comprehensive earthquake catalog.
#' @param parameters A character vector of "key=value".
makeURL <- function(parameters, method = "query") {
    baseURL <- "https://earthquake.usgs.gov/fdsnws/event/1/"
    defaultParameters <- c("format=csv", "eventtype=earthquake", "minmagnitude=4.5")
    paramURL <- paste0(c(defaultParameters, parameters), collapse = "&")
    paste0(baseURL, method, '?', paramURL)
}

DTquery <- data.table::data.table(year = 1900:2019)
DTquery[, starttime := paste0(year, "-01-01")]
DTquery[, endtime := c(starttime[-1], "2019-11-01")]
DTquery[, query := lapply(paste0("starttime=", starttime, '&', "endtime=", endtime), makeURL)]
system.time(response <- purrr::map2(DTquery$query, paste0(DTquery$year, ".csv"),
                        ~httr::GET(.x, httr::write_disk(.y, overwrite = TRUE)))) #~60s / 20 yrs

## After some postprocessing to merge the files.
DT <- data.table::fread("earthquakes-all.csv")
str(DT)
if (all(DT$type == "earthquake")) DT[, type := NULL]
DT[, time := lubridate::as_datetime(time, tz = "UTC")]
DT[, updated := lubridate::as_datetime(updated, tz = "UTC")]

## Data availability and quality by year: all good since 1973.
DT[, .(.N, autoN = sum(status == "automatic")), keyby = .(year = lubridate::year(time))
   ][, {barplot(N, names.arg=year); barplot(autoN, names.arg=year, col = "red", add = TRUE)}]
## Quality subset
DT <- DT[lubridate::year(time) >= 1973 & status == "reviewed"]
DT[, status := NULL]
if(DT[, sum(is.na(depth)) < 5]) DT <- DT[!is.na(depth)]

## Columns after cleaning
adminCols <- c("id", "net", "locationSource", "nst", "magSource", "magNst")
timeCols <- c("time", "updated")
locationCols <- c("latitude", "longitude", "depth")
locErrorCols <- c("horizontalError", "depthError", "gap", "dmin", "rms")
magnitudeCols <- c("mag", "magError", "magType")
humanCols <- c("place")
allCols <- c(adminCols, timeCols, locationCols, locErrorCols, magnitudeCols, humanCols)
data.table::setcolorder(DT, allCols)
data.table::setkey(DT, "time")

## Tapered Gutenberg-Richter distribution of earthquake magnitude: drop around M7.8
## https://nctr.pmel.noaa.gov/education/science/docs/tsun2975/tsun2975_appendixF.pdf
DT[, .N, keyby = .(mag)
   ][order(mag, decreasing=TRUE), pTail := cumsum(N) /sum(N)
     ][, {plot(mag, log(pTail), pch = 19); abline(v = 7.8)}]
DT[, mean(mag >= 7.8)] #3.5e-4

## For earthquakes with both horizontal and vertical error estimates,
## horizontal error is typically 8 km, with 90% within [3.6, 14];
## vertical error is typically 6 km, but 40% take value from {1.7, 1.8, 1.9, 2.0},
## excluding this group, 90% are within [3, 10.5].
DT[!is.na(horizontalError) & !is.na(depthError), {
    hist(horizontalError, breaks = seq(0, 100), col = "gray");
    hist(depthError, breaks = seq(0, 100), col = adjustcolor("red", alpha = 0.5), add = TRUE);
    abline(v = 7.5, lty = 2)}]

## Save data subset in use.
data.table::fwrite(DT, "earthquake.csv")
#DT <- data.table::fread("earthquake.csv")

## Depth distribution and defaults at 10/33/35 km.
DT[, plot.ecdf(depth, verticals=TRUE)]
DT[, {plot.ecdf(depth, verticals=TRUE, xlim = c(-3, 50), ylim = c(0, 0.7), xaxs = 'i', yaxs = 'i');
    abline(v = c(10, 30, 33, 35), col = "red")}]
## Color code by depth
DT[depth <= 10, color := "orange"]
DT[depth > 10 & depth <= 35, color := "blue"]
DT[depth > 35, color := "gray50"]

## WGS84 geographic coordinates to geocentric coordinates
DTtemp <- DT[, .(latitude, longitude, elevation = -depth * 1000)]
df <- sf::st_as_sf(DTtemp, coords=c("longitude", "latitude", "elevation"), crs = 4326)
df <- sf::st_transform(df, crs = "+proj=geocent +units=km")
DT[, c("x", "y", "z") := as.data.frame(sf::st_coordinates(df))]
## Robinson projection, as in "Seismicity of the Earth" map.
## May add "+ellps=WGS84", but does not differ.
DTtemp <- DT[, .(latitude, longitude)]
df <- sf::st_as_sf(DTtemp, coords=c("longitude", "latitude"), crs = 4326)
df <- sf::st_transform(df, crs = "+proj=robin +lon_0=180")
DT[, c("robin1", "robin2") := as.data.frame(sf::st_coordinates(df))]

## 2D Plot: WGS64
DT[, plot((longitude + 360) %% 360, latitude, type = 'n', xaxs = 'i', yaxs = 'i')]
DT[color == "gray50", points((longitude + 360) %% 360, latitude, col = color)]
DT[color == "blue", points((longitude + 360) %% 360, latitude, col = color)]
DT[color == "orange", points((longitude + 360) %% 360, latitude, col = color)]
#abline(h = c(-50, -5), v = c(160, 195), lty = 3)
## 2D Plot: Robinson
DT[, plot(robin1, robin2, type = 'n',
          xlim = 17005833 * c(-1, 1), ylim = 8625155 * c(-1, 1), xaxs = 'i', yaxs = 'i')]
DT[color == "gray50", points(robin1, robin2, col = color)]
DT[color == "blue", points(robin1, robin2, col = color)]
DT[color == "orange", points(robin1, robin2, col = color)]

## 3D plot
rgl::bg3d("black")
DT[, rgl::plot3d(x, y, z, size = 2, col = color, box = FALSE, axes = FALSE,
                 xlab = '', ylab = '', zlab = '')]
rgl::spheres3d(0, 0, 0, radius = 100, color = "white")
rgl::spheres3d(0, 0, 6356, radius = 100, color = "red")

## 3D plot: Fiji basin
DTfiji <- DT[(longitude > 160 | longitude < -165) & latitude > -50 & latitude < -5
             ][!(depth %in% c(10, 33, 35))]
DTfiji[, rgl::plot3d(x, y, z, size = 2, col = "gray", box = FALSE, axes = FALSE,
                     xlab = '', ylab = '', zlab = '')]

## Mean Shift Sampling on Manifolds --------------------------------------------------
## Density Ridge
set.seed(42L)
DTsubset <- DT[color == "orange"][sample(seq(.N), 2e4L), .(x, y, z, robin1, robin2)]
X <- t(as.matrix(DTsubset[, .(x, y, z)]))
M0 <- X
system.time(outliers <- outlier(X)) #60s
DTsubset[, outlier := ifelse(.I %in% outliers, TRUE, FALSE)]
Xs <- t(as.matrix(DTsubset[!(outlier), .(x, y, z)]))

## 2D Plot: Robinson
op <- par(bg = "gray90")
DTsubset[, plot(robin1, robin2, xlim = 17005833 * c(-1, 1), ylim = 8625155 * c(-1, 1),
                xaxs = 'i', yaxs = 'i')]
par(op)

## 3D plot
rgl::bg3d("black")
DTsubset[, rgl::plot3d(x, y, z, size = 4, col = "white", box = FALSE, axes = FALSE,
                       xlab = '', ylab = '', zlab = '')]
rgl::spheres3d(0, 0, 0, radius = 100, color = "white")
rgl::spheres3d(0, 0, 6356, radius = 100, color = "red")

## Convergence criterion
system.time(l <- scms(M0, Xs, d = 1, h = rep(80, 3), maxIter = 5L)) #(19906, 111s)
system.time(l <- scms(M0, Xs, d = 1, h = rep(100, 3), maxIter = 5L)) #(, s)
system.time(l <- scms(M0, Xs, d = 1, h = rep(8, 3), epsilon = 0.055)) #(, s)
system.time(l <- scms(M0, Xs, d = 1, h = rep(8, 3), epsilon = 0.45)) #(, s)

DTsubset <- cbind(DTsubset, data.table::as.data.table(t(l$M[-nrow(l$M),])))
data.table::setnames(DTsubset, paste0('V', seq(3)), paste0(c('x', 'y', 'z'), 'r'))
DTtemp <- DTsubset[, .(xr, yr, zr)]
df <- sf::st_as_sf(DTtemp, coords=c("xr", "yr", "zr"), crs = "+proj=geocent +units=km")
df <- sf::st_transform(df, crs = "+proj=robin +lon_0=180")
DTsubset[, c("robin1r", "robin2r", "robin3r") := as.data.frame(sf::st_coordinates(df))]

DTsubset[, points(robin1r, robin2r, pch = 19, cex = .2, col = "orange")]

rgl::points3d(M[1,], M[2,], M[3,], size = 4, col = "orange")


