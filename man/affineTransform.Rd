% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{affineTransform}
\alias{affineTransform}
\title{Transform data space}
\usage{
affineTransform(x, y, theta, right = 0, up = 0)
}
\arguments{
\item{x}{original x coordinate}

\item{y}{original y coordinate}

\item{theta}{rotation angle, clockwise, in radius}

\item{right}{offset rightward}

\item{up}{offset upward}
}
\value{
A list: (xt, yt), transformed x and y coordinates.
}
\description{
From a 2D histogram of grid size 0.02 / 0.03 (smoothing error < ~7.7e-3 / 1.15e-2)
}
