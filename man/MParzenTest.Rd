% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MParzen.R
\name{MParzenTest}
\alias{MParzenTest}
\title{Testing MParzen model}
\usage{
MParzenTest(y, X, model, sigma)
}
\arguments{
\item{y}{density evaluation point.}

\item{X}{matrix of training data.}

\item{model}{Top d-th order local covariance structure; output of `MParzenTrain()`.}

\item{sigma}{Maximum noise level added to the normal space of each kernel,
i.e. marginal standard deviation of an isotropic Gaussian. Should be greater than zero.}
}
\value{
Manifold Parzen estimator at the evaluation point.
}
\description{
Modified to add noise level \eqn{\min(\sigma, s_d)} to the normal directions,
instead of adding \eqn{\sigma^2} to all directions.
}
\references{
@Vincent2002
}
