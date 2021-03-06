% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MParzen.R
\name{MParzenSample}
\alias{MParzenSample}
\title{Sample from MParzen model}
\usage{
MParzenSample(N, X, model, sigma)
}
\arguments{
\item{N}{Number of samples.}

\item{X}{matrix of training data.}

\item{model}{Top d-th order local covariance structure; output of `MParzenTrain()`.}

\item{sigma}{Maximum noise level added to the normal space of each kernel,
i.e. marginal standard deviation of an isotropic Gaussian. Can be zero.}
}
\value{
A sample of size N from the Manifold Parzen model.
}
\description{
Sample from MParzen model
}
