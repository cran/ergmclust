\name{grad_HMM_stat_undir}
\alias{grad_HMM_stat_undir}
\title{
C++ function for updating the gradient function in a static binary undirected network.
}
\description{
The C++ function grad_HMM_stat_undir() updates the gradient function when the number of clusters is greater than or equal to 2 in a static binary undirected network.
}
\usage{
grad_HMM_stat_undir(theta, gamma, network, N, K)
}

\arguments{
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
}
  \item{gamma}{
  Vector of current \code{gamma} (variational parameters) values.
}
  \item{network}{
  Array of a network adjacency matrices.
}
  \item{N}{
  Number of nodes.
}
  \item{K}{
  Number of clusters.
}
}
\value{
Returns the gradient vector of ELBO with respect to theta parameters for Newton-Raphson update in the M-step.
}
\author{
Authors: Amal Agarwal, Lingzhou Xue

Maintainer: Yinqi Zhang <ybz5148@psu.edu>
}
\keyword{internal}