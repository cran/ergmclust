\name{grad_HMM_stat_dir_oe}
\alias{grad_HMM_stat_dir_oe}
\title{
C++ function for updating the gradient function in a static binary directed network.
}
\description{
The C++ function grad_HMM_stat_dir_oe() updates the gradient function when the number of clusters is greater than or equal to 2 in a static binary directed network. (The above description may need further improvement.)
}
\usage{
grad_HMM_stat_dir_oe(theta, gamma, network, N, K)
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
\keyword{internal}