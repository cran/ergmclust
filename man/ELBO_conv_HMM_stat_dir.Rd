\name{ELBO_conv_HMM_stat_dir}
\alias{ELBO_conv_HMM_stat_dir}
\title{
C++ function for updating the ELBO convergence function in a static binary directed network.
}
\description{
The C++ function ELBO_conv_HMM_stat_dir() updates the ELBO convergence function when the number of clusters is greater than or equal to 2 in a static binary directed network.
}
\usage{
ELBO_conv_HMM_stat_dir(gamma, alpha, theta, network, N, K)
}
\arguments{
  \item{gamma}{
  Vector of current \code{gamma} (variational parameters) values.
}
  \item{alpha}{
  Vector of current \code{alpha} values.
}
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
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
Returns the value of ELBO for calculating relative error from the previous iterate that is subsequently used in the stopping criteria of the algorithm.
}
\keyword{internal}