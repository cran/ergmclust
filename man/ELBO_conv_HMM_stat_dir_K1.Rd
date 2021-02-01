\name{ELBO_conv_HMM_stat_dir_K1}
\alias{ELBO_conv_HMM_stat_dir_K1}
\title{
C++ function for updating the ELBO convergence function in a static binary directed network.
}
\description{
The C++ function ELBO_conv_HMM_stat_undir() updates the ELBO convergence function when the number of cluster is equal to 1 in a static binary directed network.
}
\usage{
ELBO_conv_HMM_stat_dir_K1(theta, network, N)
}
\arguments{
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
}
  \item{network}{
  Array of a network adjacency matrices.
}
  \item{N}{
  Number of nodes.
}
}
\value{
Returns the value of ELBO for calculating relative error from the previous iterate that is subsequently used in the stopping criteria of the algorithm.
}
\keyword{internal}